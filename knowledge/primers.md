# Primer register — concepts and tools already explained in a notebook

Appended by the knowledge-keeper after every chapter from the notebook's `metadata.fluidpy.primers` and the design's
prerequisite ledger. Later chapters do not repeat a primer: they write a one-sentence reminder ("primed in Ch. 1,
P44") and point here. IDs are the notebook's `P` numbers (not in numeric order inside ch01). P numbers continue across
chapters (ch01 P01–P61, ch02 P62–P86, ch03 P87–P110, ch04 P111–P133, ch05 P134–P148, ch06 P149–P164, ch07
P165–P184, ch08 P185–P199, ch09 P200–P220 + P218a, ch10 P221–P254, ch11 P255–P279, **ch12 P280–P306**); a new chapter starts at the next free number (**ch13: P307**; total after ch12: 307 = P01–P306 and P218a; in ch12 P286 sits in the conventions cell before P280, P283 before P281, P300 before P299)
(ch12 started at P280). Total after ch11: 280 (P01–P279 and P218a — an inserted id: the ch09 lesson review moved integration by parts
before D06, and the later P218 now opens with a recap line pointing to it; ch10's P254 Taylor–Green was added by the lesson
review and sits before P237 in the notebook). In ch11 P269 sits before P268 and P263 before D10 (the lesson review moved both in front of their
first use).

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

## ch02 primers (P62–P86, notebook `ch02_cartesian_tensors`; the term is the exact `nb.primer` title)

| Term (maths / physics / Python) | Explained in (chapter · notebook section · CORE block) | One-line gist (our words) |
|---|---|---|
| **Python and tooling** | | |
| np.einsum index strings (P62) | ch02 · §2.1 · C01 | the summation convention in code: `'i,i'` sums over i; `'ik,kj->ij'` sums k and keeps i, j; the letters after `->` are the free indices |
| matrix multiplication, transpose and identity (P63) | ch02 · §2.1 · C01 | (AB)_ij = Σ_k A_ik B_kj (row of A times column of B); `A @ B`, `A.T` swaps rows and columns, `np.eye(3)` is the identity |
| plotly 3-D arrows, lines and meshes (P64) | ch02 · §2.1 · C01 | `go.Cone` for arrows, `go.Scatter3d(mode='lines')` for edges, `go.Mesh3d` for faces; all stay rotatable on the page |
| np.linalg.norm and np.linalg.qr (P67) | ch02 · §2.2 · C02 | `norm` is the length √(v_i v_i); the Q factor of a random Gaussian matrix, with det fixed to +1, is a random rotation |
| np.arctan2 (P70) | ch02 · §2.6 · C05 | four-quadrant angle of (x, y) in (−π, π]; `arctan(y/x)` cannot tell (−1, −1) from (1, 1) |
| numpy arrays with three axes and np.transpose (P73) | ch02 · §2.7 · C08 | `eps[i, j, k]` indexes a 3×3×3 array; `np.transpose(a, axes)` puts old axis `axes[p]` in position p, so `np.transpose(eps, (1, 2, 0))[i, j, k] = eps[k, i, j]` (demo on `np.arange(27).reshape(3, 3, 3)`) |
| np.meshgrid and the project grid layout (P76) | ch02 · §2.9 · C09 | coordinate arrays for every grid point; **project layout**: arrays `[k, j, i]` = (z, y, x), x on the last axis, components on axis 0, but `Grid.h` and components ordered (x, y, z) |
| numpy broadcasting (P77) | ch02 · §2.9 · C09 | an operation between a grid array and a scalar or smaller array acts at every point at once; neighbour slices `phi[:, 2:] - phi[:, :-2]` build stencils without loops |
| plt.contour, plt.quiver and plt.streamplot (P78) | ch02 · §2.9 · C09 | level curves; an arrow per point; curves tangent to a vector field (1-D x, y and fields on a `[j, i]` grid) |
| scipy.linalg.expm (P79) | ch02 · §2.10 · C12 | the matrix exponential e^{Gt} = I + Gt + (Gt)²/2 + …; the trajectory of u = G·x is x(t) = e^{Gt}x₀ |
| **Physics vocabulary** | | |
| vector area of a closed surface (P69) | ch02 · §2.6 · C05 | Σ n dA over a closed surface is zero; for the tetrahedron this gives the face areas dA_i = n_i dA |
| right-hand rule and orientation (P74) | ch02 · §2.7 · C08 | fingers along the first vector, curl toward the second, thumb = cross product; e₁ × e₂ = e₃ fixes a right-handed frame and orients loops (counterclockwise about n) |
| **Mathematics** | | |
| orthonormal basis, projection and completeness (P65) | ch02 · §2.2 · C02 | e_i·e_j = δ_ij; a component is a dot product (projection); completeness Σ_j e'_j e'_jᵀ = I says a vector is the sum of its projections |
| cosines of angles between unit vectors (P66) | ch02 · §2.2 · C02 | for unit vectors the dot product is the cosine; cos(π/2 − θ) = sin θ, cos(θ + π/2) = −sin θ; numpy works in radians (`np.deg2rad`) |
| limits and orders of smallness (P68) | ch02 · §2.6 · C05 | as an element of size h shrinks, h³ terms vanish faster than h² terms (ratio h → 0), so volume terms drop out of a face balance (ch01 only glossed this) |
| Vieta's formulas (P71) | ch02 · §2.5 · C07 | a cubic with roots λ^k is λ³ − (Σλ)λ² + (Σ_{k<l} λ^kλ^l)λ − Πλ: coefficients are sums and products of roots |
| permutations, cyclic order and parity (P72) | ch02 · §2.7 · C08 | 123, 231, 312 are even (cyclic), 132, 213, 321 odd; ε is +1 on even, −1 on odd, 0 with a repeat (`itertools.permutations`) |
| level sets and the directional derivative (P75) | ch02 · §2.9 · C09 | a level set is where φ is constant; ∂φ/∂n = ∇φ·n is the rate per metre along unit n (chain rule, ch01 P49) |
| eigenvalues and eigenvectors (P80) | ch02 · §2.11 · C13 | A·b = λb: A only stretches b; nonzero b exist when det(A − λI) = 0; `np.linalg.eigh` for symmetric A, `np.roots` for the polynomial |
| complex conjugate (P81) | ch02 · §2.11 · C13 | z̄ = a − ib, z z̄ = \|z\|² ≥ 0, z real ⇔ z = z̄; the tool that proves a symmetric tensor's eigenvalues are real |
| quadratic form and the Rayleigh quotient (P82) | ch02 · §2.11 · C13 | n·τ·n for unit n is the normal stress on the plane ⊥ n; its values lie between the smallest and largest eigenvalue (Gram–Schmidt glossed for repeated λ) |
| volume and surface integrals as midpoint sums (P83) | ch02 · §2.12 · C14 | add f × (cell volume or area) at cell centres, error ∝ h²; iterated integrals do one axis at a time; exact for linear f |
| fundamental theorem of calculus (P84) | ch02 · §2.12 · C14 | ∫_a^b f′ dx = f(b) − f(a): Gauss' theorem in one dimension |
| mean-value theorem for integrals (P85) | ch02 · §2.12 · C15 | ∭_V f dV = f(x*) V for some x* in V; as V shrinks to x₀, x* → x₀ (how the integral definitions become limits) |
| line integral of a vector field around a loop (P86) | ch02 · §2.13 · C16 | parametrise the loop, add u·t ds; circle x = c + R(cos s, sin s), rectangle side by side; extends ch01 P35 (∫p dv) to vectors |
| **ch03 — Python and numerics** | | |
| scipy.integrate.quad and dblquad (P87) | ch03 · §3.1 · N03 (→ C01) | adaptive integration: `quad(f, a, b)` returns the integral and an error estimate; `dblquad` does area integrals (inner limits as functions); used for section averages and Ex. 3.2 |
| solve_ivp options: t_eval, dense_output, events, backward integration (P94) | ch03 · §3.3 · C03 | extends P31: ask for chosen output times, get a continuous solution, stop at an event (a stagnation point), integrate backwards by giving t_span in decreasing order |
| RK4 by hand (P95) | ch03 · §3.3 · C04 | four slope samples per step weighted 1-2-2-1; error ∝ Δt⁴; the from-scratch path line of Ex. 3.1 agrees with `pathline` (extends P30's Euler) |
| np.expm1 and cancellation near zero (P107) | ch03 · §3.5 · C14 | 1 − e^{−x} for tiny x subtracts nearly equal numbers and loses digits; `-np.expm1(-x)` keeps them (Gaussian vortex near r = 0) |
| scipy.optimize.brentq (P108) | ch03 · §3.5 · C14 | root of f(x) = 0 inside a bracket [a, b] where f changes sign; guaranteed and fast; pick the bracket from the physics and exclude trivial roots (x = 0 in 1 + 2x = eˣ) |
| **ch03 — maths** | | |
| cylindrical and spherical unit vectors (P88) | ch03 · §3.1 · N05 (→ C01) | at each point the curvilinear coordinates carry their own right-handed orthonormal unit vectors, which turn as the point moves; components are projections E[k]·u |
| functions of time with parameters (P89) | ch03 · §3.2 · C01 | x = Xe^{αt} is a function of t once the label X is fixed; X picks which particle; differentiate in t with X held fixed |
| inverse functions and sympy solve (P90) | ch03 · §3.2 · C01 | to go from Lagrangian to Eulerian, answer "which particle is at x now?": solve x = Xe^{αt} for X = xe^{−αt} (`sp.solve`) and substitute |
| multivariable chain rule along a path (P91) | ch03 · §3.2 · C02 | f(t) = F(x(t), y(t), z(t), t) ⇒ df/dt = F_x ẋ + F_y ẏ + F_z ż + F_t (extends P49 to a trajectory with explicit t) |
| parametric curves, tangent vector and arc length (P92) | ch03 · §3.3 · C03 | x(s) with tangent dx/ds; with arc length s the tangent is a unit vector, so a streamline solves dx/ds = u/\|u\| |
| parallel vectors and the cross-product test (P93) | ch03 · §3.3 · C03 | a ∥ b ⇔ a = λb ⇔ a × b = 0; the cross-product form survives zero components, the ratio form dx/u = dy/v does not |
| chain rule with a moving frame (P97) | ch03 · §3.3 · C05 | if g(x, t) = f(x − Ut, t) then ∂g/∂t at fixed x = ∂f/∂t − U ∂f/∂x: "∂/∂t" depends on what is held fixed (the trap of (3.9)) |
| multivariable first-order Taylor expansion (P98) | ch03 · §3.4 · C06 | u_i(x + dx) ≈ u_i(x) + (∂u_i/∂x_j) dx_j with an O(\|dx\|²) remainder (extends P26 to several variables) |
| small-angle approximation (P100) | ch03 · §3.4 · C08 | for ε in radians, tan ε ≈ ε and cos ε ≈ 1 with errors of order ε³ and ε²; why stretching changes angles only at second order |
| linear map of a circle is an ellipse (P104) | ch03 · §3.4 · C12 | a matrix M maps the unit circle to an ellipse; for symmetric M its axes are the eigenvectors and its semi-axes the eigenvalues; SVD glossed for finite time |
| polar coordinates as a moving basis (P105) | ch03 · §3.5 · C13 | e_r = (cos θ, sin θ), e_θ = (−sin θ, cos θ) change with θ; area element r dr dθ; line element on a circle r dθ e_θ |
| substitution in an integral (P106) | ch03 · §3.5 · C14 | replace r by s = r²/σ² and convert dr too (ds = 2r dr/σ²); used for the Gaussian vortex's circulation Γ(r) |
| differentiation under the integral sign (P109) | ch03 · §3.6 · C15 | with fixed limits, d/dt ∫_a^b F dx = ∫_a^b ∂F/∂t dx; moving limits add the Leibniz end terms |
| **ch03 — physics vocabulary** | | |
| frames of reference and relative velocity (P96) | ch03 · §3.3 · C05 | a frame = axes + a clock; two frames with parallel axes separating at constant U see velocities that differ by U (Galilean) |
| material line element (P99) | ch03 · §3.4 · C07 | δx joins two nearby fluid particles and is carried with them, so D(δx)/Dt = δu; it stretches and turns with the flow |
| rigid-body velocity Ω × x (P101) | ch03 · §3.4 · C08 | a body turning at Ω about an axis through the origin moves each point at Ω × x, perpendicular to both; its curl is 2Ω |
| angular velocity of a line (P102) | ch03 · §3.4 · C10 | a segment at angle θ from +x turns at θ̇; counterclockwise positive; for a flow, θ̇ = e_θ·G·e(θ) |
| rotating frame of reference (P103) | ch03 · §3.4 · C10 | an observer on a turntable (the Earth) at Ω sees u′ = u − Ω × x at the instant the frames coincide; vorticity drops by 2Ω (first look; Ch. 4 §4.7 develops it) |
| signed swept volume of a moving surface (P110) | ch03 · §3.6 · C15 | in Δt a patch dA moving at b sweeps a prism of volume (b·n Δt) dA, positive when advancing along the outward n, negative when retreating, zero when sliding |
| **ch04 — Python and numerics** | | |
| dataclasses and named results (P111) | ch04 · §4.2 · C01 | `@dataclass(frozen=True)` holds named fields; every fluidpy budget returns one (`MassBudget`, `MomentumBudget`, `EnergyBudget`, `Scales`), so you write `budget.residual`, not "the fourth number" |
| sympy expand, series, removeO, collect and subs (P117) | ch04 · §4.4 · C05 | the approximation toolkit: multiply out, keep powers of ds up to ds¹ (`series(…, ds, 0, 2).removeO()`), group by powers, substitute; how D06 drops the (ds)² terms of the stream-tube element |
| complementary error function erfc (P123) | ch04 · §4.6 · C08 | erfc η = 1 − erf η falls from 1 at η = 0 to 0 (0.48 at 0.5, 0.16 at 1, 0.005 at 2): the shape of anything diffusing in from a suddenly changed wall (Stokes' first problem u = U erfc(y/2√(νt))) |
| **ch04 — maths** | | |
| continuity of a function and the small-ball argument (P112) | ch04 · §4.2 · C02 | continuous f with f(x₀) ≠ 0 keeps its sign in a small ball round x₀, so its integral over that ball is not 0; this is the localisation lemma "∫_V f dV = 0 for every V ⇒ f ≡ 0" |
| product rule for a divergence (P113) | ch04 · §4.2 · C02 | ∇·(ρu) = u·∇ρ + ρ∇·u (index form ∂(ρu_i)/∂x_i = u_i∂ρ/∂x_i + ρ∂u_i/∂x_i); drives D03, D07, D20 (was a ch03 gloss in D24) |
| tensor divergence over the first index (P118) | ch04 · §4.4 · C06 | the net surface force per volume is ∂τ_ij/∂x_i because the traction is f_j = n_iτ_ij; differs from ∂τ_ij/∂x_j for a non-symmetric τ (`tensor_divergence(index=0)`, ch02's default is the second index) |
| isotropic fourth-order tensor (P119) | ch04 · §4.5 · C07 | the only fourth-order tensors unchanged by every rotation are λδ_ijδ_mn + μδ_imδ_jn + γδ_inδ_jm (cited, not proved; checked with 50 random rotations) |
| deviatoric (traceless) part of a tensor (P120) | ch04 · §4.5 · C07 | A = ⅓A_mmδ + (A − ⅓A_mmδ): the average-diagonal (isotropic) part plus a traceless remainder; for S, volume change vs shape change |
| Schwarz's theorem (P121) | ch04 · §4.6 · C08 | mixed partials commute for smooth functions: ∂_i∂_j = ∂_j∂_i; lets ∂/∂x_i(∂u_i/∂x_j) become ∂/∂x_j(∇·u) (D12) and ∂_t swap with ∇ (D26) |
| curl of a curl identity (P122) | ch04 · §4.6 · C08 | ∇×(∇×u) = ∇(∇·u) − ∇²u, from ε–δ (2.19); with ∇·u = 0 the viscous force μ∇²u = −μ∇×ω |
| derivative of a rotating unit vector (P124) | ch04 · §4.7 · C09 | a unit vector fixed to a frame turning at Ω changes only direction: de′/dt = Ω × e′ (tip on a circle round the axis); the whole difference from a Galilean frame |
| product rule for a cross product (P125) | ch04 · §4.7 · C09 | d(a × b)/dt = ȧ × b + a × ḃ, keeping the order of the factors; gives the second Ω × u′ in D15 |
| chain rule for the kinetic energy (P128) | ch04 · §4.8 · C10 | D(½u_j²)/Dt = u_j Du_j/Dt: the kinetic-energy rate per mass is velocity · acceleration (D21) |
| completing the square for tensors (P129) | ch04 · §4.8 · C10 | A_ijA_ij = Σ A_ij² ≥ 0, zero only if every component is; rewrite σ:S as 2μ(dev S):(dev S) + μ_v S_mm² to show ε ≥ 0 (D23; needs δ_ijδ_ij = 3) |
| order-of-magnitude scaling (P130) | ch04 · §4.9 · C13 | replace each quantity by its typical size: ∂u/∂x ~ U/L, ∂²u/∂x² ~ U/L², ∂T/∂x ~ δT/L; ratios of terms decide what can be dropped (Boussinesq, similarity) |
| moving level set and its normal speed (P132) | ch04 · §4.10 · C14 | a surface η(x, t) = 0 has normal n = ∇η/\|∇η\| (toward increasing η) and moves along n at −(∂η/∂t)/\|∇η\|; a point riding on it keeps Dη/Dt = 0 (extends ch02 P75) |
| scaled variables and the chain rule (P133) | ch04 · §4.11 · C15 | t* = Ωt ⇒ ∂/∂t = Ω∂/∂t*; x* = x/l ⇒ ∂/∂x = (1/l)∂/∂x*, ∂²/∂x² = (1/l²)∂²/∂x*²: each term of an equation brings out its scale factor (D30) |
| **ch04 — physics vocabulary** | | |
| momentum flux through a surface (P114) | ch04 · §4.4 · C04 | fluid crossing a patch carries its momentum: mass rate ρ(u − b)·n dA times velocity u — a vector per area per time; outward positive (D05, every CV force) |
| conservative force and its potential (P115) | ch04 · §4.4 · C04 | work independent of route ⇔ g = −∇Φ ⇔ zero work round any loop; gravity Φ = gz with z up (4.18); reused for the centrifugal potential (D18) and the Bernoulli function (D24) |
| torque, moment arm and moment of inertia (P116) | ch04 · §4.4 · C04 (used in C07 D08) | torque M = r × F; angular analogue of mass; a cube of side h has I = ρh⁵/6 about an axis through its centre — why unequal τ₁₂, τ₂₁ spin it up as 1/h² |
| latitude, Earth's rotation rate and the local vertical (P126) | ch04 · §4.7 · C09 | Ω = 2π/86 164 s = 7.292×10⁻⁵ rad/s; at latitude φ its local vertical part is Ω sin φ, so f = 2Ω sin φ (named here, taught in Ch. 13); NH Ω_z > 0 deflects right |
| power of a force and heat flux through a surface (P127) | ch04 · §4.8 · C10 | a force does work at F·u; per area the stress does f·u, per volume gravity ρg·u; heat leaves through dA at q·n dA (outward positive, so the budget has −∮q·n) |
| reduced gravity and buoyancy (P131) | ch04 · §4.9 · C13 | a parcel lighter by Δρ feels g′ = gΔρ/ρ₀ upward; the field version is b = −gρ′/ρ₀; 2 K of warm water (α = 2×10⁻⁴) gives g′ ≈ 0.004 m/s² |
| **ch05 — Python and numerics** | | |
| periodic trapezoid rule on a closed loop (P136) | ch05 · §5.2 · C03 | for a smooth closed loop sampled at N equal parameter steps, the plain average is spectrally accurate (errors cancel round a periodic curve); ∮u·dx with dx/ds from an FFT derivative (`loop_circulation`); glosses: FFT derivative multiplies a wave by 2πik, `np.roll` for the next point |
| advecting many points in one `solve_ivp` call (P137) | ch05 · §5.2 · C03 | one long state vector (all x's, then all y's), a right-hand side that reshapes to (2, N), evaluates the velocity vectorised and flattens back; one adaptive step for the whole loop (`material_loop`) |
| Fourier modes and the FFT Poisson solver (P142) | ch05 · §5.5 · C07 | on a periodic box ∂/∂x → ik_x, ∇² → −\|k\|²; ∇²ψ = −ω becomes ψ̂ = ω̂/\|k\|² per wave (k = 0 set to 0); `np.fft.fft2`/`ifft2` (`velocity_from_vorticity_fft`; → Ch. 10 spectral methods, Ch. 13 PV inversion) |
| Gauss–Legendre quadrature in 3-D and a smoothed kernel (P143) | ch05 · §5.5 · C07 | n nodes per axis integrate polynomials of degree 2n − 1 exactly; a box uses all n³ products of weights (`cylinder_quadrature` does R, φ, z); near a 1/r³ kernel replace r² by r² + ε² so no node divides by zero |
| systems of ODEs for interacting bodies and their invariants (P147) | ch05 · §5.7 · C12 | N vortices = 2N coupled ODEs in one state vector; the exact motion conserves ΣΓx, ΣΓ\|x\|² and Kirchhoff's H = −(1/2π)ΣΓ_jΓ_k ln\|x_j − x_k\|, so watching them stay flat checks the integrator |
| complete elliptic integrals with the parameter m (P148) | ch05 · §5.7 · C13 | K(m) = ∫dθ/√(1 − m sin²θ), E(m) = ∫√(1 − m sin²θ) dθ over [0, π/2]; ⚠️ `scipy.special.ellipk(m)`, `ellipe(m)` take m = k², not the modulus k (a silent 5 % bug); K → ∞ as m → 1 (ring velocities) |
| **ch05 — maths** | | |
| partial integration with an unknown function (P134) | ch05 · §5.1 · C02 | integrating ∂p/∂r in r leaves a "constant" f(z) (∂f/∂r = 0), fixed by the other equation ∂p/∂z (D02, the rotating tank; sympy `dsolve` demo) |
| closed-loop integral of an exact differential (P135) | ch05 · §5.2 · C03 | ∮dF = 0 for a single-valued F; fails for a multi-valued one such as the polar angle round the origin (2π); Kelvin uses it for ½u² and for the pressure function (needs ρ = ρ(p)) |
| loop vector area ½∮x × dx (P138) | ch05 · §5.3 · C05 (reused in C11) | a vector of length = enclosed area (flat loop) along the loop's right-hand normal; for a curved loop ∫n dA over any spanning surface; 2-D = the shoelace formula (`loop_vector_area`) |
| Poisson equation and Green's function (P139) | ch05 · §5.5 · C07 | ∇²φ = q; the response to a unit point source G = −1/(4π\|x − x′\|) solves ∇²G = δ (flux of ∇(1/r) through a small sphere = −4π); superposition φ = ∫Gq d³x′ (glosses: Dirac delta as a point source, superposition) |
| gradient of 1/distance with respect to the source point (P140) | ch05 · §5.5 · C07 | ∇′(1/\|x − x′\|) = +(x − x′)/\|x − x′\|³ points from the source to the field point; ∇ with respect to x flips the sign — the book's two §5.5 slips |
| curl of a product (P141) | ch05 · §5.5 · C07 | ∇×(fA) = f∇×A + ∇f × A, order of the cross product kept (gloss: a × b = −b × a); moves the curl off ω in D11 (was a ch04 gloss for χ∇ψ) |
| improper integral as a limit (P144) | ch05 · §5.5 · C08 | integrate to a finite end, then let it run away: ∫₀^∞ e^{−x}dx = lim(1 − e^{−b}) = 1; a segment grows into an infinite line in D13 |
| Frenet frame of a curve (P145) | ch05 · §5.6 · C10 | unit tangent, principal normal (toward the centre of curvature), binormal; ⚠️ the book's e_n points **away** (= −N); helix curvature a/(a² + c²), torsion c/(a² + c²) (gloss: plotly dropdown menus for the frame figure) |
| **ch05 — physics vocabulary** | | |
| angular momentum of a spinning cylinder (P146) | ch05 · §5.6 · C10 | L = IΩ with I = ½mR² = mA/2π; no torque ⇒ L fixed, so stretching at fixed mass and volume shrinks A and spins it up ∝ length (the skater; vortex stretching D17, D18) |
| **ch06 — maths** | | |
| 2-D divergence theorem and the Dirac delta in the plane (P149) | ch06 · §6.2 · C02 | ∫_A∇²f dA = ∮∇f·n ds; δ(x)δ(y) is zero except at the origin with total 1; a function harmonic away from 0 whose flux through every circle is the same F has ∇²f = Fδ — ∇ ln r = e_r/r has flux 2π (D03: vortex and source as deltas) |
| a limit with a product held fixed (P150) | ch06 · §6.3 · C04 | ε → 0 and m → ∞ together with 2mε = \|d\| fixed: alone each limit gives 0 or ∞, together the leading term survives and the rest shrinks like ε² (source–sink pair → doublet, D06) |
| integrals of sines and cosines over a full period (P151) | ch06 · §6.3 · C06 | over 0…2π, sin, sin³, cos, sin cos, sin²cos all integrate to 0; sin² and cos² give π (mean ½); force integrals round a circle keep exactly the surviving terms (L = ρUΓ, D11) |
| Newton's method for complex zeros (P152) | ch06 · §6.3 · C07 | z ← z − f/f′, digits double near a simple root; stagnation points solve dw/dz = u − iv = 0 with f′ = d²w/dz²; deflation f/(z − z₁) finds the next root; ⚠️ the basin is the size of the flow's own length scale (m/2πU) — seed there |
| the complex plane in numpy (P153) | ch06 · §6.4 · C09 | z = x + iy = re^{iθ}; \|z\| distance, arg z = atan2(y, x); multiplying multiplies moduli and **adds angles** (rotate + stretch); z* mirrors in x; numpy `1j`, `np.abs`, `np.angle`, `np.conj` on whole arrays |
| complex derivative and analytic functions (P154) | ch06 · §6.4 · C09 | dw/dz = lim [w(z + δz) − w(z)]/δz must be the same for every direction of δz: then w is analytic (polynomials, exp, log off its cut, 1/z off 0); z* fails (quotient 1 along x, −1 along iy) — the root of Cauchy–Riemann (D14) |
| complex logarithm, powers and branch cuts (P155) | ch06 · §6.4 · C09 | ln z = ln r + iθ is many-valued; numpy picks θ ∈ (−π, π] so ln z jumps by 2πi across the negative real axis, and zⁿ = e^{n ln z} inherits the cut; rotate the cut into the body or outside the drawn fluid (`log_branch(cut=)`, `power_branch`) |
| Cauchy's integral theorem (P156) | ch06 · §6.5 · C10 | f analytic inside and on C ⇒ ∮f dz = 0; so two contours round the same singular region give the same ∮ — a contour may be squeezed or stretched without crossing a singularity; on a circle dz = iz dθ (why Blasius's contour can leave the body, D17) |
| Laurent series and residues (P157) | ch06 · §6.5 · C10 | outside all singular points f = Σc_kz^k with negative powers; term by term round a circle only z⁻¹ survives (2πi), so ∮f dz = 2πi c₋₁ (the residue); an FFT of samples on a circle reads every c_k (`laurent_coefficients`) |
| complex square roots and the quadratic formula (P159) | ch06 · §6.6 · C11 | ζ² − zζ + b² = 0 has roots ½[z ± √(z² − 4b²)] with product b² (one inside, one outside \|ζ\| = b); `np.sqrt` gives the principal root (real part ≥ 0, cut on the negative reals) and does **not** always pick the outside one — use √(z − 2b)√(z + 2b) (D21) |
| surface integrals on a sphere (P163) | ch06 · §6.9 · C15 | dA = a² sin θ dθ dφ, e_r = (sin θ cos φ, sin θ sin φ, cos θ); integrate φ first (cos φ, sin φ vanish); ∫₀^π cos²θ sin θ dθ = 2/3; `dblquad(f, a, b, c, d)` integrates f(inner, outer) (added mass, D30) |
| Green's first identity (P164) | ch06 · §6.9 · C15 | ∇·(φ∇φ) = \|∇φ\|² + φ∇²φ; with ∇²φ = 0 the kinetic-energy density is a divergence, so Gauss turns ½ρ∫\|∇φ\|²dV into surface integrals over the body (normal out of the fluid, into the body) and a far sphere (D31) |
| **ch06 — Python and numerics** | | |
| sympy for complex series and residues (P158) | ch06 · §6.5 · C10 | `sp.I` is i, `sp.expand` multiplies out a squared series, `.coeff(z, -1)` reads one coefficient, `sp.residue(f, z, 0)` returns the 1/z coefficient (the (6.61) slip check, D18) |
| iterative solvers: Jacobi, Gauss–Seidel, SOR (P160) | ch06 · §6.7 · C12 | sweeps instead of elimination: Jacobi uses last sweep's values, Gauss–Seidel the newest, SOR scales the Gauss–Seidel change by 1 < ω < 2; error shrinks by the spectral radius ρ per sweep; stop on the **residual** b − Aψ with a tight tolerance — when ρ ≈ 1 residual and change are both ≈ (1 − ρ) × the true error |
| boolean masks and scipy.sparse (P161) | ch06 · §6.7 · C12 | an L-shaped or stepped domain = a rectangular array plus `mask[j, i]` (True = unknown) and a boundary-value array; `scipy.sparse` stores the five non-zeros per row, `spsolve` solves the whole system directly on fine grids |
| collocation and the condition number (P162) | ch06 · §6.8 · C14 | N unknowns from the condition at N points; `np.linalg.cond(A)` bounds how much relative data error can grow: ~10³ harmless, > 10¹⁰ the digits are noise (`np.vander` as the classic bad matrix; the axial method reaches 10¹⁷) |
| **ch07 — maths** | | |
| phase of a wave (P165) | ch07 · §7.1 · C01 | the phase kx − ωt is the argument of the cosine in radians: 0 at a crest, π at a trough, 2π at the next crest; k [rad/m] and ω [rad/s] vs cycles (1/λ, ν) differ by 2π; a point of fixed phase moves right when ω/k > 0 (with the ⚠️ "ω is a frequency now, not a vorticity" callout) |
| Taylor transfer of a boundary condition (P166) | ch07 · §7.2 · C02 | a condition on the moving surface z = η is written on z = 0 by F(η) ≈ F(0) + ηF′(0); if F varies over 1/k and η ~ a the correction is ≈ ka × F, so gentle waves keep F(0) — and the two O(ka) approximations (drop η_xφ_x, move to z = 0) must be made together |
| separation of variables for a PDE (P167) | ch07 · §7.2 · C03 | for a linear PDE with straight boundaries try (function of z) × (wave in x): ∂²/∂x² → −k² turns Laplace into the ODE f″ − k²f = 0 (e^{λz} trial, ch01 P44); the boundary conditions choose sin or cos (contrast ch01 P42, separation in the ODE sense) |
| hyperbolic functions cosh, sinh, tanh (P168) | ch07 · §7.2 · C03 | cosh = (eˣ + e⁻ˣ)/2, sinh = (eˣ − e⁻ˣ)/2, tanh, coth, sech; cosh² − sinh² = 1, cosh² + sinh² = cosh 2x, 2 sinh cosh = sinh 2x, (tanh)′ = sech²; small x: cosh ≈ 1, sinh ≈ tanh ≈ x; large x: cosh ≈ sinh ≈ eˣ/2, tanh → 1 (tanh 2 = 0.964 — the deep-water threshold); our Fig. 7.7 |
| curvature of a plane curve (P169) | ch07 · §7.3 · C07 | for z = η(x), 1/R = η_xx/(1 + η_x²)^{3/2} ≈ η_xx for gentle slopes; under a crest η_xx < 0 and the centre of curvature is below, inside the water, so the liquid pressure is higher (p = −ση_xx, D15) |
| sum-to-product identities (P171) | ch07 · §7.4 · C08 | cos A + cos B = 2 cos((A − B)/2) cos((A + B)/2), cos A − cos B = −2 sin((A + B)/2) sin((A − B)/2): a sum of two waves becomes a slow factor times a fast one — standing waves (D17) and beats (D19) |
| Fourier integral and a packet's spectrum (P172) | ch07 · §7.5 · C09 | a group is η = ∫A(k)e^{i(kx − ω(k)t)}dk (real part understood); a spectrum of width δk around k₀ means a group of length ~1/δk, and near k₀ ω(k) may be Taylor-expanded — the continuous version of ch05 P142's FFT modes (D20) |
| first-order wave equation and characteristics (P174) | ch07 · §7.5 · C10 | ∂q/∂t + c ∂q/∂x = 0 keeps q constant for an observer moving at c: along each line dx/dt = c of the x–t plane (a characteristic); the material derivative of ch03 with u → c (crest conservation D22, the simple wave N81) |
| Snell's law for waves (P175) | ch07 · §7.5 · C10 | where the medium changes only along x, the along-y wavenumber cannot change (crests match along every line x = const); with \|k\| changing, \|k\| sin α = const (α from the x-direction) — optics' n₁ sin α₁ = n₂ sin α₂ (refraction D24) |
| complex amplitudes (P176) | ch07 · §7.7 · C13 | ζ = Re{a e^{i(kx − ωt)}} with the Re dropped during linear algebra: ∂/∂x → ik, ∂/∂t → −iω, a complex amplitude b = \|b\|e^{iφ} carries size and phase shift; ⚠️ take real parts before multiplying two fields (ch01 P45 Euler, ch06 P153) |
| operator elimination for linear PDEs (P177) | ch07 · §7.8 · C15 | apply ∂/∂t, ∂/∂z, ∇_H² to whole linear equations and substitute until one unknown is left — like eliminating variables in a linear system; allowed because derivatives of smooth fields commute (Schwarz, ch04 P121), also with N(z) (the w-equation D33) |
| mean of a product of real parts (P178) | ch07 · §7.8 · C16 | ⟨Re(Ae^{iθ})Re(Be^{iθ})⟩ = ½Re(AB*) over a period (B* conjugate, ch02 P81); fields 90° apart (B = iA) give zero mean product — why ⟨gρ′w⟩ = 0 and how E_k, E_p, F of an internal wave are computed (D37) |
| **ch07 — Python and numerics** | | |
| `scipy.optimize.minimize_scalar` (P170) | ch07 · §7.3 · C07 | the minimum of a function of one variable inside a bracket (`method="bounded"`); used to check the closed-form c_min and c_g,min by brute force |
| envelope with `scipy.signal.hilbert` (P173) | ch07 · §7.5 · C09 | `np.abs(hilbert(eta))` follows the slowly varying amplitude of a fast oscillation without fitting; tracks a packet's peak (`ch07.envelope`, `core.waves.envelope`) |
| `np.stack` and `np.c_` (P179) | ch07 · §7.1 · C01 | `np.stack([X, Y], axis=-1)` glues equal-shaped arrays on a new last axis (each grid point carries (x, y): the layout `plane_wave` expects); `np.c_[a, b]` stacks 1-D arrays as columns (two points → the x-row and y-row that `ax.plot(*…)` needs) |
| `np.sign` and `np.nonzero`: where a curve crosses zero (P180) | ch07 · §7.5 · C09 | `np.sign` → −1, 0, +1; neighbours `e[:-1]`, `e[1:]` of different sign bracket a crossing; `np.nonzero(mask)[0]` lists the indices where a boolean array is True (a tuple per axis, hence `[0]`) |
| `np.fft.rfft` and `np.fft.rfftfreq` (P181) | ch07 · §7.5 · C09 | for a real signal only the k ≥ 0 half is returned; `rfftfreq(n, d)` gives cycles per unit length, × 2π → wavenumbers; each mode then moves with its own ω(k) (`linear_evolve`; FFT itself ch05 P142) |
| `np.interp`: reading a curve between samples (P182) | ch07 · §7.5 · C09 | `np.interp(x_new, x, y)` draws a straight line between the neighbouring samples (x increasing); reads the envelope's height at a crest that sits between grid points |
| `warnings.catch_warnings(record=True)` (P183) | ch07 · §7.7 · C13 | a warning is a message without stopping; inside the `with` block every warning is collected in a list (`simplefilter("always")` so none is swallowed) — shows `interface_omega` flagging ρ₁ > ρ₂ (Rayleigh–Taylor) before returning NaN |
| `np.errstate(invalid="ignore")` (P184) | ch07 · §7.8 · C15 | silences numpy's RuntimeWarning (0/0 → NaN) for the lines inside the block only; used where the NaN is expected and harmless (at K = 0 the direction of K, hence ω = N cos θ, is undefined) |
| **ch08 — maths and physics** | | |
| diffusivity and the diffusion time L²/ν (P185) | ch08 · §8.1 · C01 | a diffusivity D [m²/s] spreads something a distance ~√(Dt) in time t, so crossing L takes ~L²/D (double L, quadruple the time); heat has κ = k/ρC_p, momentum has ν = μ/ρ — why air (ν ≈ 1.5e-5) spreads momentum 15× faster than water; the √(νt) of every later boundary and Ekman layer |
| Laplacian in cylindrical coordinates (P186) | ch08 · §8.2 · C03 | for u(R): ∇²u = (1/R)d/dR(R du/dR) (the R inside comes from the growing circumference); for a swirl u_φ(R) the φ-component of the vector Laplacian is d/dR[(1/R)d(Ru_φ)/dR] (extra −u_φ/R² because e_φ turns); both from `core.curvilinear` |
| Euler–Cauchy (equidimensional) ODE (P187) | ch08 · §8.2 · C04 | every term R^k d^ku/dR^k ⇒ try u = R^λ, get a polynomial in λ: R²u″ + Ru′ − u = 0 gives λ = ±1, u = AR + B/R; reused for Stokes' f(r) (roots 4, 2, 1, −1) |
| anisotropic scaling with two length scales (P188) | ch08 · §8.3 · C05 | long thin flows: x in units of L, y in units of h = εL, so ∂/∂x = (1/L)∂/∂x* but ∂/∂y = (1/εL)∂/∂y*; the cross velocity scale comes from continuity (v ~ εU); each term's size then sits in its coefficient — the lubrication and boundary-layer move |
| Leibniz rule with a moving upper limit (P189) | ch08 · §8.3 · C06 | ∂/∂x ∫₀^{h(x)} u dy = ∫₀^h ∂u/∂x dy + u(x, h)∂h/∂x: the extra term counts what enters because the limit moved (ch03 P109 had fixed limits); sympy demo returns 0 |
| nonlinear diffusion in flux form (P191) | ch08 · §8.3 · C08 | ∂h/∂t = ∂/∂x(D(h)∂h/∂x) with D depending on the unknown (thin film: D = ρgh³/3μ, tiny where the layer is thin); written as h_t + q_x = 0 with q = −D h_x, the total ∫h dx is conserved when nothing leaves the ends |
| Gaussian integral (P194) | ch08 · §8.4 · C09 | ∫_{−∞}^{∞}e^{−ζ²}dζ = √π, √π/2 on each side; no elementary antiderivative, which is why its running integral (erf) has a name |
| error function erf and its inverses (P195) | ch08 · §8.4 · C09 | erf(ζ) = (2/√π)∫₀^ζ e^{−ξ²}dξ rises 0 → 1, erfc = 1 − erf computed directly (no cancellation; ch04 P123); `erfinv`/`erfcinv` answer "where does the profile reach this level?" — 2 erfcinv(0.01) = 3.643 is δ₉₉ in units of √(νt) |
| exponent matching for similarity forms (P197) | ch08 · §8.4 · C10 | if c₁t^aF(η) + c₂t^bG(η) = 0 for every t and η then a = b; power laws δ = Dt^m turn every bracket into a power of t, matching gives linear equations for the exponents (sympy: the spreading bead's n = m = 1/5) |
| dominant balance (P198) | ch08 · §8.6 · C12 | choose a quantity's scale from the term it must balance: at low Re pressure balances viscous stress, p − p∞ ~ μU/L, not ρU² (10 µm particle at 1 mm/s: 0.1 Pa vs 1e-3 Pa); the wrong scale makes a term look negligible |
| the Stokes operator E² applied twice (P199) | ch08 · §8.6 · C13 | E² = ∂²/∂r² + (sin θ/r²)∂/∂θ((1/sin θ)∂/∂θ) (Ch. 6 (6.77)) turns ψ into vorticity, ω_φ = −E²ψ/(r sin θ); E²(E²ψ) is **not** the biharmonic ∇⁴ψ (a sympy test shows the difference); E² of the uniform stream r²sin²θ is 0 |
| **ch08 — Python and numerics** | | |
| `scipy.integrate.cumulative_trapezoid` (P190) | ch08 · §8.3 · C06 | `cumulative_trapezoid(f, x, initial=0)` returns the running integral at every sample — the numerical antiderivative; turns a pressure gradient into a pressure profile (slider from scratch) |
| implicit time stepping with Picard iteration (P192) | ch08 · §8.3 · C08 | backward Euler evaluates the right side at the new time (stable for any step: damping 1/(1 + 4λ) < 1 where explicit needs λ ≤ ½); with a solution-dependent diffusivity, guess it from the last iterate, solve, update, repeat; a precursor film h_min keeps D > 0 ahead of the front (`thin_film_spread` uses Newton, same fixed point) |
| Crank–Nicolson with `scipy.linalg.solve_banded` (P193) | ch08 · §8.4 · C09 | average the diffusion term between old and new time: second order in time, stable for any step; each step is a tridiagonal system that `solve_banded((1, 1), ab, rhs)` solves in O(N) from the three diagonals stored as rows (`core.diffusion.crank_nicolson_1d`; two backward-Euler start-up steps tame the jump) |
| `scipy.integrate.solve_bvp` (P196) | ch08 · §8.4 · C09 | an ODE with conditions at both ends (F(0) = 1, F(η_max) = 0 on a truncated domain), written as a first-order system with a mesh and a guess; knows nothing about erf yet lands on erfc(η/2) (3.7e-13) — the check that η_max is large enough |
| **ch09 — maths and physics** (the term is the exact `nb.primer` title) | | |
| six Reynolds numbers (which length?) (P200) | ch09 · §9.1 · C01 | Re is inertia over viscosity, but only with its length and speed: overall U∞L/ν (9.6), local Re_x = Ux/ν, plate Re_L, cylinder/sphere on the **diameter**, jet xu₀/ν and Re_h99, Falkner–Skan ax^{n+1}/ν; "Re ≫ 1" means nothing without the length |
| parabolic, elliptic and marching in x (P201) | ch09 · §9.1 · C01 | elliptic problems (full Navier–Stokes) need conditions all round — downstream talks to upstream; parabolic ones have a time-like direction and a start (heat equation; (9.9) with x as time), so they are solved by marching from an inlet profile and information flows downstream only |
| improper integral of a deficit (truncating the tail) (P202) | ch09 · §9.2 · C02 | ∫₀^∞(1 − u/U_e)dy converges because the deficit dies like a Gaussian or an exponential; integrate numerically to a y_max where it is negligible and add the tail (an exponential tail exactly) |
| control volume with a streamline as a side (P205) | ch09 · §9.2 · C02 | no fluid crosses a streamline, so as a face of a control volume it carries no mass or momentum flux and only the ambient pressure: the roof of the box over a wall layer (D04's ρU²θ = ∫τ₀dx) |
| chain rule when the similarity variable moves with x and y (P206) | ch09 · §9.3 · C03 | η = y/δ(x): ∂F/∂y = F′/δ, ∂F/∂x = −F′ηδ′/δ; a product δ(x)f(η) also needs the product rule in x — the step the book skips in every reduction (Blasius, Falkner–Skan, both jets) |
| scaling symmetry of an ODE and the Töpfer trick (P207) | ch09 · §9.3 · C04 | if λf₁(λη) solves the ODE whenever f₁ does, solve once with f″(0) = 1, read f′(∞) = 2.085 instead of 1 and rescale by λ² = 1/2.085: a BVP becomes one IVP (Töpfer 1912); reused for the wall jet's free scale f∞ |
| integration by parts (P218a) | ch09 · §9.3 · C04 (before D06) | the product rule (ab)′ = a′b + ab′ integrated: ∫a b′ = [ab] − ∫a′b moves a derivative to the other factor at the price of a boundary term; in D06 a = f, b = f′ − 1 and the wall term vanishes because f(0) = 0 (runnable `quad` demo, 0.38177 twice) |
| continuation in a parameter and a fold (saddle-node) (P209) | ch09 · §9.4 · C05 | follow a family f(η; n) by small steps, each solution the next guess; a branch can turn back at a fold where two solutions merge and vanish — beyond it `solve_bvp` fails because nothing exists (Falkner–Skan n = −0.0904); parametrise by f″(0) to pass it |
| first-order linear ODE and the integrating factor (P210) | ch09 · §9.6 · C07 | y′ + p(x)y = q(x): multiply by μ = e^{∫p} so the left side is (μy)′, then integrate; Thwaites' p = 6U_e′/U_e gives μ = U_e⁶ (sympy `dsolve` demo) |
| integrals of powers of sine (∫sin⁵ by c = cos φ) (P211) | ch09 · §9.6 · C07 (D11) | sin⁵φ dφ = (1 − cos²φ)²sin φ dφ and c = cos φ turn it into a polynomial: ∫₀^φ sin⁵ = 8/15 − c + (2/3)c³ − c⁵/5 (the cylinder's Thwaites integral) |
| inflection point (P212) | ch09 · §9.7 · C08 | where the second derivative changes sign — the curve stops bending one way; for a velocity profile u_yy = 0 with a sign change; profiles with one are the unstable ones of Ch. 11 (Rayleigh) |
| a row of vortices: the cotangent sum (P213) | ch09 · §9.8 · C10 | a point vortex's conjugate velocity (Γ/2πi)/(z − z₀) summed over a row at z₀ + na, pairing +n with −n so it converges, is (Γ/2ia)cot(π(z − z₀)/a); companion lattice sums Σ1/(z − na)² = (π/a)²/sin²(πz/a) and the alternating one |
| linear stability of a steady configuration (perturb, linearise, eigenvalues) (P214) | ch09 · §9.8 · C10 | displace slightly, keep terms linear in the displacement: ẋ = Mx; solutions grow like e^{λt}, so any eigenvalue with Re λ > 0 means unstable, all Re λ ≤ 0 stable/neutral (`np.linalg.eig`; extends ch03 P80 eigenvalues) — the method of Ch. 11 |
| momentum flux and mass flux through a cross-section (P217) | ch09 · §9.10 · C12 | per unit span ṁ = ρ∫u dy [kg/(m s)] and J = ρ∫u²dy [N/m]; a jet in still fluid at constant pressure feels no force, so J cannot change, but ṁ can — fluid enters through the sides (entrainment) |
| sech, arccosh and (tanh)′ = sech² (P215) | ch09 · §9.10 · C12 | sech x = 1/cosh x, a bell equal to 1 at 0 dying like 2e^{−\|x\|}; (tanh x)′ = sech²x = 1 − tanh²x; arccosh y = ln(y + √(y² − 1)) answers "where does sech² reach 1 %?" (extends the ch07 tanh primer) |
| the total-derivative move (look for the derivative of a product) (P216) | ch09 · §9.10 · C12 | before integrating an ODE ask whether it already is a derivative: (ff′)′ = f′² + ff″, (3f′ + f²/2)′ = 3f″ + ff′ (product rule backwards, ch01 P38); integrating is then free |
| integration by parts with a variable lower limit (P218) | ch09 · §9.10 · C13 | recap of P218a, then ∫₀^∞u G dy with G(y) = ∫_y^∞g dy′: G′ = −g, so ∫uG = [WG] + ∫Wg with W = ∫₀^y u; d/dx of a double integral with fixed limits only differentiates the integrand — the wall-jet invariant (9.80) |
| partial fractions and sympy apart; inverting an implicit solution (P219) | ch09 · §9.10 · C13 | split 1/(1 − g³) into simple pieces that integrate (`sympy.apart`); when the integral gives η(g) but we want g(η), `brentq` finds g with η(g) − η = 0 for each η (the wall jet (9.83)) |
| radial force balance in a swirl (and a thin layer) (P220) | ch09 · §9.11 · C14 | a parcel on a circle needs inward force ρu²/R; in the fast core the pressure gradient supplies it, ∂p/∂R = ρu_e²/R; across a thin layer ∂p/∂z ≈ 0 (the (9.10) argument), so slower fluid feels more push than it needs: net inward ρ(u_e² − u²)/R (1000 N/m³ at the floor for u_e = 0.2 m/s, R = 4 cm) |
| **ch09 — Python and numerics** | | |
| `scipy.integrate.simpson` and `np.trapezoid` (P203) | ch09 · §9.2 · C02 | both integrate a sampled profile: trapezoid joins samples by lines (error ∝ Δy²), Simpson fits parabolas through triples (∝ Δy⁴); numpy 2 has no `np.trapz` |
| monotone interpolation `PchipInterpolator` (P204) | ch09 · §9.2 · C02 | reads a height such as u/U = 0.99 from samples without overshoot (a cubic spline may wiggle), then `brentq` finds the crossing; also interpolates the tabulated closure l(λ), H(λ) |
| shooting versus boundary-value solving (P208) | ch09 · §9.3 · C04 | shooting guesses the missing slope, integrates and adjusts with `brentq` on f′(η_max) − 1; `solve_bvp` treats the whole interval; Blasius shoots easily, near the Falkner–Skan fold the root is ill-conditioned (sensitivity ~e^{η²/4}), so shoot only for n ≥ −0.05 and use `solve_bvp` with continuation elsewhere |
| **ch10 — maths and physics** (the term is the exact `nb.primer` title) | | |
| big-O notation and the order of accuracy (P221) | ch10 · §10.2 · C01 | O(Δx²) = "at most a constant times Δx² for small Δx" — halving Δx divides it by 4; a stencil is p-th order when its error is O(Δx^p); unlike Ch. 2's orders of smallness (P68) the leftover is kept and measured, not dropped |
| half-angle identities (P225) | ch10 · §10.2 · C04 | 1 − cos θ = 2 sin²(θ/2), sin θ = 2 sin(θ/2)cos(θ/2), so sin²θ = 4s(1 − s) with s = sin²(θ/2): a cosine running 1 → −1 becomes a square running 0 → 1 — the move behind \|G\|² (10.26) and Noye's region |
| sign of a linear function on an interval (P226) | ch10 · §10.2 · C04 | a line a + bs is ≤ 0 on an interval iff it is ≤ 0 at both ends (a line cannot bulge); one "for all s" inequality becomes two readable ones (test the limit at an open end) |
| domain of dependence and characteristics (P227) | ch10 · §10.2 · C05 | for T_t + uT_x = 0 the value at (x, t) came from x − uΔt one step earlier (characteristics, Ch. 7 P174); a scheme can only be right if its stencil contains that point — the CFL condition |
| well-posed problem (P228) | ch10 · §10.2 · C05 | Hadamard: a solution exists, is unique and depends continuously on the data; the heat and advection equations with sensible BCs are, the backward heat equation is not; the Lax theorem assumes it |
| Péclet number (global R and cell R_cell) (P229) | ch10 · §10.2 · C05 | advection over diffusion u ℓ/D: with ℓ = L the global R = uL/D (10.87) (the Reynolds number's twin for a scalar), with ℓ = Δx the cell number R_cell = uΔx/D (10.31) — how advective one grid cell is |
| test functions and the spaces H¹, S and V (P230) | ch10 · §10.3 · C06 | a test function w weights the equation; a weak statement holds "for every w"; H¹ = functions with square-integrable slope (kinks yes, jumps no); trial space 𝒮 has the Dirichlet value built in, test space V vanishes there |
| fundamental lemma of the calculus of variations (P231) | ch10 · §10.3 · C06 | if ∫f w dx = 0 for every smooth w vanishing at the ends, then f = 0 (choose w = f × a bump); "for all w" pins f down point by point |
| bilinear form a(w, v) (P232) | ch10 · §10.3 · C07 | linear in each slot separately, so sums and constants come out like from an integral (10.50) → (10.52); here not symmetric (convection differentiates only v) |
| arbitrary coefficients: every bracket is zero (P233) | ch10 · §10.3 · C07 | if Σc_AG_A = 0 for every choice of c_A, pick unit vectors: each G_A = 0 — one weak statement becomes n equations |
| affine map to a parent element (P234) | ch10 · §10.3 · C08 | x(ξ) = (h/2)ξ + x_mid sends [−1, 1] onto [x_a, x_b]; dx = (h/2)dξ (substitution P106), dξ/dx = 2/h (chain rule P49); every element integral becomes one over [−1, 1] |
| linear recurrence with constant coefficients (geometric trial) (P236) | ch10 · §10.4 · C09 | aT_{j+1} + bT_j + cT_{j−1} = 0 is the discrete twin of a constant-coefficient ODE (P44): T_j = r^j gives ar² + br + c = 0; the ends fix the two constants; a **negative root makes r^j flip sign each node** |
| Lagrange multiplier (P239) | ch10 · §10.4 · C10 | an extra unknown times a constraint, adjusting itself until the constraint holds; its value measures how hard the constraint pushes — in incompressible flow the pressure is the multiplier of ∇·u = 0 |
| the Taylor–Green vortex (an exact decaying Navier–Stokes solution) (P254) | ch10 · §10.4 · C10 | u = sin x cos y e^{−2t/Re}, v = −cos x sin y e^{−2t/Re}, p = ¼(cos 2x + cos 2y)e^{−4t/Re} in a 2π-periodic box: divergence-free, (u·∇)u balanced exactly by −∇p, pure viscous decay — the V1 test field of MAC and MacCormack (added by the lesson review: used before explained) |
| conservation (flux) form U_t + E_x + F_y = 0 (P237) | ch10 · §10.4 · C10 | stack the conserved quantities U = (ρ, ρu, ρv); each law is "rate + divergence of a flux = 0" with E = (ρu, ρu² + p, ρuv), F = (ρv, ρuv, ρv² + p) (+ viscous terms); one scheme updates all components and a flux sum telescopes (conservation) |
| predictor–corrector (Heun's second-order idea) (P238) | ch10 · §10.4 · C10 | predict with a cheap first-order step, re-evaluate the slope there, average: y* = yⁿ + Δt f(yⁿ), yⁿ⁺¹ = yⁿ + (Δt/2)[f(yⁿ) + f(y*)] — second order, like the Runge–Kutta family (P95) |
| operator splitting and the commutator [A₁, A₂] (P241) | ch10 · §10.4 · C11 | step with A₁ alone, then A₂ alone: e^{−ΔtA₂}e^{−ΔtA₁} vs e^{−Δt(A₁+A₂)} differ by (Δt²/2)[A₁, A₂] per step for matrices — O(Δt²) per step, O(Δt) overall (Marchuk–Yanenko first order); Θ-scheme cancels it at Θ = 1 − 1/√2 |
| Helmholtz–Hodge decomposition (P240) | ch10 · §10.4 · C11 | any smooth field = divergence-free part + gradient; the divergence gives ∇²φ = ∇·w (Poisson, P139), so one solve finds φ and u = w − ∇φ; the gradient part is curl-free (ch02) — the projection method in one line |
| half-index notation and staggered array shapes (P242) | ch10 · §10.4 · C12 | u_{i+1/2, j} lives on the face between cells (i, j) and (i + 1, j); code has no half indices: p[ny, nx] (centres), u[ny, nx+1] (vertical faces incl. both walls), v[ny+1, nx]; `u[j, i]` = u_{i−1/2, j} (the `[j, i]` layout, P76) |
| singular linear systems and the compatibility condition (P243) | ch10 · §10.4 · C12 | the pure-Neumann Poisson matrix sends a constant to zero: Ax = b is solvable only if Σb = 0 (net inflow zero) and then up to a constant — pin one value or ask for zero mean; `np.linalg.solve` refuses, `lstsq` gives the min-norm answer |
| null space by SVD (count the tiny singular values) (P245) | ch10 · §10.4 · C12 | A = UΣVᵀ; columns of V with (numerically) zero singular values span the inputs A cannot see; `np.linalg.svd` + a tolerance counts them (checkerboard: collocated 4, staggered 1) |
| saddle-point (KKT) matrix (P246) | ch10 · §10.4 · C13 | [[A, B], [Bᵀ, 0]]: every "minimise subject to a constraint" system (the multiplier is the second block); indefinite, zero diagonal block, invertible only if B has full column rank — no invisible multiplier pattern; `np.block` builds it, `eigvalsh` shows the sign mix |
| generalised symmetric eigenproblem scipy.linalg.eigh(A, B) (P247) | ch10 · §10.4 · C13 | A x = λ B x for symmetric A, positive-definite B: eigenvalues of A measured in units of B (P80 was B = I); the discrete inf–sup constant β_h² = λ_min(BᵀA⁻¹B, M_p) over zero-mean pressures |
| barycentric (area) coordinates on a triangle (P248) | ch10 · §10.4 · C13 | a point = weighted average of the corners with weights (ζ, ξ, η) ≥ 0 summing to 1 (ζ = 1 − ξ − η); each is 1 at its corner and 0 on the opposite side = the linear shapes (10.187); products like 4ξζ are the P2 mid-edge shapes (10.185) |
| Jacobian determinant of a 2-D map (P249) | ch10 · §10.4 · C13 | (ξ, η) → (x, y) stretches dξdη into J dξdη with J = x_ξy_η − x_ηy_ξ — the 2-D dx = (h/2)dξ; a straight triangle has J = 2 × area, a curved (isoparametric) one varies inside |
| verification versus validation (P253) | ch10 · §10.5 · C14 | verification = solving the equations right (bugs, grid and Δt errors; exact, manufactured and benchmark solutions such as Ghia's); validation = solving the right equations (experiments; Ch. 12) |
| **ch10 — Python and numerics** | | |
| floating-point round-off and machine epsilon (P222) | ch10 · §10.2 · C01 | double precision keeps ~16 digits, `np.finfo(float).eps` = 2.2e-16; a difference quotient's round-off ≈ ε\|f\|/h grows as h shrinks while truncation falls as h^p — the total error has a floor (near h ≈ ε^{1/2} for a first-order stencil, ε^{1/3} for a centred one); in float32 the forward floor is ≈ 1e-4 at h ≈ 3e-4, the centred ≈ 1e-5 at h ≈ 5e-3 |
| ghost node for a Neumann boundary (P223) | ch10 · §10.2 · C02 | invent T_{N+1} so the centred difference at the end gives the slope: T_{N+1} = T_{N−1} + 2Δx q; the ordinary FTCS formula then runs at the last node and stays second order; q = 0 is an insulated end |
| norms of an error array: rms, max, L1 (P224) | ch10 · §10.2 · C03 | one number for an error array: rms √mean(e²) (the book's (10.15)), max \|e\| (the worst point, honest near a jump or wiggle), mean \|e\|; for smooth errors all shrink at the same rate |
| scatter-add assembly (np.add.at, COO duplicates) (P235) | ch10 · §10.3 · C08 | add small blocks into a big matrix where indices repeat: `np.add.at(K, (rows, cols), block)` adds every contribution (plain `+=` keeps the last); sparse: collect triples, `coo_matrix(...).tocsr()` sums duplicates |
| scipy.sparse.diags, kron and a cached splu factorisation (P244) | ch10 · §10.4 · C12 | 1-D second difference = `sparse.diags([1, -2, 1], [-1, 0, 1])`; 2-D 5-point Laplacian = `kron(I_y, L_x) + kron(L_y, I_x)` (P161); the matrix never changes, so `lu = splu(A.tocsc())` once and `lu.solve(b)` per step |
| GMRES in one line (P250) | ch10 · §10.4 · C13 | for big nonsymmetric systems (Newton steps), GMRES builds the best answer in span{b, Ab, A²b, …} and stops at a small residual; `scipy.sparse.linalg.gmres(A, b)` → (x, info), info = 0 converged |
| reading reference data from a file and interpolating (np.interp) (P251) | ch10 · §10.5 · C14 | benchmark tables live in `reference/chNN/*.csv` with the citation in the header; read with `np.loadtxt`/`pandas.read_csv`, then interpolate **our** profile at **their** points (`np.interp(their_y, our_y, our_u)`, P182) so no benchmark value is invented |
| caching expensive runs (np.savez and a parameter key) (P252) | ch10 · §10.5 · C14 | save a long run once (`np.savez` in `outputs/chNN/`, parameters in the file name or a hash) and `np.load` it next time; `functools.lru_cache` in memory; the notebook says whether a result was loaded or computed |
| **ch11 — maths and physics** (the term is the exact `nb.primer` title) | | |
| necessary vs sufficient conditions, "for every k", and proof by contradiction (P255) | ch11 · §11.2 · C01 | sufficient: A guarantees B (Ri > ¼ everywhere ⇒ stable); necessary: no B without A, but A alone does not force B (an inflection point for inviscid instability); "stable" = σ_r ≤ 0 for **every** k, "unstable" needs **one** k; contradiction: assume a growing mode, reach something impossible |
| eigenvalue problem for a differential equation (P257) | ch11 · §11.4 · C03 | an ODE with boundary conditions has a nonzero solution only for special parameter values (−u″ = λu, u(0) = u(π) = 0 ⇒ λ = 1, 4, 9, …, u = sin nz); those values are the eigenvalues, the solutions the eigenfunctions (mode shapes) — extends the matrix case of P80 |
| Chebyshev–Gauss–Lobatto points and the differentiation matrix (P258) | ch11 · §11.4 · C03 | sample at x_j = cos(jπ/N) (crowded near the ends, `x[0]` is the right/top end) and differentiate the interpolating polynomial exactly with a matrix D (`D @ D` for the second derivative); for smooth solutions the error falls exponentially with N — digits, not an order |
| real and imaginary parts of a complex identity (P260) | ch11 · §11.4 · C03 | one complex equation = two real ones; ∫∣f∣² dz > 0 unless f ≡ 0; 1/(U − c) = (U − c*)/∣U − c∣², so Im 1/(U − c) = c_i/∣U − c∣²; (real quantity) × c_i = 0 forces c_i = 0 or the quantity to vanish — the engine of D08, D18, D19, D22 |
| even and odd functions; parity under z → −z (P261) | ch11 · §11.4 · C04 | even f(−z) = f(z) (cos, cosh), odd f(−z) = −f(z) (sin, sinh); a problem symmetric under z → −z has even or odd eigenfunctions, found separately (even W one row of cells, odd W two; sinuous vs varicose jet modes) |
| cube roots of a negative number (P262) | ch11 · §11.4 · C04 | −1 = e^{iπ} has cube roots −1 and ½(1 ± i√3); s³ = −a gives one real root and a complex-conjugate pair (the q, q* of the Bénard determinant) |
| neutral curve as the zero contour of the growth rate; critical point as its minimum (P263) | ch11 · §11.4 · C04 (before D10) | for each wavenumber find where the leading growth rate crosses zero (`brentq`, P108) → the neutral curve; minimise it over the wavenumber (`minimize_scalar`, P170) → the critical point; the same nested searches give Ra_c, Ta_c and Re_c |
| Helmholtz equation in the plane; planforms as sums of cosines (P264) | ch11 · §11.4 · C04 | linear theory fixes only ∣K∣: any f(x, y) with ∇_H²f = −K²f grows at the same rate — rolls cos Kx, squares cos Kx + cos Ky, hexagons (three wave vectors 120° apart) |
| quotient rule, and differentiating with respect to K² as the variable (P265) | ch11 · §11.4 · C05 | (f/g)′ = (f′g − fg′)/g²; a formula in K² only is minimised in x = K² (same point for K > 0) — used for (11.44) and to expose slip #3 |
| uniqueness of a linear boundary-value problem (P266) | ch11 · §11.5 · C06 | two unknowns with the same linear ODE, right side and boundary conditions are equal when the homogeneous problem has only the zero solution; (d²/dz² − K²)f = 0 with f = 0 at both walls has only f = 0 — hence T̂ ∝ ŝ in D13 |
| narrow-gap (small-curvature) approximation: expand in d/R and keep the leading order (P267) | ch11 · §11.6 · C07 | when d ≪ R₁ the curvature terms (1/R)d/dR, 1/R² are smaller than d²/dR² by d/R₁, (d/R₁)²: keep the leading order (as lubrication scaling, P188); x = (R − R₁)/d ∈ [0, 1], d/dR = (1/d)d/dx |
| singular point of an ODE (P269) | ch11 · §11.7 · C08 (before D17) | where the highest-derivative coefficient vanishes — U(z) = c in the Rayleigh and Taylor–Goldstein equations; never for a growing mode (c_i ≠ 0); a neutral mode has a kink or logarithm there (the critical layer); viscosity removes it |
| quadratic eigenvalue problem c²M₂ + cM₁ + M₀ and its companion linearisation (P268) | ch11 · §11.7 · C08 | (c²M₂ + cM₁ + M₀)v = 0 becomes an ordinary generalised eigenproblem twice as large with w = cv: [[0, I], [−M₀, −M₁]] (v, w) = c [[I, 0], [0, M₂]] (v, w); for 1 × 1 blocks it is the quadratic formula |
| mapping an infinite domain to [−1, 1] (P270) | ch11 · §11.7 · C08 | z = s tan(θξ), θ = arctan(z_max/s): crowds the Chebyshev points within ∣z∣ ≲ s and still reaches ±z_max, where the e^{−k∣z∣} tail is set to zero; derivatives by the chain rule; always check that the answer does not move when z_max (and N) double |
| integrating on a Chebyshev grid (Clenshaw–Curtis weights) (P271) | ch11 · §11.7 · C09 | ∫f dx ≈ Σw_j f(x_j) with weights exact for polynomials up to degree N — as accurate as the collocation itself; used for the integral identities and the energy budget (Simpson would be far worse on clustered points) |
| completing the square into a circle (x − a)² + y² ≤ R² (P272) | ch11 · §11.7 · C10 | x² + y² − 2ax + b ≤ 0 is the disc (x − a)² + y² ≤ a² − b, centre (a, 0), radius √(a² − b) — Howard's last step with x = c_r, y = c_i |
| the integral of an x-derivative of a periodic function over one period is zero (P273) | ch11 · §11.10 · C14 | ∫₀^λ ∂f/∂x dx = f(λ) − f(0) = 0: averaged over whole wavelengths every ∂(…)/∂x term disappears (the two ends of the control volume cancel) |
| Jacobian matrix of a nonlinear ODE system and the stability of its fixed points (P274) | ch11 · §11.14 · C15 | fixed point f(s*) = 0; nearby δ̇ = Jδ with J_ij = ∂f_i/∂s_j; stable if every eigenvalue has negative real part (P214); a complex pair crossing the imaginary axis is a Hopf bifurcation |
| Galerkin truncation: keep a few modes and project with orthogonality of sines (P275) | ch11 · §11.14 · C15 | unknown = a few fixed shapes × time-dependent amplitudes; multiply by each shape and integrate; orthogonality leaves one equation per amplitude and discards what does not fit (ch10 used hats, here sines and cosines) |
| purely imaginary roots of a cubic λ³ + a₂λ² + a₁λ + a₀: exactly when a₂a₁ = a₀ (P276) | ch11 · §11.14 · C15 | put λ = iω: ω² = a₁ and a₂a₁ = a₀; for positive coefficients all roots are stable when a₂a₁ > a₀ (Routh–Hurwitz); equality is the Hopf point |
| Lyapunov exponent: the slope of log(separation) against time (P277) | ch11 · §11.14 · C15 | ∣δ(t)∣ ≈ ∣δ₀∣e^{λt} ⇒ ln∣δ∣ is a line of slope λ; λ > 0 means chaos; time to reach an error Δ ≈ ln(Δ/δ₀)/λ (halving the initial error adds only ln 2/λ); fit only while the separation is small |
| iterated maps: fixed point, stability ∣f′(x*)∣ < 1, cobweb diagram (P279) | ch11 · §11.14 · C15 | x_{n+1} = f(x_n); x* = f(x*) attracts if ∣f′(x*)∣ < 1 (compare ∣G∣ ≤ 1 of Ch. 10); a cobweb goes up to y = f(x), across to y = x, and repeats |
| **ch11 — Python and numerics** | | |
| np.roots and np.lib.scimath.sqrt (P256) | ch11 · §11.3 · C02 | `np.roots(coeffs)` returns every root of a polynomial (highest power first, complex if needed); `np.lib.scimath.sqrt(-1.0)` returns `1j` where `np.sqrt` gives nan — what a discriminant that can turn negative needs |
| boundary-row replacement and the generalised non-symmetric eigenproblem scipy.linalg.eig(A, B) (P259) | ch11 · §11.4 · C03 | collocation gives Av = σBv; a boundary row in A states the condition and the same row of B is zero, so B is singular: `scipy.linalg.eig(A, B)` returns infinite eigenvalues (drop them) and round-off spurious ones (keep only eigenvalues that do not move when N grows); extends the symmetric `eigh(A, B)` of P247 |
| matplotlib 3-D line plots (P278) | ch11 · §11.14 · C15 | `ax = fig.add_subplot(projection="3d")`, `ax.plot(X, Y, Z)`, `ax.view_init(elev, azim)`; lighter than plotly (P41/P64) for a static page figure |
| **ch12 — maths and physics** (the term is the exact `nb.primer` title) | | |
| random variable, probability density and histogram (P280) | ch12 · §12.3 · C01 | a quantity that differs in every realization; the density says how often each value occurs (area = fraction of runs); a histogram is its estimate from samples |
| standard error of a mean (P281) | ch12 · §12.3 · C01 | the mean of N **independent** samples scatters by σ/√N; "within 5 standard errors" is the test for every sampled number; a correlated record holds only Δt/t_c independent samples |
| ergodicity (P282) | ch12 · §12.3 · C01 | for a stationary signal a long time average of one run equals the ensemble average over many — assumed every time a measurement is compared with the theory |
| Ornstein–Uhlenbeck signal (a Langevin equation) (P283) | ch12 · §12.3 · C01 | the chapter's test signal: each step keeps e^{−Δt/τ_c} of the old value and adds fresh Gaussian noise; variance σ², correlation e^{−\|τ\|/τ_c}, a cusp at zero lag (no Taylor microscale) |
| the sinc function and np.sinc (P285) | ch12 · §12.3 · C01 | sin(x)/x: 1 at 0, zero at multiples of π; `np.sinc(x)` is sin(πx)/(πx); appears whenever a wave is averaged over a window |
| covariance, the covariance matrix and its ellipse (P287) | ch12 · §12.4 · C02 | mean product of two zero-mean variables; the 2 × 2 matrix has variances on the diagonal; its eigenvectors are the axes of the scatter cloud (tilted cloud ⇔ non-zero covariance) |
| a quadratic that is never negative has discriminant ≤ 0 (P288) | ch12 · §12.4 · C02 (before D02) | aλ² + bλ + c ≥ 0 for every real λ ⇒ at most one real root ⇒ b² − 4ac ≤ 0: the whole proof of the Schwartz inequality |
| Fourier-transform pair: where the 2π sits (P290) | ch12 · §12.4 · C03 | a pair needs one factor 1/2π in total; the book puts it in the forward transform with angular frequency; for a real even function the transform is a cosine integral |
| change of variable in a density; units of a spectral density (P293) | ch12 · §12.4 · C03 | a density is "variance per unit of the axis": S(k)dk = S(ω)dω, so S(k) = U₀S(ω) for k = ω/U₀; S_e in m²/s, a wavenumber spectrum in m³/s² |
| counting the independent components of a symmetric tensor (P295) | ch12 · §12.5 · C04 | a symmetric 3 × 3 tensor has 6, a fully symmetric triple product 10 — the arithmetic of the closure problem |
| kinematic vs dynamic fluxes (P296) | ch12 · §12.5 · C04 | observers report H [W/m²] and τ₀ [Pa]; the equations carry $\overline{wT'}$ [K m/s] and u_*² [m²/s²]: H = ρC_p$\overline{wT'}$, τ₀ = ρu_*² |
| derivatives of functions of r = ∣r∣ (P297) | ch12 · §12.6 · C05 (before D07) | ∂r/∂r_j = r_j/r, ∂r_i/∂r_j = δ_ij, δ_jj = 3, r_jr_j = r² — everything needed to take the divergence of an isotropic tensor |
| inner, outer and overlap: two descriptions that must agree where both hold (P300) | ch12 · §12.9 · C10 (before D19) | near the wall one set of variables works, far away another; where both are valid the formulas must give the same answer (matching) — the separation argument with two variables |
| composite profile: inner law + outer correction (P301) | ch12 · §12.9 · C11 | add to the inner formula a correction that vanishes at the wall and reaches full size at the outer edge (`WT.coles_wake`) — usable across the whole layer |
| dividing two ODEs to eliminate time (P303) | ch12 · §12.10 · C13 (before D23) | dy/dt = f and dx/dt = g give dy/dx = f/g; if that is n·y/x the solution is the power law y ∝ xⁿ |
| the stability parameter ζ = z/L and integrating a flux–profile relation (P304) | ch12 · §12.11 · C15 | φ_m = (κz/u_*)dU/dz as a function of ζ alone (1 neutral, > 1 stable, < 1 unstable: Monin–Obukhov similarity); integrating dU/dz gives the logarithm minus a correction ψ_m |
| a double integral over a triangle (P305) | ch12 · §12.12 · C16 (before D26) | ∫₀^t dt′∫₀^{t′} g(τ)dτ covers 0 < τ < t′ < t; each τ is counted over a strip of length t − τ, so it equals ∫₀^t (t − τ)g(τ)dτ |
| proof by induction (P306) | ch12 · §12.12 · C16 | show the statement for n = 1, and that n − 1 implies n; then it holds for every n — used for the random walk $\overline{R_n^2}=nL^2$ |
| **ch12 — Python and numerics** | | |
| sliding-window average with np.cumsum (P284) | ch12 · §12.3 · C01 | the sum over a window is the difference of two cumulative sums — one pass over the data; edges within half a window are NaN |
| random-phase synthetic fields (np.fft.ifftn) (P286) | ch12 · conventions cell (front) | give every Fourier mode a chosen amplitude and a random phase and transform back: a field with a prescribed spectrum, solenoidal if projected — **kinematic, no cascade** |
| correlation by FFT with zero padding (P289) | ch12 · §12.4 · C02 | transform, multiply by the conjugate, transform back: all lags at once; pad with zeros or the end wraps onto the beginning; "unbiased" divides lag m by N − m (for a periodic record the wrap-round is the correct lag product) |
| the discrete Fourier transform as a Riemann sum; Parseval (P291) | ch12 · §12.4 · C03 | `np.fft.rfft(u)*dt` ≈ ∫u e^{−iωt}dt at ω_k = 2πk/(N dt); in the book's normalisation S = \|û\|²/(2πT) and the sum over both signs of ω times Δω is the variance |
| segment averaging (Welch) and leakage (P292) | ch12 · §12.4 · C03 | one raw periodogram is as noisy as its mean however long the record; averaging K segments cuts the noise by √K at the price of resolution; a taper (Hann) stops a line leaking into distant bins |
| a sympy averaging operator (P294) | ch12 · §12.5 · C04 (before D06) | teach sympy the rules of D01 — linear, leaves a mean alone, kills a single fluctuation, keeps a product of two as a new symbol — and let it average the equations (`ch12.rans_sympy`) |
| quadrature on a logarithmic grid (P298) | ch12 · §12.7 · C08 | a spectrum spans decades: with K = e^s, ∫E dK = ∫E K ds on an evenly spaced s = ln K |
| semi-log axes: a logarithm is a straight line (P299) | ch12 · §12.9 · C10 | on `ax.semilogx` the log law is a line rising 2.303/κ per decade (5.6 for κ = 0.41), and U⁺ = y⁺ is a curve |
| a nonlinear diffusion problem by Picard iteration on the eddy viscosity (P302) | ch12 · §12.10 · C12 | ν_T depends on the answer: freeze it from the last guess, solve the linear problem, update with under-relaxation, repeat until it stops changing (`shear_flow_eddy_viscosity_solve`) |

Reminders written in ch03 instead of new primers (point here): P13 log–log slope, P15 `assert np.allclose`, P16
animate, P17 slider_figure, P18 show_viz, P21/P22 finite differences, P25 partial derivative, P26 Taylor, P27 definite
integral, P29 lambda, P31 `solve_ivp`, P37 trapezoid, P38 product rule, P40 sympy, P41/P64 plotly 3-D, P46 `np.where`,
P49 chain rule, P62 `np.einsum`, P68 orders of smallness, P70 `arctan2`, P74 right-hand rule, P75 directional
derivative, P76 `meshgrid` and the grid layout, P77 broadcasting, P78 `streamplot`, P79 `expm`, P80 eigenvalues, P83
midpoint sums, P84 FTC, P85 mean-value theorem, P86 line integral round a loop.

Glosses in ch03 (one sentence where used): sin² + cos² = 1 to eliminate a parameter (D04), implicit differentiation
(D05), three-factor product rule (D11), Jacobi's formula det e^{Gt} = e^{t tr G} (C09, checked with `expm`), SVD for
finite-time ellipse axes (C12), Dirac delta as an infinitely concentrated finite total (N40, reminder of ch02 D21),
piecewise functions and continuity at a join (D18), maximum ⇒ derivative zero (D20), Lambert W (D20 cross-check),
Taylor in time F(x, t + Δt) ≈ F + Δt ∂F/∂t (D22), ∇·(Fu) = u·∇F + F∇·u (D24), plotly stream tubes (N15, code comment),
the ψ streamfunction (u = ∂ψ/∂y, v = −∂ψ/∂x, constant along streamlines, with a two-line sympy proof), a Python/sympy
idiom gloss in the setup cell.

Reminders written in ch02 instead of new primers (point here): P01 matplotlib, P04 f-strings, P05 stress, P06 linspace
(also `np.geomspace`), P09 Newton II, P10 random generator, P13 log–log slope (observed order), P14 tuple unpacking,
P15 `assert np.allclose`, P16 animate, P17 slider_figure, P18 show_viz, P21/P22 finite differences and np.gradient,
P23 dicts, P25 partial derivative (R02 recap of ∇), P26 Taylor, P27 definite integral, P29 lambda, P35 line integral
along a path, P37 trapezoid, P38 product rule, P40 sympy, P45 i² = −1, P47 live widgets, P49 chain rule, P53
determinants (det(AB) = det A det B was added as a gloss in D02 step 7: "volume-scale factors multiply"), P61 sympy
Matrix.

Glosses (one sentence where used, no demo) in ch02: index letters and ≡ (§2.1), bilinearity of the dot product (D01,
D02), continuity argument for det C = +1 (D02 step 8), det(AB) as volume-scale factors (D02 step 7), Rodrigues'
formula (cos θ I + sin θ [k×] + (1 − cos θ) k kᵀ, C02), "true for every n ⇒ coefficients agree" (D06 step 7), Mohr's
circle (C05/C06, taught in Ch. 4), eigen/principal frame forward gloss (C07, D18 step 8), `itertools.product`,
`np.count_nonzero`, `np.array_equal` (C08), ∂x_i/∂x_j = δ_ij (C10), einsum ellipsis `'ii...'` (C10), `ax.imshow`
(C08, C10), irrotational vortex u_θ = K/r (C11), `np.pad` (C12), Popoviciu's variance bound (D17 step 15), Riemann sum
(D26 step 7), Green's theorem as planar Stokes (D26), delta function "an infinitely concentrated source whose total
is finite" (D21), `try/except ValueError` (C13), climate hooks (Coriolis 2Ω × u, geostrophy, planetary vorticity)
in the front matter.

Reminders written in ch04 instead of new primers (point here; each "Tools from earlier chapters" 🔁 cell names the
chapter): P09 Newton II, P15 `assert np.allclose`, P16 animate, P17 slider_figure, P18 show_viz, P25 partial derivative,
P26/P98 Taylor, P27 definite integral, P28 net pressure force, P29 lambda, P31/P94 `solve_ivp`, P32 internal and kinetic
energy, P33 perfect gas, P37 trapezoid, P38 product rule, P40 sympy, P41/P64 plotly 3-D, P46 `np.where`, P49/P91 chain
rule, P50 Gibbs free energy (for Helmholtz f), P59 scaling of base units, P62 `np.einsum`, P68 orders of smallness,
P74 right-hand rule, P75 level sets, P76 `meshgrid`, P77 broadcasting, P78 contour/quiver/streamplot, P80 eigenvalues,
P84 FTC (variable upper limit for the pressure function), P85 mean-value theorem, P86 line integral round a loop, P87
`quad`/`dblquad`, P88 cylindrical unit vectors, P96 frames, P101 rigid-body velocity, P103 rotating frame (one instant;
ch04 P124 makes it time-dependent), P108 `brentq`, P110 signed swept volume.

Glosses in ch04 (one sentence where used, no demo): curl of a product ∇×(χ∇ψ) = ∇χ × ∇ψ + χ∇×∇ψ and ∇·∇× = 0,
∇×∇ = 0 (C03, D04), Newton's third law inside the drag-sign ⚠️ (N28), the side pressure force on a slowly widening tube
(D06), two-δ substitutions δ_imδ_jnS_mn = S_ij with `np.einsum` (D09, D10, D11), the Cartesian vector Laplacian
component by component (D12, D16), constant-acceleration kinematics s = ½at² (D17), the quotient rule
D(1/ρ)/Dt = −(1/ρ²)Dρ/Dt (D22), a·(a × b) = 0 (D25), "zero gradient everywhere ⇒ a function of t only" (D25, D26),
the barotropic pressure function ∫dp/ρ(p) (N87, D24), `np.cross` broadcasting over arrays of vectors (C09, C11),
`scipy.integrate.tplquad` (N07). Python idioms first met in ch04 and explained by their line comment only:
`sp.lambdify`, `np.ma.masked_where`, `contourpy`, plotly `make_subplots`, `sp.KroneckerDelta`, `sp.solve`/`.coeff`,
`np.outer` (the lesson reviewer accepted this; a ch03-style idiom gloss would be better).

Reminders written in ch05 instead of new primers (point here; each "Tools from earlier chapters" 🔁 cell names the
chapter): P29 lambda, numpy arrays, P04 f-strings, P92 parametric curves and arc length, P93 parallel vectors, P74
right-hand rule, P111 named results (NamedTuple), P76 `meshgrid` (with `indexing="ij"`, moved before first use in lesson
round 1), P83 midpoint sums, P15 `assert np.allclose`, P41/P64 plotly 3-D, log axes (P13), P01 matplotlib, P18 show_viz,
P88 cylindrical unit vectors, P124 rotating unit vector, P40 sympy, P108 `brentq`, P37 trapezoid and
`cumulative_trapezoid`, P116 torque and moment of inertia, P127 power of a force, P87 `quad`/`dblquad`, P17
slider_figure, P115 conservative force, P103 rotating frame, P70 `arctan2`/`unwrap`, P31/P94 `solve_ivp` (DOP853,
t_eval), P89 functions with parameters, P109 differentiation under the integral sign, P99 material line element, P98
Taylor, P128 chain rule for kinetic energy, P91 chain rule along a path, P95 RK4 by hand, P16 animate, P28 net pressure
force, P13 observed order, P66 `np.deg2rad` (moved before first use), P131 reduced gravity, P112 small-patch
localisation, `np.cross`/`np.linalg.norm`/`@` (added in lesson round 1), P121 Schwarz, P75 directional derivative, P123
erfc, P21/P22 finite differences and FTCS, P122 curl of a curl, P106 substitution, P10 random generator, P62
`np.einsum`, P72 parity of ε, P38 product rule, P65 projection and completeness, P79 `expm`, P44 separation of
variables, exp/log, P126 latitude and f, P125 product rule for a cross product, P22 `np.gradient` (added in lesson round
1).

Glosses in ch05 (one sentence where used, no demo unless noted): polar strain rate and polar stress divergence
(1/r²)∂(r²σ_rθ)/∂r (D03), **polar vector Laplacian (∇²u)_θ = u_θ″ + u_θ′/r − u_θ/r²** with a 2-line sympy check (added in
lesson round 1; the −u_θ/r² because e_θ turns), centre of mass of a non-uniform disc and disc I = ½MR² (D06), **moment
transfer M_G = M_O + (r_O − r_G) × F** (D06, round 1), quotient rule ∇(1/ρ) = −∇ρ/ρ² (D05 step 9, D15), smoothed step with
`np.tanh` (D07), **Taylor–Green flow** u = (sin x cos y, −cos x sin y)e^{−2νt}, the viscous twin of C03's cellular flow
(added in round 1), frozen-kernel far field with its 2a/distance estimate (D12), 1 + cot² = csc² (D13), cylindrical
Laplacian of ω_z(R) (D18), column mass Ah = const (D20), lever rule and circular motion (D21), mirror symmetry (D22),
erf = 1 − erfc (N22), limit of N filaments → a sheet (N44), comma notation and renaming the free index (D14, D15), Krasny
smoothing and the periodic sheet kernel (C14), "cat's eye" (C14), Gauss–Legendre pointer before its first use in C04
(round 1). Python idioms explained inline: `itertools.product`, `np.ndindex`, `np.roll`, `np.ptp`, `np.hypot`, symlog
axes, plotly `updatemenus`.

Reminders written in ch06 instead of new primers (the notebook's "Tools from earlier chapters" 🔁 cell and one line
where first used; point here): P25 partial derivatives, P26/P98 Taylor, P27 definite integrals, P28 net pressure force
−∮p n dA, P29 `lambda`, P37 trapezoid, P38 product rule (and the divergence form, ch04 P113), P40 sympy, P45 i² = −1 and
Euler's formula, P46 `np.where`, P49/P91 chain rule (incl. along a path), P57 `np.linalg.solve`, P63 `A @ x` (added in
lesson round 1), P67 `np.linalg.norm` (round 1), P68 orders of smallness, P70 `np.arctan2`, P71 Vieta, P75 directional
derivative and level sets, P76 `np.meshgrid` (with `indexing="ij"`, round 1), P77 broadcasting, P78 contour/streamline
plots (and `pcolormesh`, round 1), P80 eigenvalues and the spectral radius, P81 complex conjugate, P87 `quad`/`dblquad`,
P88/P105 polar and spherical unit vectors, P92 parametric curves and normals, P106 substitution, P108 `brentq`, P111
dataclasses and named results, P114 momentum flux, P121 Schwarz, P122 curl of a curl, P136 periodic trapezoid, P139
Poisson and Green's function, P140 gradient of 1/distance, P142 FFT, P143 Gauss–Legendre, P144 improper integrals,
P16 animate, P17 slider_figure, P18 show_viz, P31/P94 `solve_ivp`, P13 log–log slopes and `np.polyfit`, P15
`assert np.allclose`, P47 live widgets.

Glosses in ch06 (one sentence where used, no demo): dict comprehension `{k: f(v) for k, v in d.items()}` (added in
lesson round 1), logarithm rules ln√(1 + s) = ½ ln(1 + s) (D06), arcsin roots and the quadratic formula for r₊, r₋
(D10), odd/even functions for images (D12), derivative of arctan with a moving argument (D13), 1/i = −i (D14),
outward normal of a counterclockwise contour n = (e_x dy − e_y dx)/ds (D16), dz* = dx − i dy (D17),
e^{iθ} + e^{−iθ} = 2cos θ and foci c² = A² − B² (D20), cross products of cylindrical unit vectors (D24), flux through
a sphere (D25), the cot substitution d(cot α)/dα = −1/sin²α with the limits keeping their order (D26; corrected in
lesson round 1), "a function of x − x_s: ∂/∂x_s = −∇" and |A|² = A·A (D29), the qualitative separated-cylinder band
(C06), corner exponent n = π/α (C09), deflation (C07).

Reminders written in ch07 instead of new primers (point here; each "Tools from earlier chapters" 🔁 cell names the
chapter): P25 partial derivative, P26/P98 Taylor (one and several variables), P27 definite integral, P31/P94
`solve_ivp` (DOP853, tight tolerances for Stokes drift), P95 RK4 by hand, P40 sympy, P117 sympy expand/series, P44
linear second-order ODE (e^{λz} trial), P45 Euler's formula, P49/P91 chain rule (incl. along a path), P68 orders of
smallness, P75 level sets and the normal, P76 `meshgrid`, P77 broadcasting, P78 contour/quiver, P80 eigenvalues (the FD
sloshing problem), P81 complex conjugate, P87 `quad`/`dblquad`, P92 parametric curves, P104 ellipse geometry, P107
`np.expm1` and cancellation, P108 `brentq`, P114 momentum flux, P121 Schwarz (mixed partials), P127 power of a force,
P131 reduced gravity (with the ρ₂-vs-ρ₁ callout), P132 moving level set, P142 Fourier modes and the FFT, P149 the delta
(1-D version as a narrow tanh step), P151 integrals of sin/cos over a period, P153 the complex plane in numpy, P158 sympy
with complex symbols, P16 animate, P17 `slider_figure` (and `animate_figure`, the plotly time slider), P18 show_viz, P47
live widgets, P13 log–log slopes.

Glosses in ch07 (one sentence where used, no demo unless noted): overflow of cosh/sinh at large kH and the stable ratio
e^{kz}(1 + e^{−2k(z+H)})/(1 − e^{−2kH}) (`np.tanh` saturates safely), matching coefficients of cos(kx − ωt) (D06, D28),
inverting a monotonic function with a bracket from the two limits (brentq reminder + a 3-line Newton), time average ⟨·⟩
vs wavelength average (overbar) and ⟨cos²⟩ = ½, minimum of a function f′ = 0, f″ > 0 (minimise c² not c), derivative as
the limit of Δω/Δk (chord vs tangent), the complex-step derivative Im f(k + ih)/h (no cancellation error), quadratic
formula and the physical root (Vieta reminder, D25), sech² and its derivatives (the soliton), Lagrangian vs Eulerian
mean (ch03 recap sentence), complex 2 × 2 systems with sympy `I`, `solve`, `factor`, linearisation about a base state
(drop small × small, keep small × O(1)), angle of a vector from its components (cos θ = |k|/K, `np.arccos`,
`np.degrees`), gradient in wavenumber space ∇_K ω, the 1-D Dirac delta as the derivative of a narrow tanh step (P149
reminder), **WKB** (slowly varying wave trains, defined where first used), **Hamilton's equations** for rays
(dx/dt = ∂ℋ/∂p, dp/dt = −∂ℋ/∂x with ω as the Hamiltonian ℋ — written ℋ to avoid the depth H, lesson round 2), "take it
on trust for now (derived in C03–C04)" for the deep-water ω = √(gk) used in §7.1.

Reminders written in ch08 instead of new primers (point here; the "Tools from earlier chapters" 🔁 cells before first
use — erf and `np.trapezoid` before the entrance-length cell, `meshgrid`/contour and `observed_order` before the C06
figures, moved there in lesson round 1): dynamic vs kinematic viscosity (ch01 P05 stress), P23 dicts, P04 f-strings,
P14 tuple unpacking, P15 `assert np.allclose`, P10 random generator, P01 matplotlib, P13 log–log slopes and observed
order, P29 lambda, component-first arrays, unit normal/tangent (ch02), P108 `brentq`, P123 erfc (first look at erf),
P37 trapezoid, P25 partial derivative, P42/P167 separation of variables, P44 integrating an ODE twice, P84 FTC, P40 sympy
and `dsolve`, P57 `np.linalg.solve`, P21/P22 finite differences (tridiagonal from scratch), signed τ = μ du/dy (ch01),
P19 derivative as a slope, inequalities under division by a positive number, P26 Taylor, P87 `quad`, P17
`slider_figure`, P18 show_viz, P88 cylindrical coordinates, P36 ln R → −∞, P83/P163 area element 2πR dR, P122 curl of a
curl, centripetal acceleration (ch04), P68 orders of smallness, P22 `np.gradient`, P133 scaled variables and the chain
rule, P130 order-of-magnitude scaling, P109 differentiation under the integral sign, P132 kinematic condition at a
moving wall, P76/P78 `meshgrid`/contour/streamplot, P106 substitution, P170 `minimize_scalar`, P107 `np.log1p`/cancellation,
P47 live widgets, hydrostatics and a stress-free surface (ch01, ch04), P30 explicit stepping and its limit, P46 `np.where`,
P16 animate and `show_animation`, P49/P91 chain rule, P43 exponent rules, P17 `animate_figure`, P176 complex amplitudes,
P159 √i, P45 Euler's formula, P165 phase and time lag, ch06 spherical θ from the stream axis (with the downstream origin
flagged), P121 Schwarz, ch02 ε_ijk, P31/P94 `solve_ivp` (many tracers in one call), P35 line integral of a gradient,
P163 surface integrals on a sphere, P131 effective weight and buoyancy, P143 Gauss–Legendre, P117 series of 1 − e^{−s}.

Glosses in ch08 (one sentence where used, no demo unless noted): integrating an ODE twice with a 3-line sympy `dsolve`
demo, regularity at the axis (ln R → −∞), the area element 2πR dR, bookkeeping of a small parameter ε (multiply through
so one chosen term has coefficient 1), antiderivatives of (1 + αx/L)^{−n} and partial fractions, stress-free surface and
hydrostatic thin-layer pressure, order of a PDE vs number of conditions (two in y, one in t), the product rule read
backwards (ηF)′, complex trial Re{e^{iωt}f(y)} with √i = (1 + i)/√2, curl of a gradient is zero and curl commutes with
the Cartesian Laplacian, traction projection t_x = σ_rr cos θ − σ_rθ sin θ, line integral of a gradient to recover p,
asymptotic size of a term far away, the series 1 − e^{−s} = s − s²/2 + …; `np.diag`/`np.r_` for a hand-made
tridiagonal matrix; "the book skips this move; we add it" for inserted derivation steps.

Reminders written in ch09 instead of new primers (47 names in the 🔁 "Tools from earlier chapters" cells at the top of each
CORE block; point here): ν = μ/ρ and its diffusion time (P185), order-of-magnitude scaling and "~" (P130), orders of
smallness (P68), two-length anisotropic scaling (P188), scaled variables and the chain rule (P133), chain rule (P49/P91),
dominant balance (P198), partial derivative (P25), sympy `symbols`/`subs`/`diff`/`simplify` (P40), finite-difference residual
with `np.gradient` (P22), power laws and log–log plots (P13), dicts (P23), matplotlib (P01), f-strings (P04), tuple unpacking
(P14), lambda and functions as arguments (P29; plus the default-argument capture `lambda y_, U0_=U0_:` glossed where used),
`assert np.allclose` (P15), definite integral (P27), trapezoid (P37), `brentq` (P108), animate/`show_animation` (P16),
`slider_figure` (P17), `show_viz` (P18), product rule (P38), `solve_ivp` (P31/P94), `solve_bvp` (P196), RK4 by hand, erfc
(P123/P195), exponent rules (P43), live widgets (P47), Leibniz and differentiation under the integral (P109), Leibniz with a
moving limit (P189), FTC (P84), `cumulative_trapezoid` (P190), sign-change search with `np.sign` (P180), `np.interp` (P182),
the dp/dx sign convention (ch08), Gauss–Legendre (P143), eigenvalues and eigenvectors (P80), complex conjugate and the
complex plane in numpy, hyperbolic functions (P168), finite-difference Jacobian, the complex velocity of a point vortex
(ch05/ch06), exponent matching (P197), the free scale f → λf(λη) (P207), centripetal acceleration (ch04), broadcasting
`xs[None, :]` (P77).

Glosses in ch09 (one sentence where used, no demo unless noted): von Mises variables (ψ as the cross-stream coordinate: v
disappears and the grid follows streamlines), `CubicSpline` in a comment, Lamb's low-Re cylinder drag
C_D ≈ 8π/[Re(2.002 − ln Re)] (0.5 − γ_E + ln 8 = 2.002), the Goldstein singularity at separation (named), Hartree's
β = 2n/(n + 1), the axisymmetric form of continuity (1/R)∂(Ru_R)/∂R + ∂w/∂z = 0 in the teacup loop, "a bracket, not a number"
for inlet-sensitive marched separation, printed-vs-correct boxes for slips (and ⚠️ traps for statements that only look like
slips).

Reminders written in ch10 instead of new primers (63 "(reminder)" entries in `metadata.fluidpy.primers`; the 🔁 "Tools from
earlier chapters" cell lists them once and each block repeats what it needs; point here): first-order and multivariable
Taylor (P26/P98), finite differences and FTCS for diffusion (P21/P22), orders of smallness (P68), `np.linalg.solve` (P57),
exact rational weights (sympy `Rational`), power laws and log–log slopes / observed order (P13), diffusion spreading
s² = s₀² + 2Dt, neighbour slicing (P77), `linspace`/`logspace` (P06), matplotlib and log axes (P01), f-strings (P04), dicts
(P23), lambda (P29), `assert np.allclose` (P15), `slider_figure` (P17), `show_viz` (P18), partial derivative (P25), explicit
stepping and its limit (ch01 C12), boundary conditions, Crank–Nicolson and `solve_banded` (P193), implicit stepping (P192),
sympy `symbols`/`subs`/`series`/`removeO` (P40, P117), Fourier modes and FFT Poisson (P142), Euler's formula and complex numbers
in numpy (P45, P153), complex conjugate, seeded generator (P10), live widgets (P47), first-order wave equation and
characteristics (P174), separation of variables, `animate`/`show_animation` (P16), product rule (P38), integration by parts
(P218a), FTC (P84), `np.trapezoid` (P203), `solve_ivp` (P31/P94), `scipy.sparse` and `spsolve` (P161), substitution (P106),
chain rule (P49), Gauss–Legendre (P143), exponential trial for an ODE (P44), `np.expm1` (P107), ln and e⁻¹, e⁻², quadratic
formula, tridiagonal (Thomas) solve, RK4 by hand (P95), Schwarz's theorem (swap ∂t and ∂x), product rule for a divergence,
`expm` (P79), matrix multiplication (P63), Poisson equation (P139), `meshgrid` and `[j, i]` (P76), null space and rank (P58),
5-point Laplacian and Jacobi/Gauss–Seidel/SOR (P160), eigenvalues (P80), `np.interp` (P182), contour/quiver/streamplot (P78),
frames of reference (ch03), Newton's method (P152), `np.fft.rfft` for a dominant frequency, zero crossings with `np.sign`
(P180).

Glosses in ch10 (one sentence where used, no demo unless noted; several added by the lesson review): the Lamb wave (the
fastest external acoustic–gravity mode of the atmosphere, ≈ 310–320 m/s — the atmospheric CFL speed), Ch. 8's θ-method
(θ = 0, ½, 1) vs Glowinski's Θ-scheme, GLS = Galerkin/least-squares stabilisation, the numerical phase speed
c_num/u = −arg G/(Cθ), Gibbs-like overshoot of second-order schemes at a jump, `np.linalg.lstsq` (min-norm answer of a
singular system), `np.block` and `eigvalsh` (assembling and inspecting a saddle-point matrix), Euler's V − E + T = 1 for a
disc and 0 with one hole (P2 node counts), `tripcolor`/`tricontourf`/`triplot` (plotting on a triangle mesh), `re.findall`
(parsing a table header), Kovasznay flow (an exact steady NS solution used for FE orders), Hopf bifurcation (steady wake →
periodic shedding; → Ch. 11), the Arakawa C-grid (the MAC staggering in ocean/atmosphere models), "reduced speed of sound"
models, dynamics/physics splitting in GCMs, printed slips named "slip #1…#12" (not recap IDs), DEVIATION boxes.

Reminders written in ch11 instead of new primers (point here): P13 log–log, P15 `assert np.allclose`, P16 animate, P17
slider_figure, P18 show_viz, P25 partial derivative, P26/P98 Taylor, P31/P94 `solve_ivp`, P38 product rule, P40 sympy, P41/P64
plotly 3-D, P44 e^{mx} trial, P45 Euler's formula, P48 inequalities, P49 chain rule, P53/P56 determinants, P58 null space, P68
orders of smallness, P75 level sets, P80 eigenvalues, P81 complex conjugate, P95 RK4 by hand, P108 `brentq`, P121 commuting
operators, P129 completing the square, P133 scaled variables, P142 Fourier modes, P151 integrals of sines, P153 the complex
plane in numpy, P155 complex powers, P159 complex square roots and the quadratic formula, P162 collocation, P166 Taylor
transfer to the mean level, P168 cosh/sinh, P170 `minimize_scalar`, P177 operator elimination, P186 cylindrical Laplacian,
P188 anisotropic scaling, P203 `simpson`, P212 inflection point, P214 linear stability by eigenvalues, P215 sech, P218a
integration by parts, P233 arbitrary coefficients, P247 generalised symmetric eigenproblem, P249 Jacobian determinant, P252
caching; ch10 C07 (Galerkin with hats) and C04 (amplification factor G ↔ e^{σΔt}).

Glosses in ch11 (one sentence where used, no primer): **self-adjoint form** (pφ′)′ − qφ = 0, whose derivative term integrates
by parts into −∫p∣φ′∣² with no boundary term (C09, the only property D18 uses); **pitchfork** — one fixed point loses stability and two symmetric ones (C₊, C₋) appear (D25 step 7);
**eigenfunction** (P257); **capillary length** √(σ_s/(gΔρ)) (C02); **volume contraction rate** = trace of the Jacobian,
−(Pr + 1 + b) (C15); Reynolds stress −⟨uv⟩ (forward pointer in C02, defined in C14; → Ch. 12); Tollmien–Schlichting wave
(C13); exchange of stabilities and overstability (C01, C03); Clenshaw–Curtis vs Simpson; the compound-matrix and shooting
routes (named as independent checks, not taught); continuous spectrum vs discrete modes on unbounded profiles (N108, N109);
Routh–Hurwitz (inside P276); `np.genfromtxt`, `symlog` axes, `twinx`, `ListedColormap`, `Chebyshev.fit`, `np.polyfit`, pandas
tables (explained where they appear); Marangoni convection; thermohaline staircases; inertial instability and Rayleigh–Kuo
(pointers to Ch. 13); slips named "slip #1…#13" (S1–S13 in code).

Reminders written in ch12 instead of new primers (46 one-sentence reminders; point here): sinh (P168), seeded random
generators `np.random.default_rng` (P10), the Gaussian distribution, P13 log–log slopes and `np.polyfit`, P15
`assert np.allclose` and `np.isclose`, P16–P18 animate / slider_figure / show_viz, P25 partial derivative, the integral as
a limit of sums and the fundamental theorem of calculus, P40 sympy, P159 the quadratic formula, P106 substitution in an
integral or average, P261 even and odd functions, P26 Taylor series to second order, P38 product rule, P37 `np.trapezoid`
and P203 Simpson, P45 Euler's formula, P144 improper integrals, P142 `np.fft.rfft` and Fourier modes,
`scipy.integrate.quad`, P252 caching slow symbolic results, P119 isotropic tensors, P218a integration by parts, P22
`np.gradient`, P130 order-of-magnitude scaling, P57 `np.linalg.solve`, pandas tables, live widgets, P188 two-length
(thin-layer) scaling, P206 chain rule with a moving similarity variable, P189 Leibniz rule, P167 the separation argument,
P197 exponent matching, P194 Gaussian profile and integrals, P60 `fractions.Fraction`, P133 scaled variables, P49 chain
rule, P108 `brentq`, P190 `cumulative_trapezoid`, reading reference data files, P31/P94 `solve_ivp` and P95 an RK4 loop,
P48 inequalities under a sign change, P255 necessary vs sufficient (theorem vs observation), iterated integrals, P195 erf,
P109 differentiation under the integral sign.

Glosses in ch12 (one sentence where used, no primer): **stationary** and **homogeneous** statistics (independent of the time
or space origin); **realization** and **ensemble**; **solenoidal** (divergence-free) field; **Wiener–Khinchin** (the
spectrum is also the squared magnitude of the transform of the record; the periodic correlation is `ifft(|fft|²)/n`
exactly); **biased vs unbiased** correlation estimate and the triangle weight 1 − r/L (its kink gives a k⁻² tail);
**Lorentzian** spectrum of an exponential correlation; **frozen turbulence**; **closure problem**; **local isotropy**;
**cascade**, **inertial subrange**, **intermittency** (named); **entrainment** and **self-preservation**; **wall units**,
**buffer layer**, **wake function**, **indicator function** y⁺dU⁺/dy⁺; **van Driest damping**; **wall functions**;
**forced vs free convection**; **Monin–Obukhov similarity**; **Businger–Dyer** form (named "a commonly used form", source
unread); **Batchelor scale**; **Lagrangian** vs Eulerian integral scale; **ballistic** and **diffusive** regimes;
**Richardson's 4/3 law**; `np.fft.fft`/`ifft` (complex, full), `np.fft.ifft2`/`ifftn`, `itertools.product`, `ax.semilogx`,
`scipy.signal.welch` (cross-check only); slips named "slip #1…#18"; the library string "Γ < Γa" read as Γ_met < Γ_d.
