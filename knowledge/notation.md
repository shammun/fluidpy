# Notation register — symbols used in fluidpy code, in our own words

Filled chapter by chapter by the knowledge-keeper from the symbols the code actually uses (not a copy of the book's
nomenclature list). A symbol whose meaning changes — between chapters or inside one chapter — gets a ⚠️ and one row
per meaning. "ch01 →" means introduced in ch01 and still valid.

## Sign traps first

**⚠️ Two conventions for the lapse rate.** The physics is identical; the sign and the direction of the stability
inequality flip.

| Convention | Definition | Γ_a (dry air) | Stable when | Used by |
|---|---|---|---|---|
| Kundu (ours) | Γ ≡ dT/dz | ≈ −9.8 K/km (−g/C_p = −9.7607 K/km) | Γ > Γ_a (dθ/dz > 0) | this book, all fluidpy code (`dT_dz` arguments) and explainers |
| Meteorology (standard) | Γ ≡ −dT/dz | ≈ +9.8 K/km | Γ < Γ_a | most atmospheric-science texts, incl. the lapse-rate-feedback literature |

Converting: Γ_met = −Γ_Kundu. When reading a meteorology paper, negate its Γ before passing it to fluidpy, and flip the
inequality. `lapse_rate_convention(dT_dz, "meteorology")` returns −dT/dz; `lapse_rate_stability(dT_dz,
convention=…)` returns the verdict and the criterion written with numbers in either convention. The notebook's C54
slider and the `parcel_stability` badge always show both.

**⚠️ Three pressures with a subscript.** `p0` = pressure at z = 0 (Eq. 1.9, `core.statics`); `p_ref` = `P_REF` =
1.0e5 Pa, the reference pressure of θ and ρ_θ (Eqs. 1.31, 1.33; the book writes p_o and takes it as the surface
pressure, meteorology and our code use exactly 1000 hPa); `p_o` in the Laplace jump (1.5) is the pressure outside a
curved interface. `P_ATM` = 101 325 Pa is the standard atmosphere.

**⚠️ Heat: two quantities named q.** In §1.8 `q` is heat added per unit mass [J/kg]; in Fourier's law (1.2) **q** is a
heat-flux vector [W/m²]. Code: `process_heat_work` returns q [J/kg]; `fourier_heat_flux` returns a flux.

**⚠️ Passive vs active rotation (ch02).** The book's direction-cosine matrix is C_ij = e_i·e'_j (row = old axis,
column = new axis; its columns are the new unit vectors in old components). Components transform **passively**:
x' = Cᵀx (2.5), x = Cx' (2.7), τ' = CᵀτC (2.12). For axes turned by +θ, C *equals* the active rotation matrix R(θ) of
Wikipedia and `scipy.spatial.transform.Rotation.as_matrix()`, but it is applied as its transpose: scipy's R rotates
the arrow, Cᵀ re-describes a fixed arrow in turned axes. `rotation_matrix_2d(θ)` / `rotation_matrix_3d(axis, θ)` return
C; `transform_vector(x, C)` computes Cᵀx. To use a scipy rotation as the book's C, pass `R` as C (not `R.T`).

**⚠️ Which index gets contracted (ch02).**

| Operation | Formula | Contracted index | Code |
|---|---|---|---|
| Cauchy traction | f_i = τ_ji n_j = (n·τ)_i (2.15) | **first** | `traction(tau, n)` (never symmetrises) |
| Tensor divergence | (∇·τ)_i = ∂τ_ij/∂x_j (§2.9) | **second** | `tensor_divergence(T, h, index=1)` (default) |
| A·u vs Aᵀ·u | A_ij u_j vs A_ji u_j (2.14) | second vs first | `dot_tensor_vector(A, u, index=1 / 0)` |

The two agree only for symmetric τ (proved in Ch. 4). Tests discriminate both with a non-symmetric τ.

**⚠️ Two double dots (ch02).** Book A:B = A_ij B_ji = tr(AB). Frobenius A_ij B_ij = tr(ABᵀ). They coincide when either
operand is symmetric. `double_dot(A, B, convention="book" | "frobenius")`: always pass the convention explicitly.
E3's term bars (S:S + 2S:A + A:A = G:G) use Frobenius.

**⚠️ Two gradients of a vector (ch02).** `vector_gradient(u)[i, j] = ∂u_i/∂x_j` (velocity gradient G; Ch. 3 builds
S and R from it). `integral_gradient` and `gauss_gradient_box` return [i, j] = ∂Q_j/∂x_i, the natural order of (2.31),
which is the **transpose**.

**⚠️ Rotation tensor vs antisymmetric part: the factor-2 trap (ch02; review M1).**

| Name | Definition | Its vector (via ω_k = −½ ε_ijk R_ij) | Code |
|---|---|---|---|
| R (book's rotation tensor, §2.10; Ch. 3 (3.15), (3.17)) | R = G − Gᵀ, R_ij = ∂u_i/∂x_j − ∂u_j/∂x_i | **ω = ∇×u (vorticity)** | `rotation_tensor(G)` |
| A (antisymmetric part) | A = ½(G − Gᵀ) = ½R | **½∇×u (spin/angular velocity of the element)** | `antisymmetric_part(G)` |

G = S + ½R = S + A. R_ij = −ε_ijk ω_k and R·x = ω × x. **Never call the vector of A "ω".** For u = b × x:
vector(R) = 2b, vector(A) = b.

**⚠️ Stokes orientation (ch02).** n is the chosen normal of the open surface A. n_c is the in-surface normal to the
rim C that points **into A** (Fig. 2.10; some of our early text said "outward", which is wrong). t = n_c × n is the rim
tangent, **counterclockwise seen from the tip of n**; (n_c, n, t) is right-handed. For a disc in the x–y plane with
n = +e₃: n_c = −r̂ and t = +θ̂. Flipping n flips both sides of (2.34). Closed surfaces (Gauss) use the outward n.

**⚠️ Γ, γ: four meanings for one Greek letter (ch01–ch03; review Should-fix 3).** Say which one in every docstring,
slider label and callout.

| Symbol | Meaning | SI unit | Where | Code |
|---|---|---|---|---|
| Γ (ch01) | environmental lapse rate dT/dz (Kundu) | K/m | §1.10 | `dT_dz`, `Gamma` in `core.stratification` |
| Γ (preset rate) | velocity gradient du₁/dx₂ of the `simple_shear` preset u = (Γx₂, 0) (**= ch03's γ = 2S₁₂**); for the other presets just the rate scale (solid body: Γ = ω₀, so (∇×u)₃ = 2Γ) | 1/s | ch02 E3, ch03 E4/E5 presets | `velocity_gradient_preset(name, Gamma)` |
| Γ (ch02 Ex. 2.4) | the strain-rate element S₁₂ itself (half of the preset's Γ); the book's "2S₁₂ = Γ" line is inconsistent with its own matrix | 1/s | §2.11 | `example_2_4(Gamma)` |
| γ (ch03 §3.5) | shear rate du₁/dx₂ of the parallel shear flow; S₁₂ = γ/2, ω₃ = −γ | 1/s | §3.5, D14 | `parallel_shear_kinematics(gamma)`; E5 "Rate k" = γ in shear mode |
| **Γ (ch03 →)** | **circulation** ∮u·ds = ∫ω·n dA (3.18) | m²/s | **already in ch03** (§3.4–3.5, `core/vortices.py`), not only from ch05 | `Gamma` in `rankine_vortex`, `gaussian_vortex`, `vortex_profile`; `circulation`, `circulation_circle` |

ch02 wrote Γ_circ for circulation to avoid the clash; from ch03 on Γ is the circulation whenever a vortex or a loop is on
screen, and the preset rate should be called γ or k in new text.

**⚠️ Γ / γ in ch05 — five meanings in the project, and the book adds a per-length one.** The table above plus:

| Symbol | Meaning | SI unit | Where | Code |
|---|---|---|---|---|
| Γ (ch05) | circulation of a loop / strength of a vortex tube, filament, line or point vortex; **counterclockwise positive** | m²/s | §5.1–5.7 | `Gamma` everywhere in `core.vorticity`, `core.biot_savart`, `ch05` |
| Γ_a (ch05) | absolute circulation Γ + 2∫Ω·n dA = Γ + 2Ω·A_vec (5.33) | m²/s | §5.6 | `absolute_circulation(u, pts, Omega)` → (Γ, Γ_a) |
| γ (ch05, book writes Γ) | **vortex-sheet strength per unit length** = jump in tangential velocity, γ = u₂ − u₁ = u_below − u_above (counterclockwise circuit) | m/s | §5.8, Fig. 5.16 | `gamma` in `vortex_sheet_strength(u_above, u_below, convention="ccw")`, `diffusing_vortex_sheet`, `sheet_rollup` |

Rule: `Gamma` [m²/s] is always a circulation; `gamma` [m/s] in ch05 is always a sheet strength (ch03's γ was a shear rate,
ch04's the third isotropic coefficient, ch01's C_p/C_v). The Fig. 5.16 caption's u₁ − u₂ is the clockwise count
(`convention="caption"`, −2 for the text's +2).

**⚠️ ω is the vorticity in (5.1); the tank turns at ω/2 (ch05).** u_θ = ωr/2 with ω = ∇×u. A tank turning at Ω = 1 rad/s
has ω = 2 s⁻¹; Fig. 5.2's "2ω" label is a slip. ch03's `solid_body_rotation(r, omega0)` takes the *rate* ω₀;
`ch05.solid_body_from_vorticity(r, omega)` converts (passing ω as the rate makes the pressure 4× too big, a planted
variant caught by 5 tests). Elsewhere ω still means angular frequency in ch03 Ex. 3.1 and Ch. 7.

**⚠️ Baroclinic vector: order and angle (ch05 (5.28), E4).** The source is (1/ρ²)∇ρ × ∇p — ∇ρ first (the reversed order
fails 5 tests). E4 and the notebook write ∇p = (0, −ρ₀g) (hydrostatic) and **∇ρ = |∇ρ|(sin θ, −cos θ)**, θ = the tilt of the
isopycnals from the isobars: θ = 0 stable, no torque; (∇ρ × ∇p)_z = −|∇ρ|ρ₀g sin θ, so 0 < θ < 180° spins **clockwise**
(−0.0981 s⁻² at θ = 90°, |∇ρ| = 10 kg/m⁴, ρ₀ = 1000) and θ = 270° counterclockwise. Lock exchange with the heavy fluid on
the left: +2.42 s⁻² (counterclockwise). The design first had this angle reversed; the builders caught it by computing.

**⚠️ Natural coordinates: e_n = −(Frenet N) (ch05 (5.31), Fig. 5.9).** e_s along ω; the book's e_n points **away** from
the centre of curvature (the Frenet normal N points toward it); e_m = e_s × e_n. For a helix of radius a and pitch 2πc,
e_n is radially outward, κ = a/(a² + c²), τ = c/(a² + c²) (τ = torsion here, not stress). `ch05.helix_frame`,
`frenet_frame` return the book's frame; the Frenet variant fails 2 tests.

**⚠️ (5.14) sign and ∇′ (ch05).** The book prints u = −(1/4π)∫(∇′×ω)/|x − x′| d³x′; the right factor is **+1/(4π)** (∇²u =
−∇×ω, G = −1/(4π|x − x′|)). Its integrand rewrite flips a second sign, so (5.16) is right. ∇′(1/|x − x′|) =
+(x − x′)/|x − x′|³ (derivative with respect to the source point = minus the derivative with respect to x).
`core.biot_savart.velocity_from_curl_omega(..., sign=+1)`; `sign=-1` reproduces the printed form (reverses the swirl).

**⚠️ Elliptic integrals take the parameter m = k² (scipy).** `scipy.special.ellipk(m)`, `ellipe(m)`; the modulus k is the
square root. `ring_ring_velocity` passes m (the modulus variant misses by > 5 %). Explainers mirror it in `ellipKE`
(AGM). K(m) → ∞ as m → 1 (a point on the filament).

**⚠️ Burgers' α (ch05 N21, D18).** The book's u_z = αz, u_R = −½αR; Wikipedia writes the same flow with α_w = α/2. Core
radius √(4ν/α) in the book's α. `core.vortices.burgers_vortex(R, z, Gamma, alpha, nu)` uses the book's.

**⚠️ Book slips taught corrected (ch05, analysis §9).** (5.14) −1/(4π) → +1/(4π), with a compensating slip in the
integrand rewrite (so (5.16) is right) · (5.27) silently drops u_{j,j}(ω_n + 2Ω_n) (zero by (5.19)) · "Π" for Φ in (5.26) ·
the Lamb step writes ∇(u·u) for ∇(½u·u) (harmless under the curl) · Fig. 5.2's "2ω" (the tank turns at ω/2) · "ρ and p
single-valued" is not why ∮dp/ρ = 0 — barotropy ρ = ρ(p) is · "irrotational C ⇒ no viscous term" holds only for
incompressible constant-μ flow · Fig. 5.11's G is the centre of vorticity, not a stagnation point · Fig. 5.16 caption
u₁ − u₂ vs text u₂ − u₁ · "hyperboloids of the second degree" (Fig. 5.3) are cubic surfaces (c − z)r² = const · exercise
pointers "5.8" after (5.14) → 5.9 and "5.11" after (5.33) → 5.10.

**⚠️ Rotating tank is not geostrophic (ch05 lesson review).** Water at rest in a rotating tank has no Coriolis force
(u′ = 0); its paraboloid is gravity + centrifugal folded into an effective gravity (a geopotential surface, the reason
for the Earth's equatorial bulge). Geostrophy is Coriolis against a pressure gradient and needs motion relative to the
frame.

**⚠️ R has no ½; ω = 2 × spin (ch03 (3.13), (3.15)).** R = G − Gᵀ, its vector ω = ∇×u (vorticity), G = S + ½R. A small
element turns rigidly at **½ω** (`element_rotation_rate` = ½ω₃); a paddle wheel turns at ½ω; solid-body rotation at
ω₀ has ω = 2ω₀. Vector of the antisymmetric part A = ½R is the spin ½ω, never "ω".

**⚠️ Galilean frame sign (3.9) (ch03).** O′ moves at constant U with axes parallel to O: x = x′ + Ut + x′_o, t = t′,
u = U + u′ (so u′ = u − U). Towed cylinder: in the body frame the far fluid moves at +U e_x (steady); in the fluid frame
the cylinder moves at −U e_x (unsteady). `cylinder_flow(..., U_frame)`: observer moves at −U_frame e_x relative to the
far fluid (0 = fluid frame, U = body frame). The sum ∂u/∂t + (u·∇)u is frame-free; the local/advective split is not.
**Primes changed meaning**: ch02 x′ = Cᵀx is a *rotated* frame; ch03 x′ is a *translating* frame.

**⚠️ Rotating vs Galilean (ch03 D13).** An observer rotating at Ω about z measures ω′_z = ω_z − 2Ω (not ω − Ω); the
co-rotating frame Ω = ω_z/2 sees no spin. On the Earth: relative vorticity ζ vs absolute ζ + f (Ch. 4 §4.7, Ch. 13).

**⚠️ Signed wall term in the Reynolds transport theorem (ch03 (3.35)).** ∫_{A*} F b·n dA with the **outward** n and
the **sign** of b·n kept: an advancing wall (b·n > 0) sweeps F in, a retreating wall (b·n < 0) sweeps it out, a wall
sliding tangentially sweeps nothing (b ≠ 0 but b·n = 0 — Ex. 3.2's base). `|b·n|` is wrong (5 tests fail). In Leibniz
(3.30) the lower-limit term ȧF(a) is **subtracted** (`leibniz_terms(...).lower`).

**⚠️ (3.6) needs F (book typo).** The book prints DF/Dt = ∂F/∂t + |u| ∂/∂s; the correct streamwise form is
DF/Dt = ∂F/∂t + |u| ∂F/∂s (the printed one has units 1/s, not [F]/s). Code: `streamwise_derivative`.

**⚠️ Explainer "Rate k" (ch03 E5).** One slider drives flows whose rate means different things: shear → γ = du₁/dx₂,
solid body → ω₀, pure strain → S₁₁ = −S₂₂ = k/2, other presets → G = k × a fixed matrix. The slider is labelled
"Rate k" and a meaning line (`kMeaning`) says which, next to G. Later explainers with mode-dependent rates copy this.

**⚠️ The prime has four meanings — and ch04 uses three of them (ch04 §4.1 notation callout).**

| Where | x′, u′, p′ mean | Code names (one name per meaning) |
|---|---|---|
| ch02 (2.5) | components in **rotated axes** (x′ = Cᵀx) | `xp` in `core.tensors` |
| ch03 (3.9) | **translating (Galilean) frame** O′ moving at constant U | `xp`, `tp`, `x0p` in `galilean_transform` |
| ch04 §4.7 (4.42)–(4.45) | **noninertial frame** translating at U(t) and rotating at Ω(t) | `u_prime`, `x_prime`, `a_prime` in `core.rotating`; notebook `u_rot` |
| ch04 §4.9 (4.84)–(4.89), §4.11 | **perturbation** from the hydrostatic state, p′ = p − p_s(z), ρ′ = ρ − ρ_s(z) | `perturbation_fields` → `p_pert`, `rho_pert`; `buoyancy(rho_pert, rho0)` |
| ch04 (4.67) | a **dummy integration variable** in ∫_{p_o}^{p} dp′/ρ(p′) | inside `pressure_function` only |

**⚠️ Coriolis: acceleration term vs apparent force (ch04 (4.43) vs (4.45)).** Seen from a frame rotating at Ω, the
particle's inertial acceleration is D′u′/Dt + dU/dt + **2Ω × u′** + Ω̇ × x′ + **Ω × (Ω × x′)** (4.43) — the
*acceleration* side, "+2Ω × u′" and the *centripetal* Ω × (Ω × x′), pointing towards the axis. Moved to the right of
Navier–Stokes they become apparent body forces per mass: **−2Ω × u′** (Coriolis force) and **−Ω × (Ω × x′) = +Ω²R e_R**
(centrifugal, away from the axis), with −dU/dt and −Ω̇ × x′ (4.45). The book's prose "the Coriolis acceleration 2Ω × u
deflects a particle to the right" describes −2Ω × u. Code: `frame_acceleration_terms` returns the (4.43) acceleration
terms, `apparent_body_forces` / `coriolis_force` / `centrifugal_acceleration` the (4.45) forces (`coriolis_acceleration`
= +2Ω × u′). In the NH (Ω_z > 0) the force deflects moving parcels to the right. `rotating_frame_coriolis` labels every
bar with its signed form ("−2Ω×u′ Coriolis") and a side toggle "forces (−) / accelerations (+)" — pinned by exact-text
selftest rows. The factor 2 is two separate Ω × u′ (turning basis + d(Ω × x′)/dt, D15).

**⚠️ Rossby number Ro = U/(2Ωl) in ch04 (forward pointer).** `core.similarity.rossby_number(U, Omega, l, factor=2.0)`
= advective U²/l over the Coriolis term 2ΩU of (4.45). Ch. 13 writes Ro = U/(fL) with f = 2Ω sin φ
(`core.rotating.coriolis_parameter`), which equals U/(2ΩL) only at the pole — at 45° f = √2 Ω. Some texts use U/(ΩL);
pass `factor=1` for that. 10 m/s over 1000 km: Ro = 0.0686 (ours).

**⚠️ Tensor divergence contracts the FIRST index in Cauchy's equation (ch04 (4.20b)–(4.24)).** f_j = n_iτ_ij, so the
net surface force per volume is ∂τ_ij/∂x_i: `stress_divergence(…)` / `tensor_divergence(T, h, index=0)` — ch02's
default `index=1` is the other one. The book's prose after (4.24) and in §4.8 writes ∂τ_ij/∂x_j; harmless only because τ
is symmetric (4.25), which is proved after (4.24). A non-symmetric test τ (τ₂₁ = x₁ at rest) discriminates.

**⚠️ 2-D stream function sign (ch04 §4.3).** χ = −z ⇒ ρu = ∂ψ/∂y, ρv = −∂ψ/∂x (ρu = −e_z × ∇ψ); ψ increases to the
left of the flow; flux between two streamlines = ψ₂ − ψ₁ per unit depth (ρ = 1: m²/s). Many GFD texts use u = −∂ψ/∂y,
v = ∂ψ/∂x (ψ = geostrophic pressure/f): negate ψ when importing. Axisymmetric (Stokes) ψ: χ = −φ ⇒ ρu_R = −(1/R)∂ψ/∂z,
ρu_z = (1/R)∂ψ/∂R; flux through a ring = 2πΔψ. Uniform stream: ψ = Uy (2-D), ψ = ½UR² (axisymmetric).

**⚠️ Mean vs thermodynamic pressure; "deviatoric" σ (ch04 (4.27), (4.33)–(4.37)).** p = thermodynamic pressure (from
the equation of state); p̄ = −⅓τ_ii = mean normal stress; p − p̄ = μ_v∇·u (4.34), μ_v = λ + ⅔μ (bulk viscosity). The
Stokes assumption is **μ_v = 0, not λ = 0** (the λ = 0 variant fails 5 tests). σ_ij is called "deviatoric" but is
traceless only if μ_v = 0 or ∇·u = 0 (tr σ = 3μ_v∇·u). For plane (2 × 2) stress, p̄ needs τ₃₃ = −p + (μ_v − ⅔μ)∇·u:
`mean_pressure(tau, tau33=…)` raises without it (review Must-fix 2). Incompressible p is mechanical only, defined up to
a constant.

**⚠️ Two defaults for g (ch04).** `fluidpy.core.*` default to `G0` = 9.80665 m/s² (standard gravity, ch01); the ch04
chapter module's functions default to `ch04.G` = 9.81 (the book's value). The notebook passes `g=9.81` explicitly
when it calls a core function next to a ch04 one (lesson review Should-fix 3). `ch04.plane_poiseuille(y, G=…)` uses G
for −dp/dx, **not** `ch04.G`. `effective_gravity` uses a toy uniform g_n = 9.80 (so |g_e| = 9.783 m/s² at 45°, not
standard gravity). Air C_p: `core.bernoulli.CP_BOOK` = 1004.5 (book) vs ch01's derived `CP_AIR` = 1004.7.

**⚠️ Book slips taught corrected (ch04, analysis §9).** (4.15) ends with a spurious "= 0" (it is an identity) ·
(4.51) prints dA inside volume integrals (dV) · (4.74) gauge must be φ − ∫B dt′ (the printed + doubles B,
`gauge_absorbed_bracket`) · after (4.63) "μ, κ, k > 0" means μ, μ_v, k ≥ 0 (κ is a 4th-edition leftover) · §4.10 curve
C is ζ = x²/2R₁ **+** y²/2R₂ · Ex. 4.7 drops a minus in the separated ODE and writes "η = h/δ" for γ = h/δ · Ex. 4.2's
"(½)ρU² + gz + p/ρ" is dimensionally inconsistent (4.19 is ½U² + gz + p/ρ) · Cauchy prose ∂τ_ij/∂x_j (first index is
right) · "Section 3.6" for (3.14) is §3.4 · (4.100) ↔ (4.101) bracket reference; "(4.106), (4.107)" before (4.114) mean
(4.109), (4.112) · Ex. 4.8 total rounded after rounding the wave drag (+1.1e-3 vs ours) · Fig. 4.9 "budge" · Earth bulge truncated
(WGS-84 42.77 km) · (4.45) is derived from the incompressible constant-μ (4.39b) and needs U, Ω uniform in space ·
(4.9) is displayed twice with different content (§4.2 Dρ/Dt = 0; §4.11 the general identity with c²).

**⚠️ Grid layout: two orders (ch02, project-wide).** Arrays are indexed `[k, j, i]` = (z, y, x) in 3-D and `[j, i]` in
2-D, with x on the **last** axis (`np.meshgrid(z, y, x, indexing="ij")`; 2-D `indexing="xy"`). Vector components sit on
axis 0 (`u[c, k, j, i]`). But `Grid.h`, direction numbers d, `gradient` components and field components are ordered
**(x, y, z)**. Direction d lives on array axis ndim − 1 − d (`axis_of_direction`). Values are at nodes; periodic
grids drop the duplicate end node (h = L/n).

**⚠️ Γ clockwise vs counterclockwise inside ONE chapter (ch06).** The project's Γ (ch03 →, ch05) is counterclockwise
positive, and so are the book's (6.6), (6.8) ψ = −(Γ/2π) ln r and (6.47). But the book's (6.36)–(6.40), (6.52),
(6.61)–(6.62), (6.68) and Example 6.1 ("a vortex of strength −Γ") use a **clockwise** Γ: the flow's circulation is −Γ,
which is why L = +ρUΓ comes out positive (Wikipedia's Kutta–Joukowski takes its contour clockwise for the same reason).

| Where | Γ means | Lift for U = 10 m/s, ρ = 1.2 kg/m³, Γ = 2 m²/s | Code |
|---|---|---|---|
| (6.6), (6.8), (6.47), `Vortex`, ch05 | counterclockwise circulation | −24 N/m (a ccw vortex in a stream from the left pushes down) | `Vortex(Gamma)`, `Gamma_ccw=` |
| (6.36)–(6.40), (6.52), (6.61)–(6.62), (6.68), Ex. 6.1 | clockwise strength (circulation −Γ) | **+24 N/m** = ρUΓ | `Gamma_cw=` |

Rule: functions that follow the book take the keyword `Gamma_cw=` or `Gamma_ccw=` (never a bare `Gamma`);
`core.potential.gamma_ccw_from(Gamma_cw=…)` converts; loop circulation of the book's cylinder is −Γ_cw on every radius.
Explainers print both ("Γ_cw 2 ⇔ Γ_ccw −2 m²/s"); the planted "book Γ read as counterclockwise" variant fails the polar
velocity test.

**⚠️ The 2-D doublet vector points from the sink to the source (ch06 (6.28)–(6.29), (6.49)).** d = Σx_i m_i, so a
source at +ε and a sink at −ε give d = +2mε e_x and φ = −d·x/2πr². The cylinder needs **d = −2πUa² e_x** (pointing
upstream, "opposing the stream"); (6.49)'s scalar d is a dipole −d e_x (`Doublet.from_book_scalar(d)`). 3-D: (6.88)'s
dipole is −d e_z, the sphere d = 2πa³U; the moving sphere's (6.97) d(t) = +2πa³u_s (it follows the motion). A reversed
dipole misses the pair's far field by > 0.1 (test).

**⚠️ Three angle origins in ch06.** (1) Polar formulas (6.21)–(6.22), (6.33)–(6.39), (6.89)–(6.91): θ from **+x
(downstream)**, so the upstream stagnation point is θ = π and C_p = 1 − 4 sin²θ is symmetric anyway. (2) Fig. 6.10-style
real-vs-ideal comparisons: angle from the **upstream stagnation point** (= π − θ; `separated_cp_band(theta_front_deg)`).
(3) (6.106): θ_s from the **sphere's velocity** (θ_s = π − θ of (6.91)). Fig. 6.8's zero-C_p angle 113.2° is from +x at
the source. Always name the origin next to an angle.

**⚠️ w = φ + iψ (§6.4–6.6) vs w = z-velocity (§6.9); dw/dz = u − iv.** In §6.4–6.6 w(z) is the complex potential and
its derivative is the **conjugate** of the velocity (u + iv = conj(dw/dz)); in §6.9 (and Ch. 4, 13) w is the z-component
of velocity. `Flow.w(z)` is always the potential; velocities come from `Flow.velocity(x, y)` or `.complex_velocity`.

**⚠️ ζ is the Zhukhovsky (circle) plane in ch06.** z = ζ + b²/ζ (6.65) maps |ζ| = b to the slit [−2b, 2b] and
|ζ| = a > b to an ellipse; the inverse has two roots with ζ₁ζ₂ = b² and the physical one is **outside** |ζ| = b:
ζ = ½[z + √(z − 2b)√(z + 2b)] (`joukowski_inverse(branch="outside")`). numpy's principal ½[z + √(z² − 4b²)] returns
the inside root at every Re z < 0 point off the slit (|ζ| = 0.37 instead of 2.70 at z = −3 ± 0.5i, b = 1). ζ was the
relative vorticity in ch05 and a parcel displacement in ch01.

**⚠️ Stokes stream function [m³/s] vs plane ψ [m²/s] (ch06 §6.8).** Axisymmetric: u_R = −(1/R)∂ψ/∂z,
u_z = (1/R)∂ψ/∂R (6.75); spherical u_r = (1/(r² sin θ))∂ψ/∂θ, u_θ = −(1/(r sin θ))∂ψ/∂r (6.83); the volume flux between
two stream surfaces is **2πΔψ** (6.78), not Δψ; the field equation (6.77) is not the Laplacian (so no complex
variables). In §6.8 the symmetry axis z is **horizontal, along the stream** (z was "up" in ch01–ch05). ξ is the axial
coordinate along a line sink in §6.8 but the vector x − x_s in §6.9.

**⚠️ (6.82) is correct as printed (ch06 review M1).** The book's (1/r)∂(r²u_r)/∂r + (1/sin θ)∂(u_θ sin θ)/∂θ = 0 is
r × the Appendix-B divergence; `spherical_continuity_residual` returns the App. B normalisation. The analysis first
listed it as a slip against a hybrid form the book never prints — retracted.

**⚠️ Book slips taught corrected (ch06, analysis §9).** (6.61) 1/z² coefficient −(Ud/π + Γ²/4π²) (printed
Ud/π − Γ²/4π² and an extra outer square; harmless) · (6.104) middle bracket +u_s/a³ (printed −; `printed_bracket=True`
reproduces it, off by u_s) · (6.108) stray dφ · "(6.8)" for (6.15) in the source velocities · "Figure 6.5" for 6.3 ·
"(6.5) and (6.12)" reversed and "(6.43)" for (6.44) · "Section 3" for §6.3 · "first-order" central differences are
second-order accurate · Example 6.2 loop index and Δψ in m²/s.

**⚠️ ρ defaults per medium (ch06 review S4).** 2-D force helpers (`lift_per_span`, `surface_pressure_force`,
`blasius_force`, `cv_force_on_body`, `cylinder_circulation_state`, `laurent_contributions`, `blasius_state`,
`force_on_held_singularity`, `ComplexFlow.pressure`) default to **air 1.2 kg/m³**; Example 6.1, the §6.9 sphere family,
`flow_field_callables` and `ideal_flow_residuals` default to **water 1000 kg/m³**. Teaching code passes ρ explicitly.

**⚠️ ω is an angular frequency from ch07 on (ch07 §7.1).** In ch02–ch06 ω was the vorticity (and ch03's `omega0` a
rotation rate); in Ch. 7 ω = 2π/T [rad/s] is the (intrinsic) angular frequency of a wave, ω ≥ 0, and the **sign of k**
carries the direction (cos(kx − ωt) with k < 0 travels to −x). `phase_speed` returns the speed |ω/k| ≥ 0;
`group_velocity` and `group_velocity_numeric` return the **signed** dω/dk (negative for a left-going wave; review M1).
The cyclic frequency ν = 1/T [Hz] (`wave_parameters(nu=)`) is 2π smaller than ω and is **not** the kinematic viscosity
(which `viscous_decay(a0, k, nu, t)` does take). A probe in a current U sees ω₀ = ω + U·K (7.9); every other ω is intrinsic.

**⚠️ η is the surface *elevation* in ch07, the surface *function* in ch04.** ch04 (4.90): the surface is {η = 0} and
n = ∇η/|∇η|. ch07 (7.14): η(x, t) is the height of the surface above z = 0, the level-set function is f = z − η and the
upward normal is (−η_x, 1)/√(1 + η_x²). Call `core.interfaces.kinematic_bc_residual` with z − η (as
`ch04.linear_wave_fields` already does); passing the elevation fails a test.

**⚠️ ζ and θ each mean two things inside ch07.** ζ = a particle's vertical excursion from its mean position (§7.2, §7.6,
§7.8; (7.35b), (7.149)) **and** the interface displacement (§7.7). θ = the local phase θ(x, t) = kx − ωt (7.72) **and**
the angle of K above the horizontal (7.139), cos θ = |k|/K — which is the angle of the beam (c_g, particle motion) from
the **vertical** (Fig. 7.33's "45° with the horizontal" coincides only because 45° is symmetric). Code names
`theta_phase`, `theta_K`; explainers label both angles. Earlier ζ: parcel displacement (ch01), relative vorticity
(ch05, Ch. 13), Zhukhovsky plane (ch06).

**⚠️ Two reduced gravities (ch07 (7.117) vs ch04).** (7.117) g′ = g(ρ₂ − ρ₁)/ρ₂ divides by the **lower (heavier)**
density — it is what (7.113) gives as kH → 0; ch04's `core.similarity.reduced_gravity` divides by ρ₁ (the upper). They
differ by the factor ρ₁/ρ₂ (ρ₁ = 1025, ρ₂ = 1027: 0.01910 vs 0.01914 m/s²); `core.waves.reduced_gravity_book(rho1,
rho2, ref="lower"|"upper")` names the choice; `ch07.reduced_gravity` is ch04's ρ₁ form re-exported. Ch. 13 often uses
ρ₀. Always say which.

**⚠️ The printed internal-wave formulas assume k > 0 (ch07 (7.138), (7.145)).** ω = kN/K and
c_g = (Nm/K³)(m e_x − k e_z) are right only for k > 0; the book's own Fig. 7.29 draws K up-left (k < 0), where the printed
c_g points the wrong way. Code: ω = N|k|/K and c_g = ∇_K ω (`internal_wave_velocities`, `printed=True` keeps the
printed form as a wrong variant); θ from cos θ = |k|/K, never arctan(m/k). Horizontal parts of c and c_g share a sign;
vertical parts are opposite (phase up ⇔ energy down).

**⚠️ Complex notation from ch07 §7.7 on.** Fields are Re{A e^{i(kx − ωt)}} with the Re dropped during linear algebra
(∂/∂x → ik, ∂/∂t → −iω, ∇² → −K²); a complex amplitude carries size and phase (b = |b|e^{iφ}). **Take real parts before
multiplying**: ⟨Re(Ae^{iθ})Re(Be^{iθ})⟩ = ½Re(AB*) (P178); ρ′ lags w by 90° (factor i in (7.153)), so ⟨gρ′w⟩ = 0.
`core.waves.real_field(amp, phase)` restores the real field. Contrast ch06: z = x + iy was a *position* and w = φ + iψ a
complex potential — there the imaginary part meant ψ, here it means a quarter-period phase shift.

**⚠️ p′ has two definitions in ch07.** §7.2 (7.30): p′ ≡ p + ρgz with p **gauge** (atmosphere = 0) — the deviation from
the still-water hydrostatic −ρgz; §7.8 (7.124): p′ = p − p̄(z) about a stratified hydrostatic base p̄ with dp̄/dz = −ρ̄g.
The Boussinesq ρ′ is likewise ρ − ρ̄(z).

**⚠️ Energy per what? (ch07).** Surface and interfacial E, E_k, E_p (7.39)–(7.42), (7.96) are per unit **horizontal
area** [J/m²] (depth-integrated) and the flux F (7.44) per unit **crest length** [W/m]; internal-wave E (7.157) is per
unit **volume** [J/m³] and F (7.158)–(7.159) per unit **area** [W/m²]; the hydraulic-jump E = u²/2 + gH is per unit
**mass** of a surface particle [J/kg] (head loss ΔE/g [m]). ⟨η²⟩ (overbar) is a wavelength average, ⟨·⟩ a time average;
⟨η²⟩ = a²/2 only for a sinusoid.

**⚠️ Book slips taught corrected (ch07, analysis §9).** (7.66) envelope ½Δω **x** → ½Δω **t** (`beat_wave(printed=True)`
is frozen) · "u from (7.28)" before (7.44) means (7.27) · (7.105) e^{i(k**z** − ωt)} → e^{i(k**x** − ωt)}
(`two_layer_residuals(printed_7_105=True)` fails) · (7.98) ∂φ₁/**d**z (cosmetic) · p. 254 "y = 0" → z = 0 · the Ursell
remark's "(7.88)" means (7.87) · (7.40) stray comma · p. 288 interfacial E_p middle form over λ/2 gives ⅛Δρga² (read
/λ; ¼ is right) · p. 288 cites Exercise 7.16 for the interfacial E_k (it is 7.18) · the printed (7.138)/(7.145) k > 0
assumption (above). Not slips: Stokes' γ = 1 in (7.83) is right (a literal truncated exercise set-up gives 3/8).

## Register

| Symbol | Meaning | SI unit | Convention / sign | Chapters | Code name |
|---|---|---|---|---|---|
| **Coordinates, time, kinematics** | | | | | |
| x, y, z | Cartesian coordinates | m | **z positive upward** (§1.7, §1.10); depth = −z; in the Couette gap y runs from the fixed wall (0) to the moving plate (h) | ch01 → | `x`, `y`, `z` |
| t ⚠️ | time | s | | ch01 → | `t` |
| u | velocity component along x (Couette profile u(y, t)) | m/s | | ch01 → | `u`, `u_hist` |
| U | speed of the moving plate; mean flow speed in a Π group | m/s | | ch01 → | `U` |
| ζ | upward displacement of a parcel from its rest height z_o | m | + up | ch01 → | `zeta`, `zeta0`, `zeta_max` |
| w₀ | initial vertical velocity of a parcel | m/s | + up | ch01 → | `w0` |
| **Molecules and the continuum** | | | | | |
| n ⚠️ | number density of molecules | 1/m³ | | ch01 → | `number_density` |
| n ⚠️ | number of molecules in a volume (Eq. 1.21) | – | | ch01 | `molecular_gas_pressure(n, V, T)` |
| n ⚠️ | number of dimensional variables (§1.11) | – | | ch01 → | (len of `variables`) |
| m ⚠️ | mass of one molecule; mass of a pendulum bob (preset) | kg | m = M_w / A_o | ch01 → | `m`, `molecule_mass` |
| N̄ | mean molecule count in a sampling box | – | relative noise ≈ N̄^(−1/2) | ch01 | `density_noise_expected` |
| L ⚠️ | side of the sampling box; body size in Kn | m | | ch01 | `L`, `L_box` |
| L_flow | length over which the macroscopic density varies | m | | ch01 | `L_flow` |
| ε ⚠️ | relative amplitude of the imposed density variation (E1) | – | | ch01 | `variation`, `gradient` |
| ε ⚠️ | pipe wall roughness height (§1.11) | m | | ch01 | `eps` in `PIPE` |
| l | mean free path | m | Jennings formula for air | ch01 → | `l`, `mean_free_path_jennings`, `mean_free_path_air` |
| Kn | Knudsen number l/L | – | continuum when Kn ≪ 1 | ch01 → | `knudsen_number` |
| k_B | Boltzmann constant | J/K | exact (SI 2019) | ch01 → | `K_B` |
| N_A, A_o | Avogadro constant per mol; per kilomole (book's A_o) | 1/mol; 1/kmol | kmol throughout the thermodynamics | ch01 → | `N_A`, `N_A_KMOL` |
| M_w | molecular weight (mass of one kilomole) | kg/kmol | air 28.9644 (USSA-1976) | ch01 → | `M_w`, `M_W_AIR`, `MOLAR_MASS` |
| **Transport (§1.5)** | | | | | |
| τ ⚠️ | shear stress τ_xy = μ ∂u/∂y | Pa | signed; the stress on the moving top plate is −μ ∂u/∂y | ch01 → | `newton_shear_stress`, `wall_shear_history` |
| τ ⚠️ | pendulum period (preset variable) | s | | ch01 | `tau` in `PENDULUM` |
| τ_y | yield stress of a Bingham material | Pa | | ch01 | `tau_y` |
| μ | dynamic viscosity | Pa s | ≥ 0 (second law) | ch01 → | `mu`, `sutherland_viscosity` |
| ν | kinematic viscosity μ/ρ (momentum diffusivity) | m²/s | never interchange with μ | ch01 → | `nu`, `kinematic_viscosity` |
| D ⚠️ | diffusivity in the model PDE ∂f/∂t = D ∂²f/∂y² (ν, κ or κ_m) | m²/s | FTCS stable when r = DΔt/Δy² ≤ ½ | ch01 → | `D` (`core.diffusion`) |
| D ⚠️ | blast-front radius (Ex. 1.4); sphere diameter (drag preset) | m | | ch01 | `D` in `BLAST`, `SPHERE_DRAG` |
| r | FTCS diffusion number DΔt/Δy² | – | ≤ ½ enforced | ch01 → | (inside `ftcs_diffusion_1d`, `stable_time_step`) |
| κ_m | mass diffusivity of a species | m²/s | | ch01 → | `kappa_m` |
| κ | thermal diffusivity k/(ρC_p) | m²/s | | ch01 → | `thermal_diffusivity` |
| k ⚠️ | thermal conductivity | W/(m K) | ≥ 0 | ch01 → | `k` (`fourier_heat_flux`, `thermal_diffusivity`) |
| k ⚠️ | exponent vector of a Π group (null-space vector) | – | | ch01 | (`pi_groups` dict values) |
| Y | mass fraction of a species | – | ΣY = 1 | ch01 → | `mass_fractions`, `grad_Y` |
| J_m | species mass flux −ρκ_m∇Y | kg/(m² s) | down the gradient | ch01 → | `fick_mass_flux` |
| **q** ⚠️ | heat-flux vector −k∇T | W/m² | down the gradient | ch01 → | `fourier_heat_flux` |
| γ ⚠️ | shear strain angle of a deforming element (§1.3) | rad | | ch01 | `shear_deformation_history` |
| G ⚠️ | elastic shear modulus of a solid | Pa | τ = Gγ | ch01 | `G` |
| **Surface tension and statics** | | | | | |
| σ | surface tension (force per length = energy per area) | N/m | | ch01 → | `sigma`, `surface_tension_water` |
| R₁, R₂ | principal radii of curvature of an interface | m | signed; the concave side has the higher pressure | ch01 → | `R1`, `R2` |
| R ⚠️ | radius of a capillary tube or drop | m | | ch01 | `R` in `capillary_rise` |
| α ⚠️ | meniscus angle measured from the horizontal (Ex. 1.1) | rad | α = 90° − θ_c; α = 90° fully wetting | ch01 | `alpha`, `alpha_deg`, `alpha_from_contact_angle` |
| θ_c | contact angle (usual definition) | rad | | ch01 | `theta_c` |
| θ ⚠️ | wedge angle of the isotropy element (Fig. 1.5 idea) | rad | | ch01 | `theta` in `wedge_*` |
| h ⚠️ | capillary rise height | m | negative for non-wetting liquids | ch01 | `capillary_rise` |
| h ⚠️ | gap width between the Couette plates | m | | ch01 → | `h` |
| p | pressure | Pa | **absolute** unless the name says gauge | ch01 → | `p`, `p1`, `p2` |
| p_gauge | pressure above atmospheric | Pa | p − P_ATM | ch01 → | `gauge_pressure`, `absolute_pressure` |
| p0 ⚠️ | pressure at z = 0 | Pa | not the θ reference | ch01 → | `p0` |
| ρ | density | kg/m³ | | ch01 → | `rho`, `rho0` |
| g | gravitational acceleration | m/s² | standard g = 9.80665, acts in −z | ch01 → | `g`, `G0` |
| V | volume of a body or of gas | m³ | | ch01 → | `V`, `volume` |
| H ⚠️ | scale height RT/g of an isothermal atmosphere | m | ≈ 8.43 km at 288.15 K | ch01 → | `scale_height`, `H` in `SCALE_HEIGHT` |
| **Thermodynamics (§1.8–1.9)** | | | | | |
| T ⚠️ | absolute temperature | K | kelvin inside every function; °C only at the boundary | ch01 → | `T`, `T0`, `T1`, `T2` |
| T ⚠️ | time dimension in a dimension vector (§1.11) | – | basis (M, L, T, Θ) | ch01 → | `BASIS` |
| e | internal energy per unit mass | J/kg | state function | ch01 → | `e`, `perfect_gas_internal_energy` |
| h ⚠️ | enthalpy per unit mass e + pv | J/kg | state function | ch01 → | `enthalpy`, `perfect_gas_enthalpy` |
| v | specific volume 1/ρ | m³/kg | | ch01 → | `v`, `v1`, `v2`, `specific_volume` |
| q ⚠️ | heat added TO the system per unit mass | J/kg | path function; + in | ch01 → | `q`, `q_path` |
| w | work done ON the system per unit mass | J/kg | path function; reversible w = −∫p dv | ch01 → | `w` (in `path_heat_work_totals`) |
| s ⚠️ | entropy per unit mass | J/(kg K) | state function | ch01 → | `s`, `perfect_gas_entropy_change` |
| C_p, C_v | specific heats at constant pressure / volume | J/(kg K) | air: 1004.70, 717.64 (from γ = 1.4) | ch01 → | `cp`, `cv`, `CP_AIR`, `CV_AIR` |
| γ ⚠️ | ratio of specific heats C_p/C_v | – | air 1.4 | ch01 → | `gamma`, `GAMMA_AIR` |
| R ⚠️ | gas constant per kilogram R_u/M_w | J/(kg K) | air 287.058 (CODATA R_u) | ch01 → | `R`, `R_AIR`, `gas_constant` |
| R_u | universal gas constant per kilomole k_B A_o | J/(kmol K) | per kmol, not per mol (factor 1000 trap) | ch01 → | `R_U`, `Ru` |
| g ⚠️ | Gibbs free energy per unit mass h − Ts (device in D19) | J/kg | not gravity | ch01 | (sympy symbol in the D19 check) |
| c ⚠️ | speed of sound √((∂p/∂ρ)_s) | m/s | | ch01 → | `c`, `sound_speed_from_eos`, `perfect_gas_sound_speed` |
| α ⚠️ | thermal expansion coefficient −(1/ρ)(∂ρ/∂T)_p | 1/K | perfect gas 1/T | ch01 → | `alpha`, `thermal_expansion_coefficient` |
| K₀, n_T | bulk modulus and exponent of the Tait equation for water | Pa; – | n = 7.15 | ch01 | `K0`, `n`, `TAIT_N_WATER` |
| a, b ⚠️ | van der Waals constants (per mass) | Pa m⁶/kg²; m³/kg | | ch01 | `a`, `b`, `VDW_CO2` |
| **Stratification (§1.10)** | | | | | |
| z_o | rest height of a parcel | m | | ch01 → | `z0` |
| ρ_a | density a parcel would have after an isentropic displacement; dρ_a/dz its gradient | kg/m³; kg/m⁴ | ocean: −ρg/c²; incompressible parcel: 0 | ch01 → | `drho_a_dz`, `isentropic_density_gradient` |
| N² | squared Brunt–Väisälä frequency | 1/s² | > 0 stable, = 0 neutral, < 0 unstable | ch01 → | `N2`, `brunt_vaisala_sq*` |
| N | Brunt–Väisälä frequency | rad/s (1/s) | period 2π/N; e-folding 1/√−N² | ch01 → | `stability_timescale` |
| Γ ⚠️ | environmental lapse rate | K/m (shown in K/km) | **Kundu Γ ≡ dT/dz** in code; meteorology −dT/dz shown alongside | ch01 → | `dT_dz`, `Gamma`, `lapse_rate`, `lapse_rate_convention` |
| Γ_a ⚠️ | adiabatic lapse rate −gαT/C_p (= −g/C_p perfect gas) | K/m | negative in code; stable when Γ > Γ_a | ch01 → | `Gamma_a`, `adiabatic_lapse_rate` |
| θ ⚠️ | potential temperature T(p_ref/p)^((γ−1)/γ) | K | stable when dθ/dz > 0 | ch01 → | `theta`, `dtheta_dz`, `potential_temperature` |
| p_ref | reference pressure of θ and ρ_θ | Pa | 1.0e5 Pa (book: surface pressure) | ch01 → | `p_ref`, `P_REF` |
| ρ_θ | potential density | kg/m³ | θρ_θ = p_ref/R | ch01 → | `potential_density` |
| S ⚠️ | salinity | g/kg (used as practical salinity) | parcel keeps S | ch01 → | `S`, `S0` |
| α_T, β_S | thermal expansion and haline contraction coefficients of the linear seawater EOS | 1/K; kg/g | defaults at 10 °C, 35 g/kg | ch01 | `alpha_T`, `beta_S` |
| **Dimensional analysis (§1.11)** | | | | | |
| q_i | a dimensional variable of the problem | varies | | ch01 → | keys of `variables` |
| [q] | dimension vector of q in (M, L, T, Θ) | – | amount of substance dropped | ch01 → | `dimension_vector` |
| M, L, Θ | mass, length, temperature dimensions | – | | ch01 → | `BASIS` |
| A ⚠️ | dimensional matrix (rows = dimensions, columns = variables) | – | | ch01 → | `dimensional_matrix` |
| r | rank of A | – | largest nonzero minor | ch01 → | `rank_by_minors` |
| Π_j | dimensionless group | – | exact `Fraction` exponents; n − r of them | ch01 → | `pi_groups`, `group_value` |
| λ_M, λ_L, λ_T ⚠️ | factors scaling the base units in a change of unit system | – | | ch01 | `UNIT_SYSTEMS`, `rescale_units` |
| λ ⚠️ | light wavelength (Ex. 1.5) | m (nm in `wavelength_to_rgb`) | | ch01 | `lam`, `lam_nm` |
| λ ⚠️ | root of the characteristic equation of ζ″ + N²ζ = 0 (D18) | 1/s | ±√−N² | ch01 | (notebook, explainer D18) |
| Δp, Δx, d | pressure drop, pipe length, pipe diameter | Pa, m, m | | ch01 → | `dp`, `dx`, `d` in `PIPE` |
| E, K ⚠️ | blast energy; Taylor's constant in E = KρD⁵/t² | J; – | K = 0.856 free sphere, hemisphere ≡ sphere of 2E | ch01 | `E`, `K`, `TAYLOR_K_GAMMA14`, `geometry` |
| S ⚠️, I | scattered and incident light intensity (Ex. 1.5) | W/m² | S/I ∝ λ⁻⁴ | ch01 | `S`, `I` in `RAYLEIGH` |
| n_s | refractive index of the scatterer | – | | ch01 | `n_s` |
| φ ⚠️ | unknown function of the remaining groups | – | | ch01 → | `phi3`, `pythagoras_phi` |
| **Tensor algebra (ch02 §2.1–2.8)** | | | | | |
| i, j, k, l, m, n (subscripts) | index letters, each running over 1, 2, 3 (Python 0, 1, 2) | – | repeated in a term = summed (dummy); once = free; three times = error | ch02 → | strings in `expand_indices`, `classify_indices` |
| λ^k, b^k, x^(k) | superscript k is a **label** (k-th eigenvalue/eigenvector), not a power | – | | ch02 → | `lam[k]`, `B[:, k]` |
| e_i, e'_j | unit vectors of the old and the rotated (primed) frame | – | right-handed, same origin | ch02 → | `unit_vectors`, `E_old`, `E_new` |
| x_i, x'_j | components of the position vector in the old / new frame | m | x' = Cᵀx | ch02 → | `x`, `xp` |
| C_ij ⚠️ | direction cosine e_i·e'_j (row old, column new) | – | **passive**; orthogonal, det +1; equals active R(θ) applied as Rᵀ | ch02 → | `C`, `direction_cosines`, `rotation_matrix_2d/3d` |
| C ⚠️ | closed boundary curve of an open surface A (§2.13) | – | orientation from n (see Stokes trap) | ch02 → | `Loop` |
| θ ⚠️ | rotation angle of the axes (Ex. 2.1: polar angle) | rad | + counterclockwise about e₃; `_deg` only at the interface | ch02 | `theta`, `theta_deg`, `rotation_angle` |
| u_r, u_θ | polar components of a plane vector | same as u | (2.5) with j ∈ {r, θ} | ch02 | `polar_components` |
| δ_ij | Kronecker delta | – | δ_ii = 3 | ch02 → | `kronecker_delta` |
| ε_ijk | alternating (Levi-Civita) tensor | – | +1 cyclic, −1 anticyclic, 0 repeated; pseudotensor under reflections | ch02 → | `levi_civita`, `permutation_sign` |
| A:B ⚠️ | double dot | product of the units | **book A_ij B_ji**; Frobenius A_ij B_ij on request | ch02 → | `double_dot(A, B, convention=)` |
| I₁, I₂, I₃ | invariants: trace, ½(I₁² − A_ij A_ji), det | powers of A's unit | coefficients of λ³ − I₁λ² + I₂λ − I₃ | ch02 → | `invariants`, `characteristic_polynomial` |
| λ ⚠️ | eigenvalue of a symmetric tensor (principal value) | unit of the tensor | `principal_axes` returns them ascending; `example_2_4` keeps (Γ, −Γ) | ch02 → | `lam` |
| b ⚠️ | eigenvector (principal direction), unit length | – | columns of B, det B = +1 | ch02 → | `B` |
| b ⚠️ | constant vector of the solid-body rotation u = b × x (Ex. 2.3) | 1/s | ∇×u = 2b | ch02 | `b` in `solid_body_rotation_field` |
| a ⚠️ | shear-stress magnitude of Ex. 2.2 (Pa); constant of u = a x in Ex. 2.3 (1/s) | Pa; 1/s | | ch02 | `a` |
| **Stress and traction (ch02 §2.4, §2.6)** | | | | | |
| τ_ij ⚠️ | stress tensor: force per area in direction j on the face with normal e_i | Pa | tensile positive; +e_i face → positive along +e_j; pressure τ = −pδ | ch02 → | `tau`, `STRESS_CUBE_FACES` |
| n ⚠️ | unit normal of a plane or surface | – | outward on closed surfaces | ch02 → | `n` |
| f | traction (force per area) on the plane with normal n | Pa | f_i = τ_ji n_j (first index) | ch02 → | `traction`, `traction_2d` |
| σ_n, τ_s | normal and shear parts of the traction (σ_n = f·n, τ_s = \|f − σ_n n\|) | Pa | σ_n tensile positive; τ_s ≥ 0 (direction returned separately) | ch02 → | `sigma_n`, `tau_s`, `normal_shear_stress` |
| φ ⚠️ | angle of a cutting plane's normal in the x₁–x₂ plane (E2, Ex. 2.2) | rad | measured from e₁ | ch02 | `phi`, `stress_vs_angle` |
| dA, dA_i | surface element and its vector n dA; tetrahedron face areas n_i dA | m² | | ch02 → | `tetrahedron_face_areas` |
| h ⚠️ | size of a shrinking element (tetrahedron, box, square loop) | m | face terms ∝ h², volume terms ∝ h³ | ch02 → | `h` in `integral_*` |
| **Fields and operators (ch02 §2.9–2.14)** | | | | | |
| h, (hx, hy, hz) ⚠️ | grid spacing | m | ordered (x, y, z); nodes, h = (b − a)/(n − 1); periodic h = L/n | ch02 → | `Grid.h`, `spacing` |
| d | coordinate direction 0, 1, 2 = x₁, x₂, x₃ | – | lives on array axis ndim − 1 − d | ch02 → | `direction`, `axis_of_direction` |
| φ ⚠️ | scalar field whose gradient is taken (potential) | varies | | ch02 → | `phi`, `ScalarField`, `potential_field` |
| ∇φ, ∂φ/∂n | gradient; directional derivative ∇φ·n | [φ]/m | ∇φ ⊥ level sets, points uphill | ch02 → | `gradient`, `directional_derivative` |
| u | velocity vector field (a generic vector field in ch02) | m/s | components on array axis 0 | ch02 → | `u`, `VectorField` |
| Q | generic scalar/vector/tensor field in Gauss' theorem (2.30) | varies | | ch02 | `Q_fn` |
| G_ij ⚠️ | velocity gradient ∂u_i/∂x_j | 1/s | **row = component, column = derivative direction** | ch02 → | `G`, `vector_gradient` |
| S_ij | symmetric part ½(G + Gᵀ) (strain-rate tensor) | 1/s | Ex. 2.4 Γ ≡ S₁₂ (book's "2S₁₂ = Γ" is inconsistent) | ch02 → | `symmetric_part`, `strain_rate_tensor` |
| A_ij ⚠️ | antisymmetric part ½(G − Gᵀ) = ½R | 1/s | its vector is ½∇×u; not the dimensional matrix of ch01 | ch02 → | `antisymmetric_part` |
| R_ij ⚠️ | book's rotation tensor G − Gᵀ = −ε_ijk ω_k | 1/s | vector = ∇×u; not a gas constant or a radius here | ch02 → | `rotation_tensor`, `antisymmetric_from_vector` |
| ω ⚠️ | vector of R = vorticity ∇×u | 1/s | ω_k = −½ ε_ijk R_ij; the element spins at ½ω (ch03); **ch03 Ex. 3.1 and ch07 also use ω for an angular frequency** | ch02 → | `omega`, `vector_from_antisymmetric(rotation_tensor(G))` |
| ½ω₃ = A₂₁ | spin rate of a fluid element in a plane flow (E3 readout) | 1/s | never written "ω" | ch02 | `vector_from_antisymmetric(antisymmetric_part(G))` |
| Γ ⚠️ | rate of a linear-flow preset: `simple_shear` u₁ = Γx₂ ⇒ Γ = du₁/dx₂ = 2S₁₂ (= ch03 γ); Ex. 2.4's Γ is S₁₂ itself (half as big) | 1/s | **ch01 Γ = lapse rate; ch03 Γ = circulation** — see the Γ trap table above | ch02, ch03 presets | `Gamma` in `velocity_gradient_preset`, `example_2_4` |
| Γ_circ ⚠️ | circulation ∮u·t ds | m²/s | + counterclockwise about n; ch02 wrote Γ_circ to avoid the clash; **ch03 already calls it Γ** (not only ch05) | ch02 → (ch03 calls it Γ) | `circulation` |
| K ⚠️ | strength of the irrotational vortex u_θ = K/r | m²/s | Γ_circ = 2πK round the core; not Taylor's K | ch02 | `K` in `irrotational_vortex_field` |
| m ⚠️ | 2-D point-source strength (outflux per unit depth) | m²/s | not a molecule mass | ch02 | `m` in `point_source_field` |
| V, A ⚠️ | volume of a region; area of an open surface (or of a loop) | m³; m² | A here is an area, not a tensor or matrix | ch02 → | `volume`, `Loop.area` |
| n_c | in-surface normal to the rim C, pointing **into A** | – | t = n_c × n | ch02 → | `boundary_tangent(n_c, n)` |
| t ⚠️ | unit tangent of a curve (not time here) | – | counterclockwise about n | ch02 → | `Loop.tangents` |
| s ⚠️, ds | arc length along a curve (not entropy here) | m | | ch02 → | `Loop.ds` |
| u_i,j | comma notation for ∂u_i/∂x_j (2.36) | 1/s | a comma index transforms as a vector index | ch02 → | `comma_to_partial` |
| singular_at | point where a test field is not differentiable | m | Stokes/Gauss report `hypothesis_ok = False` | ch02 → | `VectorField.singular_at` |
| **Particles, fields and flow lines (ch03 §3.1–3.3)** | | | | | |
| u(x, t) | velocity field callable | m/s | x of shape (d,) or (d, N), coordinates on axis 0; returns the shape of x | ch03 → | `u` (kinematics convention), `as_coord_field`, `from_coord_field` |
| F(x, t) | any scalar (or vector) field followed by D/Dt | [F] | | ch03 → | `F` |
| r(t; r_o, t_o) | trajectory of the particle that was at r_o at time t_o (a label, not a variable) | m | labels constant along the path | ch03 → | `r_of_t`, `pathline` |
| X ⚠️ | Lagrangian label of D01 (x = Xe^{αt}) | m | not a grid array | ch03 | `lagrangian_map_example(X, t, alpha)` |
| α ⚠️ | stretching rate of D01's map | 1/s | not ch01's thermal expansion; not Fig. 3.11's angle | ch03 | `alpha` |
| D/Dt | material derivative ∂/∂t + u·∇ (3.5) | [F]/s | local = ∂F/∂t, advective = u·∇F | ch03 → | `material_derivative(_terms)` (local, advective, total) |
| s ⚠️, e_u | arc length along a streamline; unit vector along u | m; – | (3.6): DF/Dt = ∂F/∂t + \|u\| ∂F/∂s (book drops F) | ch03 | `streamwise_derivative`, `streamline(s_max=…)` |
| t′ ⚠️, t_o, ξ_o | Ex. 3.1: drawing instant; release time (streak-line label); orbit amplitude | s; s; m | t′ ≠ t_o; primes here are NOT a moving frame | ch03 | `example_3_1(t_prime, xi0, omega)`, `streakline(t_release=…)` |
| ω ⚠️ | Ex. 3.1's angular frequency of the uniform oscillating flow (its vorticity is 0) | rad/s | **not vorticity** in `example_3_1`, `unsteady_flow_preset(omega=…)` | ch03 (ch07 again) | `omega` |
| U ⚠️ | constant velocity of the moving frame O′ (Galilean); free stream past the cylinder | m/s | x = x′ + Ut + x′_o, u′ = u − U | ch03 → | `U` in `galilean_transform`, `cylinder_flow` |
| x′, u′, t′ ⚠️ | position, velocity and time in the translating frame O′ | m, m/s, s | ch02's primes were a rotated frame | ch03 | `xp`, `tp`, `x0p` |
| U_frame | observer speed on E3's slider (0 = fluid frame, U = body frame) | m/s | observer moves at −U_frame e_x relative to the far fluid | ch03 | `U_frame` in `cylinder_flow`, `frame_acceleration_terms` |
| a ⚠️ | cylinder radius (§3.1, §3.3); lower limit a(t) (Leibniz); ellipse semi-axis a(t) | m | | ch03 | `a` |
| ψ | streamfunction, u = ∂ψ/∂y, v = −∂ψ/∂x | m²/s | constant along streamlines | ch03 → (Ch. 4, 6) | `cylinder_streamfunction` |
| T ⚠️ | temperature field of the thermal front (E2) | K | southerly wind + north–south gradient = warm advection | ch03 | `thermal_front(..., grad_K_per_m, heating_K_per_s, front_speed)` |
| **Coordinates (ch03 §3.1)** | | | | | |
| R ⚠️, φ, z | cylindrical radius (capital R), azimuth, axial coordinate | m, rad, m | φ = atan2(y, x); φ = 0 on the axis | ch03 → | `cylindrical_from_cartesian` |
| r ⚠️, θ ⚠️, φ ⚠️ | spherical radius, **polar angle from +z**, azimuth | m, rad, rad | θ ∈ [0, π]; physics convention | ch03 → | `spherical_from_cartesian` |
| r, θ ⚠️ | plane polar radius and angle from +x (§3.5) | m, rad | counterclockwise positive | ch03 → | `polar_from_cartesian`, `core.vortices` |
| e_r, e_θ, e_φ, E[k] | local unit vectors (move with the point) | – | returned as rows E[k] (E = Cᵀ of ch02) | ch03 → | `unit_vectors_*` |
| u_r, u_θ; u_R, u_φ, u_z; u_r, u_θ, u_φ | velocity components in the local bases | m/s | projections E[k]·u | ch03 → | `velocity_components` |
| **Strain and rotation (ch03 §3.4)** | | | | | |
| dx, du | separation of a neighbour and its relative velocity du = G·dx (3.10) | m, m/s | | ch03 → | `relative_velocity(G, dx)` |
| δx, δV | material line element and material volume (carried by the flow) | m, m³ | D(δx)/Dt = δu | ch03 → | `measured_strain_rates` |
| n, n₁, n₂ ⚠️ | unit direction(s) of material lines | – | n₁ ⟂ n₂ required for the shear rate | ch03 | `linear_strain_rate(G, n)`, `shear_strain_rate(G, n1, n2)` |
| α ⚠️, β | Fig. 3.11 angles: α clockwise from the vertical side, β counterclockwise from the horizontal side | rad | S₁₂ = ½D(α + β)/Dt; spin ½D(−α + β)/Dt | ch03 | (D09, D12) |
| γ ⚠️ | shear rate du₁/dx₂ of the parallel shear flow | 1/s | S₁₂ = γ/2, ω₃ = −γ (clockwise) | ch03 | `parallel_shear_kinematics(gamma)` |
| θ ⚠️, θ̇ | angle of a material line from +x; its turning rate | rad; rad/s | shear: θ̇ = −γ sin²θ; counterclockwise positive | ch03 | `material_line_rotation_rate(G, theta)` |
| ½ω₃ | element rotation (spin) rate = average of two perpendicular lines | rad/s | = R₂₁/2; a paddle wheel's rate | ch03 → | `element_rotation_rate`, `perpendicular_pair_rotation_rate` |
| Ω ⚠️ | angular velocity of a rigid motion or of a rotating observer | rad/s | rigid motion U + Ω × x has ω = 2Ω; rotating frame: ω′ = ω − 2Ω | ch03 → (Ch. 4, 13 Earth) | `rigid_body_velocity(U, Omega, x)`, `vorticity_in_rotating_frame(omega, Omega)` |
| φ ⚠️ | velocity potential, u = ∇φ (3.17) (simply connected regions only) | m²/s | line vortex: φ = Bθ multivalued | ch03 → (Ch. 6) | `potential_velocity`, `velocity_potential_2d` |
| λ ⚠️, S̄_αα | principal strain rates (eigenvalues of S) | 1/s | `principal_strain_rates` ascending; Greek index α = no sum | ch03 → | `principal_strain_rates`, `strain_velocity_principal` |
| h ⚠️, ht | finite-difference step in space; its own step in time | m; s | never reuse h as a time step (review Should-fix 5) | ch03 → | `h`, `ht` in `material_derivative_terms`, `velocity_gradient_at` |
| **Vortices (ch03 §3.5)** | | | | | |
| ω₀ | angular velocity of solid-body rotation u_θ = ω₀r (3.22) | rad/s | ω_z = 2ω₀ | ch03 → | `omega0` in `solid_body_rotation` |
| B | line-vortex strength, u_θ = B/r (3.25) | m²/s | Γ = 2πB (3.26); book's Fig. 3.16 writes C; ch02 wrote K | ch03 → | `B` in `line_vortex` |
| ω_z | vorticity normal to the plane (polar form (3.23)) | 1/s | counterclockwise positive | ch03 → | `polar_vorticity_z`, second output of `rankine_vortex`/`gaussian_vortex` |
| Γ ⚠️ | circulation (total circulation of a vortex; of a loop) | m²/s | see the Γ trap table | ch03 → | `Gamma`, `circulation_circle` |
| σ ⚠️ | vortex core radius (Rankine, Gaussian) | m | **not surface tension (ch01)**; Lamb–Oseen σ² = 4νt | ch03 → | `sigma` |
| x*, r_max | x* = r_max²/σ² = 1.2564312 (root of 1 + 2x = eˣ); r_max = 1.1209064σ | –; m | r_max = σ√x*, not σx* | ch03 → | `gaussian_vortex_max_radius(sigma, method)` |
| **Transport (ch03 §3.6)** | | | | | |
| V*(t), A*(t) | control volume and its closed surface (may move and deform) | m³, m² | 2-D: area and arc length | ch03 → (Ch. 4) | `ControlVolume` shapes |
| b ⚠️ | velocity of the control surface | m/s | only b·n matters; b = u for a material volume; b = 0 for a fixed CV; **not** ch02's eigenvector or b × x | ch03 → | `surface_nodes(...)` fourth output |
| n ⚠️ | outward unit normal of A* | – | b·n signed | ch03 → | `surface_nodes` |
| a(t), b(t), ȧ, ḃ ⚠️ | Leibniz limits and their speeds (3.30) | m, m/s | lower term subtracted | ch03 | `leibniz_terms(F, dFdt, a, b, dadt, dbdt, t)` |
| Δt, ΔV, T1–T4 | time step of the definition (3.31); signed swept volume; the four terms of (3.32) (T4 = ∫_ΔV Δt ∂F/∂t = O(Δt²)) | s, m³ | | ch03 | `swept_terms`, `swept_terms_sphere` |
| h ⚠️, r_o, ṙ, θ | Ex. 3.2 cone height, base radius, its growth rate, half-angle | m, m, m/s, rad | b·n = 0 on the base | ch03 | `example_3_2(h, r0, rdot)`, `GrowingCone` |
| **Conservation laws: budgets (ch04 §4.1–4.4)** | | | | | |
| V(t), A(t) | material volume and its surface (move with the fluid, b = u) | m³, m² | sealed-balloon picture | ch04 → | `material=True` in the budgets; `material_interval`, `material_mass` |
| V*(t), A*(t), b | control volume, its surface, surface velocity (any motion) | m³, m², m/s | outward n; flux uses the **relative** velocity (u − b)·n, u and b in the same frame | ch04 → | `ControlVolume` shapes (ch03) passed to `mass_budget`, `momentum_budget`, `energy_budget` |
| M ⚠️ | mass in the CV ∫ρ dV (`MassBudget.dM_dt` its rate) | kg | | ch04 | `MassBudget` |
| M ⚠️ | rocket mass M(t) (Ex. 4.4) | kg | dM/dt = −ρ_eV_eA_e | ch04 | `rocket_trajectory(M0, mdot, Ve, …)` |
| M ⚠️ | torque ∫r × f dA + … (4.64) | N m | about the CV origin | ch04 | `sprinkler_torque`, `angular_momentum_budget` |
| M, Ma ⚠️ | Mach number U/c (4.111) | – | incompressible regime M < 0.3 (§4.2) | ch04 → (Ch. 15) | `mach_number`, `is_incompressible_regime` |
| P | momentum in the CV ∫ρu dV | kg m/s | vector | ch04 | `MomentumBudget.dP_dt` |
| H ⚠️ | angular momentum ∫r × ρu dV (4.64) | kg m²/s | | ch04 | `angular_momentum_budget` |
| H ⚠️ | height of the wake CV (Ex. 4.1) | m | drag independent of H once H covers the wake | ch04 | `wake_drag_per_span(…, H=)` |
| H_c | c²/g, the depth over which compressibility of a hydrostatic column matters | m | air 11.8 km; Boussinesq needs L ≪ c²/g | ch04 | `boussinesq_validity` |
| ṁ, Q | mass flow rate; volume flow (jet sheets Q₁,₂ = Q(1 ± cos θ)/2) | kg/s; m³/s | | ch04 | `cv_scenario("jet")`, `orifice_mass_flow` |
| U(y), U_∞ | wake profile, free stream (Ex. 4.1) | m/s | F_D/l = ρ∫U(U_∞ − U)dy (deficit weighted by U, not U_∞) | ch04 | `gaussian_wake`, `wake_drag_per_span` |
| F_D, F_L | drag (force on the body, + downstream), lift | N (per span N/m) | force on the fluid is −F_D (Newton III) | ch04 → | `wake_drag_per_span`, `drag_coefficient`, `lift_coefficient` |
| h_in, h_out ⚠️ | depths ahead of and behind a bore (Ex. 4.3) | m | U = √(g h_out(h_in + h_out)/(2h_in)) → √(gh) | ch04 | `bore_speed(h_in, h_out, g)` |
| V_e, A_e, ρ_e, F_S | rocket exhaust speed (relative), exit area, exit density, support force | m/s, m², kg/m³, N | Tsiolkovsky Δb = V_e ln(M₀/M₁) | ch04 | `rocket_delta_v`, `rocket_closed_form` |
| b(t) ⚠️ | rocket speed (Ex. 4.4) — the CV velocity of an accelerating CV | m/s | not ch02's eigenvector | ch04 | inside `rocket_*` |
| θ ⚠️ | angle between jet and plate (E1) | rad | θ = π/2: normal plate, F = ρV²A | ch04 | `jet_plate_force(rho, V, A, theta)` |
| a, α ⚠️, A | sprinkler arm length, nozzle angle, nozzle area (Ex. 4.6) | m, rad, m² | M = 2aρAU² cos α | ch04 | `sprinkler_torque(a, rho, A, U, alpha)` |
| Φ ⚠️ | force potential per unit mass, g = −∇Φ (4.18); gravity Φ = gz (z up) | J/kg (m²/s²) | **not** the Ch. 1 unknown function φ, not (4.99)'s scaling function | ch04 → | `gravity_potential`, `body_force_from_potential` |
| f ⚠️ | surface force (traction) per area, f_j = n_iτ_ij | Pa | on the fluid inside, outward n | ch04 → | `traction` (ch02), budgets' `traction=` |
| **Streamfunctions (ch04 §4.3)** | | | | | |
| ψ ⚠️ | 2-D stream function: ρu = ∂ψ/∂y, ρv = −∂ψ/∂x; flux between streamlines ψ₂ − ψ₁ | m²/s (ρ = 1) or kg/(m s) | see the sign trap above (GFD often uses the opposite sign) | ch03 (gloss) → ch04 → Ch. 6, 13 | `velocity_from_streamfunction_2d`, `flux_between_streamlines`, `streamfunction_preset` |
| ψ ⚠️ | axisymmetric (Stokes) stream function: ρu_R = −(1/R)∂ψ/∂z, ρu_z = (1/R)∂ψ/∂R | m³/s (ρ = 1) | ring flux 2πΔψ | ch04 → Ch. 6, 8 | `velocity_from_streamfunction_axisym` |
| χ, ψ | two stream functions (stream surfaces) of 3-D steady flow, ρu = ∇χ × ∇ψ (4.12) | – | 2-D: χ = −z; axisymmetric: χ = −φ | ch04 | `mass_flux_from_stream_functions`, `stream_surface_check`, `stream_tube_mass_flux` |
| Ψ ⚠️ | vector potential of the mass flux, ρu = ∇ × Ψ (4.12) | kg/(m s) | **(4.99)/(4.106) also use Ψ for a scaling function** | ch04 | `mass_flux_from_vector_potential` |
| **Stress and the Newtonian law (ch04 §4.5–4.6)** | | | | | |
| τ_ij | total stress = −pδ_ij + σ_ij (4.27) | Pa | tensile positive; symmetric (4.25); traction contracts the first index | ch02 → | `newtonian_stress`, `total_stress`, `static_stress` |
| σ_ij ⚠️ | viscous ("deviatoric") stress | Pa | traceless only if μ_v = 0 or ∇·u = 0; **σ also = surface tension (ch01, §4.10) and vortex core radius (ch03)** | ch04 → | `viscous_stress`, `newtonian_viscous_stress_field` |
| K_ijmn | fourth-order coefficient tensor of the linear law σ_ij = K_ijmnS_mn (4.28) | Pa s | isotropic: λδ_ijδ_mn + μδ_imδ_jn + γδ_inδ_jm (4.29) | ch04 | `isotropic_fourth_order(lam, mu, gam)`, `linear_stress(K, S)` |
| λ ⚠️ | second (Lamé-type) viscosity coefficient in (4.31) | Pa s | λ = μ_v − ⅔μ; **not** bulk viscosity; ch02/ch03 λ = eigenvalue | ch04 → | `lam` in `newtonian_stress`, `lam_from_bulk` |
| γ ⚠️ | third coefficient of (4.29); only μ + γ acts on a symmetric S, "γ = μ" (4.30) names it | Pa s | **γ also: ratio of specific heats (ch01), shear rate (ch03), Ex. 4.7's ζ/δ** | ch04 | `gam` in `isotropic_fourth_order` |
| μ_v | bulk viscosity λ + ⅔μ | Pa s | Stokes assumption μ_v = 0 (4.36); ≥ 0 by the second law | ch04 → Ch. 15 | `mu_v`, `bulk_viscosity`, `stokes_assumption_holds` |
| p̄ | mean (mechanical) pressure −⅓τ_ii (4.33) | Pa | p − p̄ = μ_v∇·u (4.34); plane stress needs τ₃₃ | ch04 → | `mean_pressure(tau, tau33=)`, `thermodynamic_pressure_from_stress`, `pressure_difference` |
| S_mm | trace of S = ∇·u (volumetric strain rate, (3.14) in §3.4) | 1/s | book cites "Section 3.6" | ch04 | inside `newtonian_stress`, `deviatoric_part` |
| dev S | S − ⅓S_mmδ, the shape-changing part | 1/s | | ch04 → | `deviatoric_part` |
| G ⚠️ | velocity gradient ∂u_i/∂x_j (input of the stress functions) | 1/s | **`plane_poiseuille(G=)` is −dp/dx [Pa/m]; `ch04.G` is g** | ch02 → | `newtonian_stress(G, p, mu, …)`, `stress_on_plane(G, …)` |
| α ⚠️ (cube) | spin-up rate of a cube with τ₁₂ ≠ τ₂₁: 6(τ₁₂ − τ₂₁)/(ρh²) | rad/s² | diverges as h → 0 unless τ symmetric | ch04 | `cube_spin_acceleration(tau12, tau21, rho, h)` |
| h ⚠️ | cube side (D08); stencil step; channel half-gap/gap in exact solutions | m | | ch04 | `h` |
| **Navier–Stokes and exact solutions (ch04 §4.6)** | | | | | |
| ω, ∇×ω | vorticity and its curl; viscous force −μ∇×ω when ∇·u = 0 (4.40) | 1/s, 1/(m s) | solid body: ω ≠ 0 but ∇×ω = 0 (no viscous force) | ch03 → | `viscous_force_forms` (laplacian / div2S / curl) |
| exact solutions | Couette(–Poiseuille), plane and pipe Poiseuille, Stokes' first problem, Taylor–Green, Lamb–Oseen, ideal cylinder, solid body with gravity | – | residual of (4.39b) < 1e-6 of the largest term | ch04 → Ch. 8 | `exact_solution(name, x, t, **p)`, `EXACT_SOLUTIONS`, `ns_terms_preset` |
| η ⚠️ (Stokes) | similarity variable y/(2√(νt)) of u = U erfc η | – | **η also = surface function (4.90)** | ch04 | `stokes_first_problem(y, t, U, nu)` |
| **Noninertial frames (ch04 §4.7)** | | | | | |
| Ω, Ω̇ ⚠️ | angular velocity of the frame and its rate | rad/s, rad/s² | Earth 7.292115e-5 rad/s (WGS-84); NH Ω_z > 0; **§4.11 Ω = imposed frequency** | ch03 → | `Omega`, `dOmega_dt` in `core.rotating`; `OMEGA_EARTH` |
| U(t), dU/dt ⚠️ | origin velocity and acceleration of the noninertial frame | m/s, m/s² | a vector; a scalar non-zero dU/dt raises | ch04 | `dU_dt` in `frame_acceleration_terms`, `apparent_body_forces` |
| x′, u′, a′ ⚠️ | position, velocity, acceleration in the rotating frame (components on the turning basis e′_i) | m, m/s, m/s² | see the prime trap | ch04 | `x_prime`, `u_prime`, `a_prime`; `rotating_basis`, `inertial_velocity` |
| e′_i | turning basis vectors, de′_i/dt = Ω × e′_i | – | | ch04 | `rotating_basis`, `basis_rate(_exact)` |
| g_n, Φ_n | gravitation alone (without the centrifugal part) and its potential | m/s², J/kg | g = g_n − Ω × (Ω × x) (effective gravity) | ch04 → Ch. 13 | `effective_gravity(lat, g_n=9.8)` |
| R ⚠️ | distance from the rotation axis (cylindrical radius) | m | centrifugal Ω²R e_R, potential −½Ω²R²; **R also: principal radii R₁, R₂ (§4.10), gas constant** | ch03 → | `centrifugal_acceleration`, `centrifugal_potential` |
| f ⚠️ | Coriolis parameter 2Ω sin φ (named only) | 1/s | NH > 0; **f also traction, Helmholtz free energy** | ch04 → Ch. 13 | `coriolis_parameter(lat_rad)` |
| φ ⚠️ | latitude (in f, g_e) | rad | `_deg` only at interfaces; **φ also velocity potential (4.73), azimuth** | ch04 → | `lat_rad` |
| ζ ⚠️ | relative vorticity (glossed in D14 "what it means"; ζ + f absolute) | 1/s | **ζ also: parcel displacement (ch01), cap height and meniscus height (§4.10)** | ch04 (gloss) → Ch. 13 | — |
| Ro | Rossby number U/(2Ωl) (forward pointer) | – | see the Ro trap above | ch04 → Ch. 13 | `rossby_number(U, Omega, l, factor=2.0)` |
| **Energy (ch04 §4.8)** | | | | | |
| E | total energy per mass e + ½u_j² | J/kg | (4.53) | ch04 | `energy_budget`, `total_energy_residual_sym` |
| **q** | heat-flux vector −k∇T (4.60) | W/m² | outward q·n is a loss | ch01 → | `energy_budget(q=)` |
| ε ⚠️ | viscous dissipation rate per mass (1/ρ)σ_ijS_ij = 2ν(dev S)² + (μ_v/ρ)S_mm² ≥ 0 (4.58) | W/kg | ρε per volume [W/m³]; **ε also = alternating tensor ε_ijk (4.40), ch01 relative amplitude and roughness** | ch04 → Ch. 12, 13 | `dissipation_rate(G, rho, mu, mu_v, form)` |
| s, Ds/Dt | entropy per mass and its rate; production k\|∇T\|²/(ρT²) + ε/T ≥ 0 (4.63) | J/(kg K), W/(kg K) | requires μ, μ_v, k ≥ 0 | ch01 → | `entropy_terms`, `entropy_production` |
| ΔT_max | peak viscous heating of Couette flow μU²/(8k) | K | 1 m/s, 1 mm, water: 2.08e-4 K | ch04 | `couette_heating`, `couette_heating_transient` |
| κ ⚠️ | thermal diffusivity k/(ρC_p) (4.89) | m²/s | the book's "κ" after (4.63) means μ_v | ch01 → | `kappa` |
| **Bernoulli (ch04 §4.4, §4.9)** | | | | | |
| B ⚠️ | Bernoulli function ½\|u\|² + ∫dp/ρ + Φ (4.69) | J/kg (m²/s²) | constant on streamlines and vortex lines (4.71); everywhere if ω = 0 (4.72); **not the ch03 line-vortex strength B** | ch04 → | `bernoulli_function`, `bernoulli_head`, `bernoulli_along_line`, `rankine_bernoulli` |
| B(t) | the time-only function in unsteady potential flow (4.74) | J/kg | absorbed by φ_new = φ − ∫B dt′ | ch04 | `unsteady_bernoulli_B`, `gauge_absorbed_bracket` |
| ∫dp/ρ | pressure function of a barotropic fluid (4.67) | J/kg | kinds: constant, isothermal, isentropic (= C_p(T − T_o)) | ch04 | `pressure_function(p, p_o, kind)` |
| u × ω | Lamb vector; (u·∇)u = −u × ω + ∇(½u²) (4.68) | m/s² | solid body: +2Ω²(x, y, 0); ω × u is the sign slip | ch04 → Ch. 5 | `lamb_vector`, `lamb_identity_terms/sym` |
| φ ⚠️ | velocity potential u = ∇φ (4.73) | m²/s | irrotational, simply connected | ch03 → | `unsteady_bernoulli_pressure`, `accelerating_sphere_fields` |
| h ⚠️ | enthalpy per mass in the energy Bernoulli h + ½\|u\|² + gz (4.78); stagnation T₀ = T + U²/2C_p | J/kg | **h also: depth, heads, meniscus height** | ch01 → | `stagnation_enthalpy`, `stagnation_temperature` |
| p₀, ½ρU² | stagnation and dynamic pressure; pitot speed √(2(p₀ − p)/ρ) | Pa | | ch04 → Ch. 14, 15 | `stagnation_pressure`, `dynamic_pressure`, `pitot_speed(_from_heads)` |
| C_c | contraction coefficient of a sharp orifice | – | 0.611 (V5) | ch04 | `orifice_mass_flow(…, Cc)`, `tank_drain` |
| L, h₀ (U-tube) | column length, initial offset | m | L dU/dt + 2gh = 0 | ch04 | `u_tube_column(t, L, h0, g)` |
| **Boussinesq (ch04 §4.9)** | | | | | |
| p_s(z), ρ_s(z) | hydrostatic base state | Pa, kg/m³ | dp_s/dz = −ρ_sg | ch04 → Ch. 7, 13 | `perturbation_fields(p, rho, z, rho_s)` |
| p′, ρ′ ⚠️ | perturbations p − p_s, ρ − ρ_s | Pa, kg/m³ | primes = perturbations here | ch04 → | `p_pert`, `rho_pert` |
| ρ₀ | constant reference density | kg/m³ | replaces ρ everywhere except next to g | ch04 → | `rho0` |
| b ⚠️ | buoyancy −gρ′/ρ₀ | m/s² | upward positive; **b also CV velocity, rocket speed, eigenvector** | ch04 → Ch. 7, 13 | `buoyancy(rho_pert, rho0, g)` |
| g′ | reduced gravity gΔρ/ρ₀ | m/s² | | ch04 → | `reduced_gravity` |
| α ⚠️, δT | thermal expansion coefficient, temperature contrast | 1/K, K | validity αδT ≪ 1 (threshold 0.1) | ch01 → | `boussinesq_validity(alpha, dT, L, U, c, …)`, `boussinesq_density` |
| **Boundary conditions and surface tension (ch04 §4.10)** | | | | | |
| η(x, t) ⚠️ | surface function, surface = {η = 0} (4.90) | m (or –) | n = ∇η/\|∇η\| toward increasing η; **not the Stokes similarity variable** | ch04 → Ch. 7 | `kinematic_bc_residual(eta, u, x, t)`, `surface_preset` |
| u_s, (u − u_s)·n | surface velocity; relative normal velocity | m/s | only the normal part of u_s is defined; = 0 ⇔ no flux | ch04 → | `surface_normal_speed`, `relative_normal_velocity`, `interface_mass_flux` |
| l (pillbox) | pillbox thickness → 0 | m | side and volume terms vanish ∝ l | ch04 | `pillbox_limit` |
| R₁, R₂ | principal radii of curvature | m | higher pressure on the concave side; Δp = σ(1/R₁ + 1/R₂) (1.5) re-derived | ch01 → | `laplace_jump_from_balance`, `cap_pressure_force`, `cap_surface_tension_force` |
| ζ ⚠️ | height of the small cap (§4.10) / meniscus height above the free level (Ex. 4.7) | m | cap z = x²/2R₁ + y²/2R₂ (book prints −) | ch04 | `cap_*(…, zeta)`, `meniscus_profile_x(zeta, theta)` |
| ℓ_c, δ ⚠️ | capillary length √(σ/(Δρg)) (water 2.7 mm) = Ex. 4.7's δ; a fully wetting wall (θ = 0) lifts the meniscus √2 δ | m | Bo = (l/ℓ_c)²; Ex. 4.7 γ = ζ/δ | ch04 → | `capillary_length(sigma, rho, g, rho_other)` |
| θ ⚠️ | contact angle at the wall (Ex. 4.7) | rad | h² = 2σ(1 − sin θ)/(ρg) | ch04 | `meniscus_height(theta, …)` |
| f, F ⚠️ | Helmholtz free energy per mass e − Ts; of a system (4.94)–(4.96) | J/kg, J | **f also traction**; (∂f/∂v)_T = −p | ch04 → Ch. 15 | `core.thermo.helmholtz_free_energy` |
| **Dimensionless groups (ch04 §4.11)** | | | | | |
| l, U, Ω (scales) ⚠️ | reference length, speed, imposed frequency | m, m/s, 1/s | t* = Ωt (4.100) or Ut/l (4.109); **Ω here is not the frame rotation** | ch04 → | `Scales(l, U, rho, mu, g, Omega, …)` |
| u*, p*, t*, ∇* | scaled variables | – | p* = (p − p_∞)/ρU² (dynamic scaling) | ch04 → | `Scales.nondimensionalise/redimensionalise` |
| St, Re, Fr | Ωl/U, ρUl/μ, U/√(gl) | – | Fr is a square root of a force ratio | ch04 → | `strouhal_number`, `reynolds_number`, `froude_number` |
| Fr′, Ri, Ri_g | internal Froude U/√(g′l) (= U/(Nl)), Richardson g′l/U² = 1/Fr′², gradient Ri N²/(dU/dz)² | – | N² computed in Kundu Γ ≡ dT/dz (ch01); met convention shown alongside | ch04 → Ch. 11, 13 | `internal_froude_number`, `richardson_number`, `gradient_richardson_number` |
| C_p ⚠️ | pressure coefficient (p − p_∞)/(½ρU²) (4.106) | – | **C_p also = specific heat at constant pressure** | ch04 → Ch. 6, 14 | `pressure_coefficient` |
| C_D, C_L | drag and lift coefficients F/(½ρU²A) (4.107)–(4.108) | – | A = frontal, plan or wetted area (say which) | ch04 → | `drag_coefficient`, `lift_coefficient`, `reference_area`, `sphere_drag_coefficient` (Morrison) |
| Ec, Pr | Eckert U²/(C_pδT), Prandtl ν/κ = μC_p/k | – | air 0.71, water ≈ 7 at 20 °C; monatomic 2/3 | ch04 → | `eckert_number`, `prandtl_number`, `eucken_prandtl`, `prandtl_of` |
| We, Bo, Ca | ρU²l/σ, ρgl²/σ, μU/σ (force ratios, not per volume) | – | Ca = We/Re | ch04 → | `weber_number`, `bond_number`, `capillary_number` |
| λ ⚠️ (model scale) | l_m/l_p (Ex. 4.8 uses 1/25) | – | Froude matching U_m = U_p√λ; wave drag × λ⁻³; Re ratio λ^{3/2} | ch04 | `froude_scaled_speed`, `model_prototype`, `ship_drag_extrapolation` |
| **Vortex lines, tubes and the basic vortices (ch05 §5.1)** | | | | | |
| ω ⚠️ | vorticity ∇×u; in (5.1) u_θ = ωr/2 the fluid turns at ω/2 | 1/s | counterclockwise positive in the plane; **ch03 `omega0` = rotation rate ω/2; ω = angular frequency in ch03 Ex. 3.1 and Ch. 7** | ch02 → | `omega`; `solid_body_from_vorticity(r, omega)` |
| Ω ⚠️ (tank) | rotation rate of a tank = ω/2 | rad/s | **Ω also the frame rotation (ch04, §5.6) and ch04's imposed frequency scale** | ch05 | `Omega_tank` in `rotating_tank_free_surface` |
| dx/ω_x = dy/ω_y = dz/ω_z | vortex line (5.3); parametric dx/ds = ω/\|ω\| | – | traced both ways from x₀ | ch05 → | `vortex_line(omega, x0, s_max, both)` |
| Γ (tube) | tube strength ∮u·dx = ∫ω·n dA, the same at every section (5.4) | m²/s | loop counterclockwise about the section normal | ch05 → | `vortex_tube_strength` → `TubeStrength`, `tube_flux_budget` → `TubeFlux(lower, side, upper, total)` |
| a ⚠️ | core radius (Rankine, rotating cylinder, ring core, Gaussian tube a₀) | m | **also ring radius R vs core a; Hill's sphere radius a** | ch05 → | `a`, `a0` |
| p_o, p_∞ | pressure on the axis at z = 0 (tank); far-field pressure (line vortex) | Pa | gauge-like references, default 0 | ch05 | `p_o`, `p_inf` |
| σ_rθ ⚠️ | viscous shear stress in polar coordinates = μ[(1/r)∂u_r/∂θ + r∂(u_θ/r)/∂r]; line vortex −μΓ/πr² | Pa | net force per volume (1/r²)∂(r²σ_rθ)/∂r (the r² trap); **σ = surface tension (ch01), core radius (ch03, `sigma_core` here)** | ch05 | `polar_viscous_stress`, `line_vortex_viscous_stress`, `polar_net_viscous_force` |
| torque per length | 2πr²σ_rθ = −2μΓ at every r outside a rotating cylinder | N m/m | on the fluid inside radius r | ch05 | `torque_per_length`, `edge_line_force` (jump −μΓ/πa² at r = a) |
| B (vortex) | u_θ²/2 + gz + p/ρ; grows as ω²r²/4 in the tank, uniform for the line vortex | m²/s² | B − B(0) | ch04 → | `bernoulli_across_vortex(kind, r)` |
| **Circulation and Kelvin (ch05 §5.2–5.3)** | | | | | |
| C, x(s, t) | material loop with fixed particle labels s ∈ [0, 1) | m | points advected in one `solve_ivp` call (DOP853) | ch05 → | `material_loop(u, pts0, t_eval)` → (n_t, 3, N) |
| A_vec | vector area ½∮x × dx (planar loop: area × normal) | m² | right-hand rule with the loop direction | ch05 → | `loop_vector_area(pts)` |
| 𝒫 | barotropic pressure function ∫dp/ρ(p) | m²/s² | ∮d𝒫 = 0 only if ρ = ρ(p) | ch04 → | `core.bernoulli.pressure_function` |
| δ ⚠️ | thickness of the smoothed (tanh) density interface in the lock exchange | m | rate ∝ 1/δ | ch05 | `delta` in `lock_exchange_*` |
| ρ₁, ρ₂ | light and heavy densities | kg/m³ | heavy fluid on the left ⇒ counterclockwise | ch05 | `rho1`, `rho2` |
| θ ⚠️ (tilt) | angle of the isopycnals from the isobars in E4 | rad | ∇ρ = \|∇ρ\|(sin θ, −cos θ); **θ also polar angle, segment end angles** | ch05 | `tilt` in `baroclinic_element_scenario` |
| x_G, I_G | centre of mass of the disc (offset R²∇ρ/4ρ₀) and its moment of inertia ½πρ₀R⁴ | m, kg m²/m | torque about G | ch05 | `pressure_torque_on_element`, `baroclinic_element_scenario` |
| **Vorticity equation (ch05 §5.4, §5.6)** | | | | | |
| (ω·∇)u | stretching + tilting term | 1/s² | = Gω with G[i, j] = ∂u_i/∂x_j | ch05 → | `vorticity_terms`, `stretching_tilting_split(omega, G)` |
| 2Ω, ω + 2Ω | planetary and absolute vorticity | 1/s | planetary term 2(Ω·∇)u = 2Ω∂u/∂z for Ω = Ωe_z | ch03 → | `absolute_vorticity`, `planetary_vorticity_terms(G, Omega)` |
| ∇ρ × ∇p/ρ² | baroclinic source (5.28) | 1/s² | order ∇ρ × ∇p | ch05 → | `baroclinic_term(rho, p, x)`, `baroclinic_rate_2d` |
| u_{i,j} | comma notation ∂u_i/∂x_j | – | free index renamed n → i at the end of D15 | ch02 → | (sympy in D14, D15) |
| α ⚠️ | strain rate of the Burgers / uniform-strain flows (u_z = αz) | 1/s | book's α = 2 × Wikipedia's (Burgers); **α also thermal expansion (ch01), E4 cube spin (ch04)** | ch05 | `alpha` in `burgers_vortex`, `strain_preset(name, rate)` |
| s ⚠️ | shear rate of the `shear_tilt` preset (w = s x) | 1/s | **s also arc length along a vortex line** | ch05 | `rate` in `strain_preset("shear_tilt")` |
| γ (sheet, diffusing) | strength of a diffusing vortex sheet (Ex. 5.6) | m/s | ω = γ/(2√(πνt))e^{−y²/4νt}, u(∞) = −γ/2 | ch05 | `diffusing_vortex_sheet(y, t, gamma, nu)` |
| A (Hill) ⚠️ | Hill's vortex constant, ω = AR e_φ; U = 2Aa²/15 | 1/(m s) | **A also area** | ch05 | `hill_spherical_vortex(R, z, A, a)` |
| **Natural coordinates (ch05 (5.31)–(5.32))** | | | | | |
| e_s, e_n, e_m | unit tangent to the vortex line, normal **away** from the centre of curvature (= −Frenet N), e_m = e_s × e_n | – | see the sign trap | ch05 | `helix_frame(s, a, c)`, `frenet_frame(curve, s)` |
| κ, τ ⚠️ | curvature, torsion of a vortex line | 1/m | helix κ = a/(a² + c²), τ = c/(a² + c²); **τ = stress elsewhere, κ = thermal diffusivity** | ch05 | `curvature`, `torsion` keys |
| ω_s, ω_n, ω_m | vorticity components in the natural frame (ω_n = ω_m = 0 at the point) | 1/s | D/Dt of each = ω ∂u_{s,n,m}/∂s | ch05 | `stretching_tilting_split` → `stretching`, `tilting` |
| **Rotating frame and columns (ch05 (5.33))** | | | | | |
| ζ ⚠️ | relative vorticity ω_z of a column | 1/s | cyclonic > 0 in the NH; **ζ = parcel displacement (ch01), cap/meniscus height (ch04)** | ch04 (gloss) → Ch. 13 | `zeta`, `zeta0` in `column_relative_vorticity` |
| f | local planetary vorticity 2Ω sin φ | 1/s | NH > 0 | ch04 → | `coriolis_parameter`, `f=` |
| h ⚠️ | column height | m | (ζ + f)/h conserved; **h also vortex spacing (§5.7), wall distance** | ch05 → Ch. 13 | `h`, `h0` |
| φ | latitude | ° at interfaces | `*_deg` names; radians inside | ch04 → | `lat_deg`, `lat0_deg`, `lat1_deg` |
| **Biot–Savart and filaments (ch05 §5.5)** | | | | | |
| G(x, x′) | Green's function of ∇² in 3-D, −1/(4π\|x − x′\|) | 1/m | ∇²G = δ | ch05 → | `poisson_green_3d(x, xp)` |
| x, x′ | field point and source point | m | ∇′ acts on x′ | ch05 → | `x`, `xp`, `nodes` |
| e_ω, dl | unit vector along a filament, element length | –, m | du = (Γdl/4π)e_ω × (x − x′)/\|x − x′\|³ | ch05 → | `filament_velocity(x, polyline, Gamma, closed)` |
| d, θ_a, θ_b ⚠️ | perpendicular distance to a segment; angles at its ends | m, rad | (Γ/4πd)(cos θ_a − cos θ_b) | ch05 → | `segment_speed(d, theta_a, theta_b, Gamma)` |
| ε (kernel) ⚠️ | smoothing length r² → r² + ε² | m | ≈ grid spacing inside a core; **ε = Levi-Civita / dissipation elsewhere** | ch05 | `eps` in `biot_savart_*`, `point_vortex_velocity` |
| R, a (ring) | ring radius, core radius | m | self-speed Γ/4πR[ln(8R/a) − C], C = ¼ uniform, ½ hollow, 0.558 Gaussian (a = √(4νt)) | ch05 → Ch. 14 | `ring_self_velocity(R, a, Gamma, core)`, `RING_CORES` |
| m ⚠️ | elliptic-integral parameter k² | – | **not the modulus k; m = mass elsewhere** | ch05 | inside `ring_ring_velocity` |
| M ⚠️ | number of segments of a polygon filament | – | error ∝ 1/M² at a ring centre | ch05 | `M` in `filament_preset` |
| **Point vortices, images, sheets (ch05 §5.7–5.8)** | | | | | |
| x_k, Γ_k | point-vortex positions and strengths | m, m²/s | counterclockwise positive; self excluded | ch05 → | `point_vortex_velocity(x, xv, Gamma)`, `point_vortex_evolve` |
| h ⚠️ | vortex separation (pairs); distance to a wall; channel position | m | V = Γ/2πh (pair), Γ/4πh (wall) | ch05 | `vortex_pair(G1, G2, h)`, `vortex_near_wall_speed(Gamma, h)` |
| G ⚠️ | centre of vorticity ΣΓ_kx_k/ΣΓ_k | m | **not a stagnation point** (Fig. 5.11 caption); **G = velocity gradient, g in ch04** | ch05 | `centre_of_vorticity` |
| P, I, H ⚠️ | linear impulse ΣΓx, angular impulse ΣΓ\|x\|², Kirchhoff energy −(1/2π)ΣΓ_jΓ_k ln\|x_j − x_k\| | m³/s, m⁴/s, m⁴/s² | conserved; **H also channel width** | ch05 | `point_vortex_invariants` |
| H ⚠️ (channel) | channel width | m | drift (Γ/4H)cot(πh/H) | ch05 | `channel_image_velocity(h, H, Gamma)` |
| u₁, u₂ | tangential velocity above / below a sheet | m/s | γ = u₂ − u₁ | ch05 | `vortex_sheet_strength(u_above, u_below)` |
| N ⚠️ | number of filaments in a discrete sheet | – | L1 error ∝ 1/N; **N² = buoyancy frequency elsewhere** | ch05 | `discrete_sheet_u(x, y, gamma, N)` |
| **Ideal flow (ch06)** | | | | | |
| φ | velocity potential, u = ∇φ; multivalued round a vortex (jumps by Γ across the cut) | m²/s (2-D and 3-D) | harmonic, ∇²φ = q (6.11) with sources | ch03 → | `Flow.phi`, `velocity_potential`, `AxisymFlow.phi` |
| ψ | plane stream function, u = ∂ψ/∂y, v = −∂ψ/∂x; ω_z = −∇²ψ (6.4) | m²/s | Δψ = volume flux per depth | ch04 → | `Flow.psi`, `stream_function` |
| ψ ⚠️ (Stokes) | axisymmetric stream function, u_R = −(1/R)∂ψ/∂z, u_z = (1/R)∂ψ/∂R (6.75) | **m³/s** | flux between surfaces 2πΔψ (6.78); field equation (6.77) ≠ ∇²ψ | ch04 (gloss) → ch06 | `AxisymFlow.psi(R, z)`, `stokes_operator_residual` |
| w ⚠️ | complex potential φ + iψ (6.42) | m²/s | analytic; **w = z-velocity in §6.9 and other chapters** | ch06 | `ComplexFlow.w(z)` |
| z, ζ ⚠️ | complex coordinate x + iy = re^{iθ} (6.43); ζ = Zhukhovsky circle plane, z = ζ + b²/ζ | m | outside root \|ζ\| ≥ b is physical; **ζ = relative vorticity (ch05), displacement (ch01)** | ch06 | `z` (complex arrays), `joukowski(zeta, b)`, `joukowski_inverse(z, b, branch)` |
| dw/dz | complex velocity u − iv (6.45) | m/s | velocity vector = conj(dw/dz) | ch06 | `ComplexFlow.dwdz`, `.complex_velocity` |
| U, α | free-stream speed and its direction (from +x) | m/s, rad | `Uniform(U, V)` takes the two components (U cos α, U sin α); tilted stream in E5: D − iL = −iρUΓe^{−iα} | ch03 → | `U`, `alpha` |
| m ⚠️ | 2-D source strength (volume flux per depth) | m²/s | sink m < 0; φ = (m/2π) ln r (6.15); **m = elliptic parameter (ch05), mass** | ch06 | `Source(m, z0)` |
| Q ⚠️ | 3-D point-source strength (volume flux) | m³/s | φ = −Q/4πr; **Q = flow rate of Example 6.2 (ours Q = 1)** | ch06 | `PointSource3D(Q, z0)`, `example_6_2(Q=)` |
| Γ ⚠️ | vortex circulation — **counterclockwise in code**, clockwise in half the book's equations (see the trap) | m²/s | L = ρUΓ_cw | ch03 → | `Vortex(Gamma)`, `Gamma_cw=`, `Gamma_ccw=`, `gamma_ccw_from` |
| d ⚠️ | doublet (dipole) vector Σx_i m_i, **from sink to source**; scalar d of (6.49) = dipole −d e_x | m³/s (2-D), m⁴/s (3-D) | cylinder d = −2πUa² e_x; sphere d = 2πa³U; **d = segment distance (ch05)** | ch06 | `Doublet(d_vec)`, `Doublet.from_book_scalar(d)`, `Doublet3D` |
| A, n ⚠️ | corner flow w = Azⁿ (6.46): walls θ = 0, π/n; n = π/α for a wedge of angle α; speed ∝ r^{n−1} | A: m^{2−n}/s | n < 1 re-entrant (infinite speed), n > 1 stagnation; **n also unit normal, A also area** | ch06 | `Corner(A, n, cut_angle, rotate)`, `corner_speed_exponent` |
| a ⚠️ | cylinder/sphere radius; half-body stagnation distance m/2πU; airship line-sink length; Zhukhovsky circle radius (> b) | m | | ch06 | `a` |
| b ⚠️ | Zhukhovsky constant (slit half-length 2b, foci ±2b) | m | **b = CV velocity (ch03/04), Hill constant** | ch06 | `b` in `core.conformal` |
| C_p | pressure coefficient (p − p∞)/(½ρU²) = 1 − \|u\|²/U² (6.32) | – | cylinder 1 − 4 sin²θ ∈ [−3, 1]; sphere 1 − (9/4) sin²θ | ch04 → | `pressure_coefficient`, `Flow.cp`, `*_surface_cp` |
| D, L | drag and lift per unit depth **on the body** | N/m | D − iL = (iρ/2)∮(dw/dz)²dz (6.60); F of (6.54) is on the fluid | ch06 | `blasius_force` → `BlasiusForce(D, L, …)`, `lift_per_span` |
| c_k | Laurent coefficients of dw/dz outside the body | m^{1−k}/s … | c₀ = U, c₋₁ = iΓ_cw/2π (residue), c₋₂ = −Ua² (cylinder) | ch06 | `laurent_coefficients(flow, R, …)` |
| k ⚠️ | half-body ratio sin θ/(π − θ) in C_p = −(2k cos θ + k²); line-sink strength per length; axial segment strengths k_n | –, m²/s | **k = wavenumber (Ch. 7), thermal conductivity (ch01)** | ch06 | `half_body_surface_cp`; `airship(U, Q, a)` (k = Q/a); `axial_singularity_solve` → k |
| ψ_{i,j}, Δx, Δy | grid values and spacings of the Laplace grid (6.70)–(6.73) | m²/s, m | `mask[j, i]` True = unknown, x on the last axis; average rule needs Δx = Δy | ch06 | `core.laplace_solvers` |
| ρ_J, ρ_GS, ω_SOR ⚠️ | spectral radii of the Jacobi / Gauss–Seidel sweep; SOR factor | – | 4-point grid ρ_J = ½, ρ_GS = ¼, ω_opt = 1.0718; **ρ = density, ω = vorticity elsewhere** | ch06 | `jacobi_spectral_radius`, `optimal_sor_omega` |
| R, z ⚠️ (§6.8) | cylindrical radius from the symmetry axis and the **horizontal** axial coordinate along the stream | m | z up in ch01–ch05 | ch06 | `AxisymFlow.psi(R, z)` |
| ξ ⚠️ | §6.8: axial source coordinate along a line sink; §6.9: vector x − x_s from the sphere's centre | m | ∂/∂x_s = −∇ for fields of ξ | ch06 | `xi_nodes`; `e_xi` |
| x_s, u_s | sphere centre and velocity (u_s = dx_s/dt) | m, m/s | arbitrary path | ch06 | `moving_sphere_*(x, xs, us, …)` |
| M ⚠️ | added mass (2π/3)ρa³ = ½ displaced mass (sphere); ρπa² per depth (cylinder) | kg, kg/m | F_s = −M du_s/dt (6.108); **M = segment count (ch05), Mach number** | ch06 | `added_mass_sphere(a, rho)`, `cylinder_added_mass` |
| θ_s ⚠️ | polar angle from the sphere's velocity in (6.106) | rad | = π − θ of (6.91); see the three-origins trap | ch06 | (inside `moving_sphere_surface_pressure`) |
| N ⚠️ | number of axial segments / panels | – | even N for fore–aft symmetric bodies (odd N singular); cap 40 | ch06 | `axial_singularity_solve(N=)`, `source_panels` |
| λ_j, S_j | source-panel strength and panel length | m/s, m | Σλ_jS_j = 0 for a closed body; self term ½ | ch06 | `source_panels`, `panel_geometry` |
| **Gravity waves (ch07)** | | | | | |
| x, z ⚠️ | horizontal (along propagation) and vertical coordinates | m | **z up**, still surface z = 0, flat bottom z = −H; §7.7 interface at z = −H | ch07 → | `x`, `z` |
| a | wave amplitude (crest height above the mean level) | m | crest +a, trough −a | ch07 → | `a` |
| k ⚠️ | wavenumber 2π/λ | rad/m | **sign gives the direction**; \|k\| in every ω(k); **k = thermal conductivity (ch01), half-body ratio and line-sink strength (ch06)** | ch07 → | `k`, `k0` |
| λ ⚠️ | wavelength | m | **λ = bulk-viscosity / model scale / eigenvalue elsewhere** | ch07 → | `lam` |
| ω ⚠️ | angular frequency 2π/T (intrinsic) | rad/s | ω ≥ 0; **vorticity in ch02–ch06** (see the trap) | ch07 → | `omega` |
| ω₀ | frequency seen by a fixed probe in a current U | rad/s | ω₀ = ω + U·K (7.9) | ch07 → | `doppler_frequency(omega, U, K)` |
| T ⚠️, ν ⚠️ | period 2π/ω; cyclic frequency 1/T | s, Hz | **T = temperature, ν = kinematic viscosity elsewhere** | ch07 | `T`, `wave_parameters(nu=)` |
| c | phase speed ω/k = λν | m/s | `phase_speed` returns \|c\| ≥ 0 | ch07 → | `phase_speed(k, H, g, sigma, rho)` |
| c_g | group velocity dω/dk | m/s | **signed** (sgn k); c/2 deep, c shallow, 3c/2 capillary | ch07 → | `group_velocity`, `group_velocity_numeric` |
| K = (k, l, m), K, e_K | wavenumber vector, its magnitude, unit vector | rad/m | c = (ω/K)e_K (7.8); **m = vertical wavenumber here (mass, source strength, elliptic parameter elsewhere)** | ch07 → | `K`, `plane_wave(X, K, omega)`, `phase_velocity_vector` |
| c_x, c_y, c_z | trace velocities ω/k, ω/l, ω/m | m/s | each ≥ c; **not components of c** (1/c² = Σ1/c_i²) | ch07 | `trace_velocities(K, omega)` → tuple |
| θ ⚠️ (phase) | local phase θ(x, t), k = ∂θ/∂x, ω = −∂θ/∂t (7.72) | rad | see the ζ/θ trap | ch07 | `theta_fn` in `local_wavenumber_frequency` |
| H ⚠️ | still-water depth; §7.7 upper-layer thickness | m | `H = np.inf` = deep water; **H = scale height (ch01), channel width (ch05)** | ch07 → | `H` |
| kH | depth parameter | – | deep kH > 2 (H > 0.318λ, error 1.815 %), shallow H < 0.07λ (kH < 0.44, error 3.039 %) | ch07 → | `depth_regime(k, H)` → label + errors |
| η ⚠️ | surface elevation η(x, t) | m | **ch04's η was the level-set function** (see the trap) | ch07 → | `eta` |
| f ⚠️ (surface) | level-set function z − η | m | n = ∇f/\|∇f\| (7.14); **f = Coriolis parameter, Helmholtz free energy, traction elsewhere** | ch07 | inside `surface_normal(eta_x)` |
| U_s | velocity of the surface point at fixed x, η_t e_z (7.15) | m/s | only its normal part matters | ch07 | `surface_velocity(eta_t)` |
| φ, ψ | velocity potential; stream function with **u = ∂ψ/∂z, w = −∂ψ/∂x** | m²/s | ch04's planar convention with y → z; GFD's opposite sign is the Ch. 13 trap | ch03/ch04 → | `wave_fields(...)["phi"]`, `["psi"]` |
| p, p′ ⚠️ | gauge pressure; perturbation p + ρgz (§7.2) or p − p̄(z) (§7.8) | Pa | see the p′ trap | ch07 → | `linear_bernoulli_pressure`, `wave_fields(...)["p_prime"]` |
| ξ, ζ ⚠️ | horizontal and vertical particle excursions from the mean position (x₀, z₀) | m | orbits clockwise for a right-going wave; **ζ = interface displacement in §7.7** | ch07 → | `orbit_linear` → (xi, zeta) |
| x₀, z₀ | mean position of a particle (integration constants zero by definition) | m | `particle_path(start="mean")` releases at the mean depth | ch07 | `x0`, `z0` |
| A, B ⚠️ (orbit) | semi-axes a cosh k(z₀ + H)/sinh kH and a sinh k(z₀ + H)/sinh kH | m | focal distance a/sinh kH at every depth; **A, B, C also the constants of (7.23)–(7.25) and (7.104)–(7.109)** | ch07 | `orbit_semi_axes(z0, a, k, H)` |
| E, E_k, E_p ⚠️ | wave energy, kinetic and potential parts | J/m² (surface), J/m³ (internal), J/kg (jump) | see the energy trap; E = ½ρga² (7.42) | ch07 → | `wave_energy`, `wave_energy_density`, `internal_wave_energy` |
| F ⚠️ | energy flux E c_g | W/m (surface), W/m² (internal) | vector for internal waves; **F = force elsewhere** | ch07 → | `energy_flux`, `internal_wave_energy(...)["F"]` |
| σ ⚠️ | surface tension | N/m | `omega_capillary_gravity` default 0.0727, `phase_speed` default 0; teaching value 0.07274 (IAPWS 20 °C) with ρ = 998.2; **σ_x = packet width in `gaussian_packet`, core radius (ch03)** | ch01 → | `sigma` |
| R ⚠️ (curvature) | radius of curvature of the surface, 1/R = η_xx/(1 + η_x²)^{3/2} | m | centre below a crest (η_xx < 0) ⇒ liquid pressure higher, p = −ση_xx (7.54) | ch07 | `curvature(eta_x, eta_xx, linear=)`, `capillary_surface_pressure` |
| c_min, λ_m, k_m | minimum phase speed (4gσ/ρ)^{1/4} and where it occurs 2π√(σ/ρg) | m/s, m, rad/m | c_g = c there; clean water 23.12 cm/s at 1.712 cm | ch07 | `capillary_minimum(sigma, rho, g)` |
| c_g,min | minimum group speed (σk²/ρg = 2/√3 − 1) | m/s | 17.76 cm/s at 4.354 cm (computed, print 4 s.f.) | ch07 | `min_group_velocity` |
| L, n ⚠️, b ⚠️ | basin length, seiche mode number (n = 0 fundamental), basin width | m, –, m | λ = 2L/(n + 1) (7.64); **n = unit normal, b = Zhukhovsky constant / interface amplitude elsewhere** | ch07 | `seiche_modes(L, H, n)`, `basin_modes(L, b, H, m, n)` |
| Δk, Δω, A(k), δk, σ_x, ω″ | beat differences; packet spectrum, its width, the packet's length; dispersion d²ω/dk² | rad/m, rad/s, m | order-2 packets keep ω″ (spreading) | ch07 | `beat_wave(x, t, k1, k2)`, `gaussian_packet(..., sigma_x, order)` |
| α ⚠️ (ray) | angle between k and the shore normal | rad (° at interfaces) | \|k\| sin α = const on straight contours (Snell); **α = thermal expansion, free-stream angle elsewhere** | ch07 | `snell_ray_plane_beach(x, alpha0, x0, T, slope)`, `refraction_state(T, alpha0, H0, H)` |
| s ⚠️ | beach slope dH/dx | – | 1:50 in the ray tests | ch07 | `slope` |
| H₁, H₂, u₁, u₂, Q ⚠️, Fr₁ | upstream/downstream depths and speeds of a jump, discharge per width u₁H₁, upstream Froude u₁/√(gH₁) | m, m/s, m²/s, – | H₂/H₁ = ½(−1 + √(1 + 8Fr₁²)); Fr₁ ≥ 1 (second law); **Q = 3-D source strength [m³/s] in ch06** | ch07 | `hydraulic_jump(H1, u1=, Fr1=)`, `jump_state(H1, Fr1)` |
| ū_L | Stokes drift (Lagrangian mean velocity) | m/s | a²ωke^{2kz₀} deep (7.85): decays as e^{2kz}; Eulerian mean 0 | ch07 → Ch. 13 | `stokes_drift(z0, a, k, H)`, `stokes_drift_numeric` |
| γ ⚠️ (Stokes) | amplitude-dispersion coefficient in c² = (g/k)(1 + γk²a²) | – | γ = 1 (consistent 3rd order); 3/8 from a truncated set-up; **γ = sheet strength (ch05), heat-capacity ratio (ch01)** | ch07 | `stokes_expansion_sympy` |
| c₀, Ur | shallow-water speed √(gH); Ursell number aλ²/H³ | m/s, – | KdV linear speed c₀(1 − (kH)²/6); crossover 4π²/9 | ch07 | `kdv_linear_phase_speed`, `ursell_number(a, lam, H)` |
| m ⚠️ (cnoidal) | elliptic parameter of a cnoidal wave (m → 1: solitary wave) | – | as ch05's `ellipk(m)` (m = k²) | ch05 → | `cnoidal_wave(x, t, H, height, m)` |
| ρ₁, ρ₂ | upper (lighter) and lower (heavier) densities | kg/m³ | ρ₁ > ρ₂ ⇒ Rayleigh–Taylor (NaN + warning) | ch05 → | `rho1`, `rho2` |
| ε² ⚠️ | density ratio (ρ₂ − ρ₁)/(ρ₂ + ρ₁) | – | ω = ε√(gk) (7.95); **ε = dissipation rate (ch04), Levi-Civita, kernel smoothing (ch05)** | ch07 | `eps2_density(rho1, rho2)` |
| ζ ⚠️ (interface), b ⚠️ | interface displacement; its complex amplitude in the two-layer problem | m | barotropic b = ae^{−kH}; baroclinic η/ζ = −((ρ₂ − ρ₁)/ρ₁)e^{−kH} < 0 (7.114) | ch07 | `two_layer_modes(k, H, rho1, rho2)` |
| g′ ⚠️ | reduced gravity g(ρ₂ − ρ₁)/ρ₂ (7.117) | m/s² | **ρ₂ in the denominator; ch04's `reduced_gravity` uses ρ₁** (see the trap) | ch04 → | `reduced_gravity_book(rho1, rho2, g, ref="lower")` |
| N, ρ₀, ρ̄(z), ρ′ | buoyancy frequency, reference density, base density, perturbation | rad/s, kg/m³ | N² = −(g/ρ₀)dρ̄/dz; ρ′ = N²ρ₀ζ/g | ch01 → | `N`, `rho0` |
| ∇_H² | horizontal Laplacian ∂²/∂x² + ∂²/∂y² | 1/m² | commutes with ∂/∂t, ∂/∂z and N(z) | ch07 | (sympy in `boussinesq_linear_sympy`) |
| θ ⚠️ (K angle) | angle of K above the horizontal, cos θ = \|k\|/K | rad | ω = N cos θ (7.139); = the beam's angle from the vertical | ch07 → | `beam_angle(omega, N)`, `internal_wave_omega(k, m, N, l)` |
| ŵ ⚠️ | complex amplitude of w in the plane internal wave | m/s | û = −mŵ/k, p̂ = −ωmρ₀ŵ/k², ρ̂ = iN²ρ₀ŵ/(ωg) (7.153); **`w0` = ch01's initial parcel velocity** | ch07 | `w0` in `internal_wave_fields`, `internal_wave_energy` |
| i, Re{} | imaginary unit; real part restored at the end | – | see the complex-notation trap | ch06 → | `real_field(amp, phase)` |

## Coordinate and sign conventions per chapter
| Chapter | Axes (which is "up") | Origin / reference level | Stress / pressure sign | Reference scales (L, U, T) | Dimensional or non-dimensional code |
|---|---|---|---|---|---|
| ch01 | z up (§1.7, §1.10); Couette y from fixed wall (0) to moving plate (h) | p0 at z = 0; θ reference p_ref = 1000 hPa; parcel rest height z_o | pressure absolute and isotropic, acts along the inward normal; τ_xy = +μ ∂u/∂y signed; lapse rate Kundu dT/dz; first law q in / w on | none fixed (Π groups carry their own; E2 uses the clock h²/ν) | dimensional SI throughout; only Π groups are dimensionless |
| ch02 | right-handed x₁x₂x₃ (no preferred "up"); rotated frame shares the origin; angles counterclockwise about e₃; grid arrays `[k, j, i]` = (z, y, x), components and `h` in (x, y, z) | origin of both frames; boxes/loops centred at x0 | τ_ij tensile positive, +e_i face → +e_j; traction f = n·τ (first index); ∇·τ on the second index; pressure τ = −pδ; passive C (x' = Cᵀx); book A:B = A_ij B_ji; R = G − Gᵀ ↔ ω = ∇×u, A = ½R ↔ ½∇×u; Stokes n_c into A, t counterclockwise about n; outward n on closed surfaces | none (pure mathematics) | dimensional where physical (Pa, 1/s, m); most results unit-agnostic |
| ch03 | right-handed Cartesian (no preferred "up"); plane polar (r, θ from +x), cylindrical (R, φ, z), spherical (r, θ from +z, φ); angles and rotation counterclockwise positive (shear spins clockwise: ω₃ = −γ); field callables `u(x, t)` with coordinates on axis 0; Galilean frame O′ at constant U (x = x′ + Ut + x′_o, u′ = u − U); rotating frame u = Ω × x + u′ at the coinciding instant | cylinder centre at the origin at t = 0 for every observer (E3); Ex. 3.1 port at the origin; vortices centred at the origin; RTT shapes with explicit reference geometry (E7 interval [1, 3] + ȧt, ḃt; ellipse a = 2 + ȧt, b = 1 + ḃt) | R = G − Gᵀ (no ½), ω = ∇×u, spin ½ω; γ = 2S₁₂; RTT outward n, signed b·n; Leibniz lower term subtracted | none (E6 draws in r/σ; its real-vortex modes use metres) | dimensional SI throughout |
| ch04 | right-handed Cartesian, **z up**, g = −g e_z, Φ = gz (4.18); cylindrical (R, φ, z) for rotating flows and Ex. 4.5; noninertial frame O′ translating at U(t) and rotating at Ω(t) with basis e′_i; NH Ω_z > 0; latitude φ | CVs with explicit geometry (E1 boxes around the wake, bore, jet, rocket, balloon); hydrostatic base state p_s(z), ρ_s(z) for Boussinesq; free surface η = 0 | τ = −pδ + σ, tensile positive; traction f_j = n_iτ_ij and Cauchy's divergence on the **first** index; outward n, signed (u − b)·n; drag on the body +x, on the fluid −F_D; acceleration terms +2Ω × u′, +Ω × (Ω × x′) (4.43) vs forces −2Ω × u′, −Ω × (Ω × x′) (4.45); Stokes assumption μ_v = 0; 2-D ψ: u = ∂ψ/∂y; primes: rotating frame (§4.7), perturbation (§4.9), dummy (4.67) | `core.similarity.Scales` holds one set per use with its time scale (1/Ω (4.100) or l/U (4.109)) and pressure scale (ρU², μU/l or ρgl); Ro = U/(2Ωl) forward pointer | dimensional SI in functions; non-dimensional via `Scales` and the `nondimensional_*` coefficient routines; g default 9.81 in `ch04`, 9.80665 in `core` |
| ch05 | right-handed Cartesian, **z up**, g = 9.81 (`core.thermo.G_BOOK`); plane polar (r, θ) and cylindrical (R, φ, z) for vortices and rings; material loops with fixed labels s ∈ [0, 1); rotating frame Ω = Ωe_z (NH > 0); natural frame (e_s, e_n = −Frenet N, e_m) on vortex lines; wall at y = 0 with the fluid above; sheet along x with u₁ above, u₂ below | vortices and rings centred on the axis; tank p_o on the axis at z = 0; line vortex p_∞ far away; lock-exchange interface at x = 0 (heavy on the left); column undisturbed depth h₀ | counterclockwise-positive ω_z, Γ, point-vortex strengths and sheet strength γ = u₂ − u₁; ∇ρ × ∇p order; (5.14) with +1/(4π); ω = vorticity (tank turns at ω/2); σ_rθ viscous stress with the polar metric; Γ_a = Γ + 2Ω·A_vec | none fixed (every function dimensional; the E-series use lab or planetary numbers directly) | dimensional SI throughout; latitudes in degrees only at `*_deg` interfaces |
| ch06 | 2-D flows in (x, y) with θ from +x (downstream); complex z = x + iy, ζ the Zhukhovsky circle plane; §6.8 cylindrical (R, φ, z) with **z horizontal along the stream**, spherical (r, θ from +z); §6.9 ξ = x − x_s; Laplace grids `mask[j, i]` with x on the last axis | bodies centred at the origin; half-body source at 0 (stagnation at −m/2πU); p∞ far upstream; ideal-flow p measured from hydrostatic; Example 6.1 wall at x = 0 with the vortex at (h, 0) at t = 0 | Γ counterclockwise in code, **clockwise** in (6.36)–(6.40), (6.52), (6.61)–(6.62), (6.68), Ex. 6.1 (`Gamma_cw=`); D, L on the body, (6.54) F on the fluid; ccw contour, outward n; 2-D dipole from sink to source; (6.82) = r × App. B | none fixed (U, a or the body length set the scale in each function; C_p is the only non-dimensional output used throughout) | dimensional SI throughout; ρ default 1.2 (2-D forces) or 1000 (Ex. 6.1, §6.9) — pass it explicitly |
| ch07 | 2-D waves in the (x, z) plane, **z up**, x along propagation; still surface z = 0, flat bottom z = −H (H = ∞ deep); §7.7 origin at the mean free surface, interface at z = −H; §7.8 (x, y, z) with z up, K = (k, l, m); rays in (x, y) with α from the shore normal | still-water level; mean particle position (x₀, z₀) | p **gauge**; p′ = p + ρgz (§7.2) or p − p̄(z) (§7.8); ψ with u = ∂ψ/∂z; ω ≥ 0, direction in sgn k; c_g signed; clockwise orbits for +x waves; sheet γ = u_below − u_above | g = 9.81 (`G_BOOK`) in ch07 and `core.waves`; ρ = 1000; clean water σ = 0.07274 N/m, ρ = 998.2 in teaching numbers; energy per area (surface), per volume (internal) | dimensional SI throughout; complex amplitudes (Re dropped) from §7.7; the KdV solver is dimensional |
