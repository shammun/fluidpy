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

**⚠️ Pressure-gradient sign (ch08 vs ch04).** The book (and every ch08 function) passes **dp/dx** — favourable when
negative, driving +x flow; ch04's `plane_poiseuille(y, G, h, mu)` and `exact_solution(..., G=)` take **G = −dp/dx**. ch08
functions accept both: `channel_flow(y, h, U, dpdx=…)` or `channel_flow(..., G=…)` (G wins the sign flip internally);
pipe functions take `dpdz` (or `G = −dp/dz`). A planted "G passed as dp/dx" flips the parabola and fails a test. Say the
sign of the number every time: "dp/dx = −0.5 Pa/m (favourable)".

**⚠️ η = y/√(νt) in ch08, y/(2√(νt)) in ch04 and on the book's own figures.** (8.25) defines η = y/√(νt), so
u/U = 1 − erf(η/2) and δ₉₉ = 3.643√(νt) (η₉₉ = 3.643); ch04's register row "η (Stokes)" and the book's Figs. 8.13–8.14
plot y/(2√(νt)), where the 1 % point is at 1.82. `similarity_variable(y, t, nu, half=False)` returns the ch08 η;
`half=True` the figures' axis. Label every axis with the variable it plots. (η was also ch04's level-set function and
ch07's surface elevation.)

**⚠️ ω twice inside ch08.** ω (ω_z, ω_φ) is the **vorticity** in §8.1, Example 8.5 and §8.6, but the **oscillation angular
frequency** [rad/s] in §8.5 (the ch07 meaning). Code: `stokes_second_problem(y, t, U, omega, nu)` takes the frequency;
vorticity functions are named `*_vorticity`, `vorticity_content`. **Stokes-layer depth**: the literature's e-folding depth
δ_e = √(2ν/ω) vs the book's "δ ~ 4√(ν/ω)" = 2√2 δ_e (5.91 % amplitude left); `stokes_layer` returns both (`delta_e`,
`delta_book`).

**⚠️ θ from the downstream axis in §8.6 (sphere); θ an azimuth in Example 8.6.** Stokes' and Oseen's solutions measure θ
from the +x axis, the direction of the stream U e_x: the **rear** stagnation point is θ = 0, the front θ = π (as ch06
(6.91); ch03's spherical θ is from +z; ch06 §6.9's moving sphere uses θ_s = π − θ). `frame="body"` = sphere at rest in
the stream; `frame="fluid"` = stream subtracted (sphere moving to −x, Figs. 8.19–8.20). Example 8.6 (line vortex) uses
plane polar (r, θ) with θ the azimuth. Always pass θ with its frame.

**⚠️ Four Reynolds numbers in ch08 (and a radius/diameter trap).** Pipe Re = Ud/ν (diameter, mean speed; laminar below
≈ 2000); lubrication Re_L = ρUL/μ (passage length — the governing group is **ε²Re_L**); generic Re = ρUL/μ (8.40);
**sphere Re = 2aU/ν (diameter)** in (8.52)–(8.53) and Oseen's C_D = (24/Re)(1 + 3Re/16); the literature (Wikipedia Oseen,
Proudman–Pearson) and the far-field ratio use the **radius** Re_a = ρUa/μ = Re/2, where Oseen's coefficient is 3/8 and the
ratio → ½Re_a r/a; Example 8.5 has Re_x = Ux/ν. Every drag function names its Re; `proudman_pearson_drag_coefficient`
converts.

**⚠️ Two pressure scales in the lubrication scaling (ch08 (8.14)).** The book scales p* = p/P_a with **absolute
atmospheric** P_a, which makes the bearing number Λ = μUL/(P_a h²) — 49.35 for our engine film, ≈ 10³ for thinner or
longer films, not "near unity"; the natural (viscous) scale is μUL/h² (5 MPa for the film), the one with Λ = 1.
`lubrication_term_magnitudes(p_scale="viscous" | "atm")`. At low Re (§8.6) the pressure scale is μU/L (dominant
balance), at high Re ρU² (4.100) — three pressure scales in one chapter.

**⚠️ Corrected book forms coded in ch08 (analysis §9 R6–R15), printed forms kept as options a test must fail.** (8.13b)
−∂p/∂**y** (printed ∂p/∂x; `lubrication_nondim_sympy(printed_8_13b=True)`) · (8.17a) and Example 8.2 with **ν** · (8.19)
with U₀(1 − y/h) (printed + U₀ gives u(h) = U_h + U₀; `lubrication_velocity(form="book")`) · Example 8.1 with (1 + αx/L)
integrands and a **squared** final denominator (`slider_bearing(model="book")`) · channel V = Q/h (printed middle form
lacks 1/h) · ∫₀^∞ω dy = **+U** (printed −U) · Example 8.5 η₉₅ = ±2.772 (printed ±2.76) · rear minimum p − p∞ =
**−**3μU/2a (printed without sign; (8.50) itself is right) · Oseen's equation with **−**∂p/∂x_i · power into the fluid
−2πR₁σ_Rφu_φ > 0 (printed (2πR₁)τ_Rφu_φ < 0) · (8.44) is (E²)²ψ = 0, not ∇⁴ψ = 0 · cross-references "(9.63)" → (8.43),
"(9.68)" → (8.48), "(8.33) into (8.20)" → (8.35), Example 8.7's "Example 8.2" → 8.3, Example 8.2's "y" → z, §8.1 "μ the
kinematic viscosity" → dynamic.

**⚠️ Earlier chapters' functions take different arguments for the same ch08 flow.** ch05 `rotating_cylinder_flow(r, a,
omega)` takes the cylinder's **vorticity** ω = 2Ω₁ (not its rotation rate); ch05 `diffusing_vortex_sheet(y, t, gamma,
nu)` takes γ = u_below − u_above, so Example 8.5 (u = ±U) is γ = **−2U**, while ch08 `vortex_sheet_diffusion(y, t, U, nu)`
takes U; `core.vortices.gaussian_vortex(r, Gamma, sigma)` equals the Lamb–Oseen decay with σ = 2√(νt). Parity tests pass
the converted arguments (a factor-2 variant is caught).

**⚠️ `slider_bearing_state["inlet_backflow"]` means "recirculation at the wide end".** After the ch08 fix it is True for
α > 1 (wide end x = L) **and** for α < −½ (wide end x = 0; for U > 0 that is the exit, so "inlet" misleads); `backflow_x`
gives the station and `backflow_any` equals it. Criterion: wide/narrow gap ratio > 2 (h > 1.5h_m, h_m = 2(1 + α)h₀/(2 + α)).
In `slider_gap_velocity` y = 0 is the moving pad side; the state docstring's η runs from the floor.

**⚠️ λ twice inside ch09.** In §9.6 (C07) λ = (θ²/ν)dU_e/dx is the **Holstein–Bohlen/Thwaites parameter** (dimensionless;
λ < 0 adverse; separation at −0.090 with the book's table, −0.068148 with our exact Falkner–Skan closure). In §9.8 (C10, D14)
λ is the **eigenvalue** of the linearised Kármán street, λ = (πΓ/2a²)(±γ ± iσ) [1/s], growth rate Re λ (the notebook uses λ
after lesson round 3; the code docstring of `karman_street_spectrum` still calls it σ). λ was the wavelength in ch07, the
second viscosity coefficient in ch04 and a trial exponent (e^{λt}, R^λ) in ch01/ch08; Töpfer's scale factor in D06 and the
wall-jet gauge f → λf(λη) (D21) are a third and fourth λ. Say which λ in every sentence.

**⚠️ γ and σ in D14 are dimensionless street coefficients, not a growth rate.** γ = ½ − sech²(πb/a) (zero at the marginal
spacing; −0.307 at b/a = 0.15) and σ = sinh(πb/a)/cosh²(πb/a); the growth rate is (πΓ/2a²)|γ| (0.48 Γ/a² at b/a = 0.15).
Reading notes once wrote "growth rate γ = 0.48" (lesson round 2, fixed). γ is also ch05's sheet strength, ch01's c_p/c_v and
ch08's ansatz variable; σ is ch07's surface tension and ch03's core radius.

**⚠️ η and δ in ch09 depend on the flow.** Blasius η = y/δ(x), δ = √(νx/U) (9.26) — identical to ch08's y/√(νt) with t = x/U
(factor 1), while ch08's figures plot y/(2√(νt)) (factor 2); Falkner–Skan δ = √(νx/U_e(x)) (9.34) (η ∝ x^{(n−1)/2}y); free jet
δ = (Cρν²x²/J)^{1/3} ∝ x^{2/3} with η ∈ (−∞, ∞) (9.59)–(9.64); wall jet δ ∝ x^{3/4} (9.82). The book's Fig. 9.7 axis is
½√(n + 1)η (E3 offers both). δ is the **similarity length, not δ₉₉** (δ₉₉ = 4.910δ for Blasius); δ̄ (§9.1) is only an
order of magnitude; δ* is the displacement thickness (the star also marks scaled variables in (9.6)–(9.8)).

**⚠️ The Falkner–Skan exponent is the book's n, the code's `m`.** `falkner_skan(m)`, `falkner_skan_state(m)`,
`falkner_skan_separation()["m_sep"]` = −0.090429; β = 2m/(m + 1) (Hartree; the literature's tables, Belden et al.) =
−0.198838 at separation. Thwaites' "m = −λ" (Fig. 9.8 axis) is not used in code. ch08's n and m are the similarity exponents.

**⚠️ Angles on bodies from the FORWARD stagnation point in ch09.** `phi_deg` (separation 82°, 125°, Thwaites 103.1°) is
measured from the upstream stagnation point, as the book does in §9.7–§9.9 and Exercise 9.21; ch06 measured θ from +x
(downstream: C_p = 1 − 4 sin²θ is the same function either way), ch08's sphere θ from the downstream axis. In D13 the outward
normal's x-component is n_x = −cos φ.

**⚠️ Sign of the wall curvature under an adverse gradient (R17).** (9.9) at the wall (u = v = 0) gives
μ(∂²u/∂y²)_wall = dp/dx: **positive** for an adverse gradient (dp/dx > 0, Falkner–Skan n < 0: f‴(0) = −n > 0), negative for a
favourable one. The book's sentence says "< 0 for n < 0" — a slip; `wall_curvature` follows (9.9). The (9.8) sign was also
coded negated once (review Must-fix 1): ∂p*/∂y* = −(1/Re)(u*v*_x* + v*v*_y*) + (1/Re²)v*_x*x* + (1/Re)v*_y*y*
(`bl_dpdy_scaled`), pinned by a manufactured-field sign test.

**⚠️ One side or two (R5, a trap, not a slip).** (9.33) C_D = 1.328/√Re_L and F_D are per unit width for **one** face, as the
book says; a plate wetted on both faces has twice the drag (`blasius_drag(..., sides=1)`, `plate_drag_coefficient(sides=)`).

**⚠️ Strouhal number with the cyclic frequency.** St = fd/U ≈ 0.2 with f in Hz (the book's Ω in (4.102) is the cyclic
frequency); with the angular 2πf the same street would give 1.26. `shedding_frequency(U, d, St)` returns both f and
`omega_rad`.

**⚠️ Two separation criteria in Thwaites' method.** `LAMBDA_SEP_BOOK` = −0.09 (Table 9.1: l(−0.09) = 0) and `LAMBDA_SEP_FS` =
−0.0681483 (computed lazily from the Falkner–Skan fold; l = 0 of our exact closure). `thwaites(closure="falkner_skan")` defaults
to the second, `closure="white"` to the first; `example_9_2()["criteria"]` lists the book's first (x/L = 0.158, then 0.126).
Always show both, the book's first.

**⚠️ Printed ch09 forms coded corrected (analysis §9), printed forms kept as options a test must fail.** (9.7) with ∂x*², ∂y*²
(`bl_nondim_sympy(printed_9_7=True)`) · η₉₉ = 4.910, not 4.93 (`blasius_delta99(printed=True)`) · wall-jet ODE
**4**f‴ + ff″ + 2f′² = 0 (`wall_jet_ode_solve(printed=True)` / `coeff=1.0`, `similarity_reduce_sympy("wall_jet", printed=True)`) ·
wall-jet separation of variables with f^{1/2} · h₉₉ coefficient 7.3319, not 5.6152 (the 4 % point;
`free_jet_halfwidth(printed=True)`) · (9.56) with the kinematic stress τ/ρ · Ψ of (9.85) in m⁵/s³ · Magnus "Re > Re_cr" for the
positive effect (`magnus_sign`) · no attached bounded Falkner–Skan solution below the fold (R11) · turbulent jets are Ch. 12,
not 13 · the wall-curvature sign (R17).

**⚠️ ch10: one letter, many meanings — the notebook's choices (conventions table, cell 6).** The book reuses letters freely
in §10.2–10.5; the notebook keeps the book's symbol where its equations are quoted and renames only where two meanings meet:

| Book letter | Meanings in ch10 | What we write / code name |
|---|---|---|
| i | grid index **and** √−1 in (10.20)–(10.24) | upright $\mathrm i$ (`1j`) for √−1; grid index j in our own derivation lines |
| α, β | FTCS numbers α = uΔt/(2Δx), β = DΔt/Δx² (10.11); Glowinski Θ-scheme weights (10.129)–(10.133); FE time weights (10.163) | FTCS `alpha`, `beta`; **Courant C = 2α** (`cfl_number`); α_Θ, β_Θ (`theta_split`); α_t, β_t (`alpha_t`, `beta_t`) |
| θ | Fourier angle θ = kπΔx (10.26); Θ-scheme fraction; Ch. 8's θ-method in `crank_nicolson_1d`/`FEM1.solve_transport(theta=)` | θ = Fourier angle (rad per cell, `theta`); **Θ** = Glowinski's fraction (1 − 1/√2); the θ-method weight is named where used |
| g | Dirichlet value T(0) = g (10.2); body force per unit mass **g** (10.79), (g_x, g_y) = the book's (f, g) (10.119)–(10.120) | `g=` Dirichlet value; `body=` force (g_x, g_y); Fourier amplitude ξ̂ⁿ (book gⁿ(k)) |
| D | diffusivity (10.1); strain-rate tensor D[u] (10.136); cavity/block side | `D` [m²/s]; **D**[u] bold with brackets (write 2**D**[u], never "2D"); sides L (cavity), d (block, cylinder) |
| M, K, S | mass matrix and Mach number; error constant (10.15) and stiffness matrix; trial space and Strouhal number | **M** and **Ma** (the book prints M in (10.151)–(10.155)); K_e and **K**; 𝒮 and St |
| n | time level; number of cells/elements; shedding frequency in S = nd/U | superscript n; n cells/elements (n_el in element loads); f_s |
| R | global Péclet uL/D (10.87); cell Péclet uΔx/D (10.31) | `R`, `R_cell` (`cell_peclet`) |
| δ | layer thickness (10.89); Kronecker delta | δ; Kronecker always with two indices δ_AB |
| F, G | FE load vector and y-flux of U_t + E_x + F_y = 0; amplification factor and weak-form residual row | the text says which **F**; G(θ) vs G_A |
| p, p′ | pressure; "order p"; p′ = Newton correction (10.164) (ch04/ch07: perturbation pressure) | "observed order p" always next to its error; p′ named "Newton correction" |
| c | sound speed (10.95), (10.99); coefficients c₁–c₅ (10.109); FE coefficients c_A (10.45) | `c` sound speed; `ns_coefficients()["c1"…]`; c_A in D11 only |
| λ | eigenvalue (y′ = λy test, split system, β_h² = λ_min) | ch09's Thwaites λ and street eigenvalue, ch07 wavelength — say which |
| R1…R12 | the analysis's printed slips | **"slip #1…#12" in the notebook** so they never read as recap IDs R01…R13 |

**⚠️ ch10 Fourier convention.** The book writes a mode as e^{iπkx_i} (10.20), so k counts *half*-waves per unit length and
θ = kπΔx (10.26); the zigzag 1, −1, 1, … is θ = π. Every `FD` function takes θ directly; `fourier_mode(x, k,
convention="book" | "standard")` gives both (standard e^{ikx}, θ = kΔx).

**⚠️ ch10 upwind side and CFL follow the sign of u (slip R11).** (10.29) T_i − T_{i−1} and (10.30) uΔt/Δx ≤ 1 assume u > 0;
the code takes the upwind side from sign(u) and the condition |u|Δt/Δx ≤ 1 (`scheme="upwind_printed"` keeps the book's
stencil and explodes for u < 0). `cfl_number(u, dt, dx)` is **signed**. With diffusion: upwind |C| + 2β ≤ 1, Lax–Wendroff /
MacCormack C² + 2β ≤ 1. The book's "forward difference" after (10.93) is the **backward** (upwind) difference (slip R3);
`steady_cd_fd(scheme="forward")` is the true downwind one (root 1/(1 − R_cell): wiggles only for R_cell > 1).

**⚠️ ch10 truncation-error sign.** E sits on the left of (10.16): E = (exact solution put into the scheme's difference
operator) − (the PDE) = "exact minus scheme" divided by Δt; (10.17) E = +(Δt/2)T_tt + u(Δx²/6)T_xxx − D(Δx²/12)T_xxxx. The
forward stencil's leading error is **−**(Δx/2)T_xx (D01). Upwind's numerical diffusivity is +uΔx/2 (steady), +|u|Δx(1 − C)/2
(unsteady; negative — anti-diffusion — for C > 1).

**⚠️ ch10 three grids.** Node-based FD grid x_i = iΔx (FD; MCK with ρ, ρu, ρv collocated at nodes, walls on node lines,
periodic grids do not repeat the end node); element meshes (FEM1: nodes x_0 … x_n, element e = [x_{e−1}, x_e], local a = 1, 2
↔ global A = e − 1, e; FEM2: vertices first then edge mid-nodes, local 1–3 vertices, 4 = mid(1,2), 5 = mid(2,3), 6 = mid(3,1),
unknown vector [u (N), v (N), p (V)]); staggered cell grid (MAC: `p[j, i]` at ((i + ½)Δx, (j + ½)Δy), `u[j, i]` at x = iΔx,
`v[j, i]` at y = jΔy; the book's u_{i+1/2, j} ↔ `u[j, i+1]`, i.e. `u[j, i]` = u_{i−1/2, j}). Wall-face normal velocities are
data, never corrected; tangential no-slip by ghost values u_g = 2U − u₀ (linear, default) or quadratic (DEVIATION: the book
leaves it open).

**⚠️ ch10 weak-continuity sign.** (10.159)/(10.162) carry a minus sign chosen so the pressure blocks are B and Bᵀ (a symmetric
saddle point) — do not "fix" it. Slip R6: the second sum of (10.172) must be v_{A′} (`assemble_newton_system(printed_10_172=
True)` gives a 1.7e55 Poiseuille error); R7: (10.186) prints v′ for the pressure expansion p′; R10: (10.166)'s star on
β_t ∂v/∂t(t_n) is spurious (known data).

**⚠️ ch10 non-dimensionalisations change inside the chapter.** §10.2–10.3 dimensional (u [m/s], D [m²/s]); §10.4 NS
non-dimensional with Re (10.81) (lengths L, speed U, time L/U, pressure ρU²); MacCormack NS (10.103)–(10.109) dimensional (μ, c);
the cavity algorithm non-dimensional with p = ρ/Ma², Re = ρ₀UL/μ; block lengths in block sides; FE cylinder in diameters d
(t̄ = tU/d); force coefficients per span on ½ρU²d, torque on ½ρU²d². St with the **cyclic** frequency (as ch09).

**⚠️ ch10 fitted vortex minimum.** ψ_min of the cavity is negative (clockwise primary eddy); `primary_vortex_centre` fits a
parabola through the discrete minimum and its neighbours — the vertex value is f₁ − (f₂ − f₀)²/(8(f₂ − 2f₁ + f₀)) per direction (x and y corrections add),
**always ≤ the grid minimum** (the code once added the correction with the wrong sign; review M1).

**⚠️ ch11: σ is the growth rate; two eigenvalue conventions.** §11.4, §11.6, §11.14 write a mode as e^{σt} with
σ = σ_r + iσ_i (σ_r > 0 grows; σ_i ≠ 0 oscillates or travels); §11.3 and §11.7–11.11 write e^{ik(x − ct)} with c = c_r + ic_i.
The two are the same thing: **σ = −ikc**, so σ_r = kc_i (growth rate) and σ_i = −kc_r (a wave moving to +x has σ_i < 0)
(`ST.sigma_from_c`, `c_from_sigma`). In ch11 code `sigma` is always a growth rate; surface tension is `surface_tension`
(σ_s in text; ch01/ch07 called it σ), the ch03 vortex core radius is still `sigma` there, ch10's MacCormack safety factor σ
is unrelated. "Stable" means σ_r ≤ 0 **for every k**.

**⚠️ ch11: Γ has three temperature-gradient conventions and is also a circulation.**

| Convention | Definition | Heated from below | Where |
|---|---|---|---|
| Kundu Ch. 1 (the code's lapse rate) | Γ ≡ dT/dz | negative | `core.stratification`, every earlier chapter |
| meteorology (shown alongside) | Γ ≡ −dT/dz | positive | notebook tables, slider legends |
| Bénard, (11.21) and (11.24) (slip S10) | Γ = −dT̄/dz = ΔT/d | positive | §11.4 only; **no ch11 function takes a Γ** — `rayleigh_number` takes `dT` = T_bottom − T_top, `gamma_conventions(dT, d)` returns all three numbers |

In §11.6 Γ = 2πrU_θ is the circulation of a fluid ring (`ring_interchange_energy`, `rayleigh_circulation_criterion`).

**⚠️ ch11: the Rayleigh number changes sign between §11.4 and §11.5.** (11.21) Ra = gαΓd⁴/(κν) > 0 when heated from below.
§11.5 defines Ra ≡ gαd⁴(dT̄/dz)/(νκ) with the gradient itself, so it is **negative** when heated from below and (11.45)
carries −Ra (`thermal_rayleigh_signed`); Rs = gβd⁴(dS/dz)/(νκ_s) uses κ_s, Rs′ uses κ. The criterion (11.46) reads
Rs − Ra = 27π⁴/4 at the marginal state.

**⚠️ ch11: three non-dimensionalisations and a Reynolds number per flow.** Bénard: lengths by d, time by d²/κ (σ in κ/d²),
w left dimensional until W ≡ (Γd²/κ)ŵ. Taylor: gap d, x = (R − R₁)/d ∈ [0, 1], σ in ν/d². Parallel flows: L and U₀ **per
flow** — plane Poiseuille half-width and centreline speed (Re_c 5772), Blasius δ* and U∞ (Re_c 519), tanh layer and Bickley
jet their L and U₀, pipe U_max d. Never compare two critical Reynolds numbers without naming the length.

**⚠️ ch11: three stream-function sign conventions (slip S9).** §11.7 u = ∂ψ/∂z, w = −∂ψ/∂x (as ch07); §11.8 u = ∂ψ/∂y,
v = −∂ψ/∂x (as ch04), so û = φ′, v̂ = −ikφ; §11.14 and our Bénard rolls u = −∂ψ/∂z, w = ∂ψ/∂x. Each function that returns
velocities states its convention.

**⚠️ ch11: one letter, many meanings.**

| Letter | Meanings in ch11 | What we write / code name |
|---|---|---|
| K, k, κ | K horizontal wavenumber magnitude of a Bénard cell (non-dimensional); k streamwise or axial wavenumber; κ thermal diffusivity, κ_s salt diffusivity | `K`, `k`, `kappa`, `kappa_s` |
| α | thermal expansion coefficient (§11.4–11.5); the book's α = Ω₂/Ω₁ − 1 (§11.6); Orszag's α = streamwise wavenumber | `alpha` (thermal); `mu` = Ω₂/Ω₁ (**not viscosity**); `k` |
| β | haline contraction coefficient (not ch10's FTCS β, not the planetary β of Ch. 13) | `beta_S` |
| μ | speed ratio Ω₂/Ω₁ of the cylinders (dynamic viscosity never appears in ch11 code; ν does) | `mu` |
| R | radius (§11.6, capital, as ch08); shear/density thickness ratio of the tanh layer, N² = J sech²(Rz) (`tg_growth(k, J, R)`); the book's nonlinearity parameter in §11.14 | `R`, `R1`, `R2`; `R` |
| r | Lorenz's reduced Rayleigh number Ra/Ra_c(k); radius of a fluid ring | `r`; `r1`, `r2` |
| J | bulk Richardson number of the tanh/sech² layer (Ri at the centre) — **not** ch09's jet momentum flux, **not** ch10's FE Jacobian (the Lorenz Jacobian matrix is `lorenz_jacobian`) | `J` |
| φ | velocity potential (§11.3); ψ̂/(U − c)^{1/2} (§11.7, Miles–Howard); Orr–Sommerfeld / Rayleigh amplitude (§11.8–11.9) | `phi` (the docstring says which) |
| U₁ | upper-stream speed (§11.3); speed at the inflection point (§11.9, Fjørtoft) | `U1`; `U_I` |
| Pr | Prandtl number ν/κ; in the Lorenz system it plays the role the literature calls σ | `Pr` |
| b | Lorenz's geometric factor 4π²/(π² + k²) = 8/3; half-width of the sin y channel (`sin_profile_max_growth(b)`) | `b` |
| N | buoyancy frequency (N², may be negative); polynomial degree of the Chebyshev grid | `N2`; `N` |
| S | salinity (§11.5; ch10's 𝒮 trial space, ch02's strain rate elsewhere) | `S`, `dSdz` |

**⚠️ ch11: cross-reference slip S13.** The book's p. 503 prints "(7.128)" beside $N^2\equiv-\frac{g}{\rho_0}\frac{d\bar\rho}{dz}$;
Chapter 7 defines the buoyancy frequency as **(7.127)** — cite (7.127). Printed slips are S1–S13 in code and analysis and
"slip #1…#13" in the notebook (SKIP ids S01, S02 and recap ids R01–R24 are different things).

**⚠️ ch11: what a zero or a NaN means in a table.** `tg_growth` / `tg_growth_map`: a 0 within 0.006 of the neutral curve
J = k(1 − k) means "kc_i < 0.004", not "stable". `rayleigh_spectrum_table`: c_i = 0 means c_i ≤ 1e-4. `bickley_neutral_curve`:
NaN in `k_lower` means "unstable down to k = 0.02 (edge not resolved)", never "stable"; do not join `k_lower` across NaN.
`gradient_richardson` returns ±∞ at U′ = 0 with the sign of N² (nan for 0/0).

**⚠️ ch12: lower-case letters are fluctuations.** In this chapter $\tilde u_i=U_i+u_i$ (12.24): tilde = total field,
capital or over-bar = mean, **lower-case u, v, w and T′ = fluctuation** (so `u` is not "the velocity"). The over-bar is an
**ensemble** average (12.1) unless a function name says `time_` (12.2) or `volume_` (12.3). Ensembles carry the member index on
axis 0 (`samples[n, ...]`).

**⚠️ ch12: lapse rate and "temperature" — the rule with consequences for Ch. 13.** Compute with Kundu's Γ ≡ dT/dz
(Γ_a = −g/C_p ≈ −9.8 K/km; stable when dT/dz > Γ_a) and always show the meteorological Γ_met ≡ −dT/dz (Γ_d ≈ +9.8 K/km;
stable when Γ_met < Γ_d) alongside. T̄ and T′ in every buoyancy term of the chapter — (12.30), (12.47), (12.106),
$\mathrm{Ri}=N^2/(dU/dz)^2$ (12.108), (12.112) — are **potential** temperature; with a thermometer gradient
N² = gα(dT/dz − Γ_a) = gα(Γ_d − Γ_met). **`gradient_richardson_thermal(dTdz, dUdz, alpha, g, *, Gamma_a, convention=)`
requires `Gamma_a`**: its former default 0.0 called the standard atmosphere unstable (Ri = −2.213 instead of +1.110). A
θ-gradient is passed with `Gamma_a=0.0` on purpose; parity rows pass `ch12.adiabatic_lapse_rate()`. **The verdict string
of ch01's `lapse_rate_stability` writes "Γ < Γa" in its meteorological line; read it as Γ_met < Γ_d** (its "Γa" is the
positive +9.8 K/km) — explained wherever the string is shown; changing the string across chapters is an open decision.
The verdict strings start with the criterion ("stable ⇔ …") even for an unstable layer: take the word from `["verdict"]`.

**⚠️ ch12: signs of the fluxes.** Heat flux positive **upward**: `wT` = $\overline{wT'}$ [K m/s], `H` = ρC_p·`wT` [W/m²].
Upward flux ⇒ $\mathrm{Rf}=\frac{-g\alpha\overline{wT'}}{-\overline{uw}(dU/dz)}$ (12.107) < 0 and
$L_M\equiv-u_*^3/(\kappa\alpha g\overline{wT'})$ (12.110) < 0 (unstable); L_M > 0 stable; `inf` neutral. The Reynolds stress is
−ρ₀$\overline{u_iu_j}$; **dictionary keys `uv_plus` hold MINUS the correlation, −$\overline{uv}$/u_*² ≥ 0 (alias
`minus_uv_plus`)**, while `mean_energy_budget` and `tke_budget` take the true $\overline{uv}$ (< 0 where dU/dy > 0).
Slip #17: the page writes "$\overline{uw}=u_*^2$" where −$\overline{uw}$ = u_*² is meant.

**⚠️ ch12: spectra are two-sided.** $S_e(\omega)=\frac1{2\pi}\int R_{11}e^{-i\omega\tau}d\tau$ (12.20) with the 1/2π in the
forward transform, angular frequency [rad/s] or wavenumber [rad/m], and $\int_{-\infty}^{\infty}S\,d\omega=$ variance
(12.22). Functions return S on ω ≥ 0 only (variance = 2∫₀^∞S dω, `spectrum_variance`); `one_sided=True` doubles the density.
`inertial_spectrum_1d(C1=…)` always takes the **one-sided** constant and the default `two_sided=True` halves it (one-sided
C₁ = 0.491, two-sided 0.2455 for C = 1.5). E(K) is defined on K ≥ 0 with ∫E dK = ē (the book's letter is S(K)).
$\overline{u^2}$ in §12.6 is **one** component (ē = (3/2)$\overline{u^2}$). The book's skewness and kurtosis are
un-normalised central moments (`statistics(normalized=True)` gives the usual ones).

**⚠️ ch12: what a classifier does not say.** `surface_layer_regime(z, L_M)` returns "forced convection" for any height
well below |L_M|, **stable or unstable** — it says nothing about stability (use `surface_layer_state(...)["verdict"]`).
`dispersion_regime` has inclusive boundaries on the ratio t/Λ_t (≤ 0.3 ballistic, ≥ 3 diffusive). The log-linear wind
$U=\frac{u_*}\kappa[\ln\frac z{z_0}+5\frac z{L_M}]$ returns NaN where 1 + βz/L_M ≤ 0: in unstable air it **ends at
z = |L_M|/5**. `WT.LOG_LAW_CONSTANTS["channel_dns_lee_moser_2015"]["B"]` is `None` (the source quotes no additive
constant). A fitted (κ, B) belongs to its window (`fit_log_law(window=(y⁺_min, y/δ_max))`).

**⚠️ ch12: one letter, many meanings** (the notebook's conventions cell shows this table; `ch12.conventions()`).

| Letter | Meanings in ch12 (and earlier) | What we write / code name |
|---|---|---|
| κ | von Kármán constant in $U^+=\frac1\kappa\ln(y^+)+B$ (12.88) and in (12.110); thermal diffusivity in (12.31) and (12.112); κ_m, κ_T nearby | κ = von Kármán `kappa` (0.41 in §12.9–12.10, 0.4 in the surface-layer block); κ_th molecular `kappa_th`; κ_T eddy `kappa_T` |
| k, K | thermal conductivity in (12.32); wavenumbers k₁ and K; the "k" of k–ε (the book's ē); kurtosis K | `k_th`; `k1`, `K`; `e` ("k" only in the name k–ε) |
| e, E, ε | ē = ½$\overline{u_i^2}$ (ch01: internal energy); Ē = ½U_i² (ch11: disturbance energy); ε̄ dissipation, ε̄_T its thermal twin; E(K) 3-D spectrum | `e`, `E_mean`, `eps`, `eps_T` |
| λ, Λ | Taylor microscales λ_t, λ_f, λ_g, λ_T (earlier: wavelength, eigenvalue, Thwaites parameter, Lyapunov exponent); integral scales Λ_t, Λ_f, Λ_g (ch11: Λ = dissipation) | always with the subscript; `lambda_t`, `lambda_f`, `lambda_g`; `Lambda_t`, `Lambda_f`, `Lambda_g` |
| η | Kolmogorov length (ch07 surface elevation; ch08–ch09 similarity variable); η_T Batchelor scale | `eta`, `eta_T` |
| f, g, F, G | correlation coefficients f(r), g(r) (12.38); wall function in (12.80); Darcy f; gravity g; tensor functions in (12.40); jet profiles in (12.56)–(12.57); defect function in (12.84) | f(r), g(r) `f`, `g_corr`; f_w(y⁺); f_D `fD`; F_R, G_R; F, G jet; F_d defect |
| R_ij, r | correlation tensor (ch02–ch03: rotation tensor); r separation or radius; r₁₁, r_α correlation coefficients | `R`, `r` |
| τ | time lag in (12.17); decay time (Ex. 12.1); stress τ̄, τ₀; τ_c memory time of our test signals; τ_η Kolmogorov time | `lag`, `tau`, `tau0`, `tau_c` |
| δ, θ, α | δ jet width, layer thickness, δ_ij; θ momentum thickness (§12.9) vs potential temperature (§12.11, ours); α thermal expansion vs the no-sum index of §12.12 | `delta`; `theta_m` vs `theta`; `alpha` |
| S, N, L, σ, Π, U₀ | spectra S_e, S₁₁, S_T vs strain rates S̄_ij, S′_ij vs skewness; N members vs buoyancy frequency; L outer scale, walk step, L_M; σ Gaussian width with **σ² = 2νt** (ch03 core radius: 4νt) and σ_e, σ_ε; Π wake strength (ch01: Π groups); U₀ probe speed vs nozzle speed | `n_members`, `N2`; `L`, `step`, `L_M`; `Pi` |
| h, d | **h = full channel height** in $dP/dx=-2\tau_0/h$ (12.90) (ch08 used half-heights in places); d slot width, nozzle or pipe diameter | `h`, `delta` = h/2, `d` |
| C₁…C₉ | jet constants (slip #4: the text's "C₃ and C₄" are C₃ and C₅); C₁ also the one-dimensional Kolmogorov constant; C_μ, C_ε1, C_ε2 | `C3`, `C5`, `C6`; `C1`; `K_EPSILON_CONSTANTS` |

**⚠️ ch12: 18 printed slips** ("slip #1…#18", `ch12.book_slips()`): #1 the exponent of (12.54) is −5/3; #6 ∂P/∂x_i in
(12.97); #12 the condition of (12.129) is $t\gg\Lambda_t$; #13 plume width ∝ x near the source, ∝ √x far; #15 the ½ in the
molecular transport of (12.112); #16 the vanishing flux term is ρ₀U$\bar v$; #17 −$\overline{uw}$ = u_*²; #18 V does not
vanish at the jet edge. #5 is worded "the printed exponential family fails; exponentials as such are not excluded". The
one-sided Schwartz inequality (12.16) is a trap, not a slip.

**⚠️ ch13: the sign of f, and every direction word.** $f=2\Omega\sin\theta$ (13.8) is positive in the northern hemisphere
and negative in the southern. The book writes its formulas for f > 0 ("to the right", "clockwise", "coast on the right",
√f, $e^{-fy/c}$). Every `core.gfd` function accepts either sign: abs(f) for scales, sign(f) for directions, and a
`ValueError` at f = 0 where f is in a denominator (`inertial_period(0)` is the one documented `inf`). `*_printed` forms
refuse f ≤ 0. Latitudes are in radians (`lat_rad`); degrees only in `*_deg` names. The rotation rate is the **sidereal**
`OMEGA_EARTH`; the book's one turn per solar day (`OMEGA_SOLAR_DAY`) is 0.27 % smaller (trap T1). **The term
−2Ωu cos θ of $2\boldsymbol\Omega\times\mathbf u\cong(-fv,\ fu,\ -2\Omega u\cos\theta)$ (13.7) is an acceleration term on the
left of the vertical equation; the Coriolis force per unit mass is +2Ωu cos θ** (trap T18; ch04's convention callout).

**⚠️ ch13: where z = 0 is changes from section to section** (trap T4). Sea surface with the ocean in z < 0: §13.6 (surface
Ekman layer), §13.9 (vertical modes, nodes from z[0] = −H to z[-1] = 0), §13.14. Solid surface with the fluid in z > 0:
§13.7 (bottom Ekman layer). Flat bottom: §13.8 (shallow water). Lower lid: §13.17 (Eady). Each function's docstring and
each notebook block says which.

**⚠️ ch13: lapse rate (the standing rule) and "density".** Compute with Kundu's Γ ≡ dT/dz (Γ_a = −g/C_p ≈ −9.8 K/km; stable
when dT/dz > Γ_a) and always show the meteorological Γ_met ≡ −dT/dz (Γ_d ≈ +9.8 K/km; stable when Γ_met < Γ_d) beside it:
`ch13.lapse_rate_table(dT_dz, *, Gamma_a)` returns both rows with the same verdict; `Gamma_a` comes from
`STRAT.adiabatic_lapse_rate()` and is passed by keyword. The section itself quotes rates of decrease in words
(meteorological magnitudes). N² in this chapter is always from **potential** density (or potential temperature). From the
thin-shell equations on, **p and ρ are perturbations with the primes dropped** (trap T3). "Close to neutral" in §13.2
compares with the moist adiabat; the standard troposphere is stable to dry displacements (trap T2).

**⚠️ ch13: stream function, Rossby number, Rossby radius, spectrum — each differs from an earlier chapter.**
Stream function: **u = −∂ψ/∂y, v = ∂ψ/∂x** (§13.5, §13.16; `geostrophic_streamfunction`, `barotropic_velocity`) — the
opposite sign to ch04, ch06 and ch11 (u = ∂ψ/∂y); the eigenvalue c of the stability problem is unchanged (trap T14).
Rossby number: $\mathrm{Ro}=\frac U{fL}$ (13.13) with abs(f) (`GFD.rossby_number(U, f, L)`); ch04's `SIM.rossby_number` is
U/(2ΩL) — they differ by sin θ (trap T17). Rossby radius Λ: external c/abs(f); internal NH/(nπ abs(f)) or √(g′H₁)/abs(f);
**the Eady radius Λ_E = NH/abs(f) has no π** (trap T11; `rossby_radius_internal(with_pi=)`). Spectrum of §13.18:
**one-sided in K with mean(u²) = ∫₀^∞ S dK and no factor ½** (ch12: two-sided, ∫S = variance; trap T16). Ekman thickness
$\delta=\sqrt{2\nu_v/f}$ (13.29) is an e-folding scale; the oceanographic "Ekman depth" is πδ
(`ekman_depth(convention="efold"|"pi")`, trap T6).

**⚠️ ch13: signs returned by the wave functions.** `rossby_omega` returns **signed ω for signed k** (a wave with ω > 0 has
k < 0; "maximum phase speed" is a maximum of magnitude; trap T8). `shallow_water_omega(k, l, c, f0, beta)` returns the
three roots of $\omega^3-c^2\omega K^2-f_0^2\omega-c^2\beta k=0$ (13.76): two fast roots of **opposite sign** and the slow
root, or **NaN where the discriminant is negative** (very long waves: "all roots real" is a β-plane scaling statement,
trap T9). `kelvin_decay_side(f, direction)` says which direction of travel is trapped. Modes: `Modes.c` is in decreasing
order; **with a rigid lid index 0 is the first baroclinic mode** (the book's n = 1); ψ_n(0) = 1; the modal amplitudes have
different units — u_n, v_n [m/s], p_n [m²/s²], w_n [1/s], ρ_n [kg/m²] (trap T12). A wind is named by where it comes from,
a current by where it goes (`wind_from_to`). The "orbit" figure of a Poincaré wave is a velocity hodograph (trap T7).

**⚠️ ch13: argument order of sibling functions.** `eady_max_growth_rate(f, N, dUdz)` and `eady_time_scale(f, N, dUdz)`
take (f, N) while `eady_alpha(k, l, N, f)`, `eady_growth_rate(k, l, N, f, H, U0)` take (N, f); `poincare_omega(K, f, c)`
against `shallow_water_omega(k, l, c, f0, beta)` and `rossby_omega(k, l, beta, f0, c, U)`. **Call them by keyword.**

**⚠️ ch13: one letter, many meanings** (the notebook's conventions block shows this table; `ch13.conventions_table()`).

| Letter | Meanings in ch13 (and earlier) | What we write / code name |
|---|---|---|
| f | Coriolis parameter, signed (earlier: a generic function, a frequency in Hz, the Blasius function, ch12's correlation function) | f, f₀ at the central latitude; a frequency is always ω; `f`, `f0` |
| Ω, ω | Earth's rotation rate; wave frequency; ω_x, ω_y horizontal vorticity in the Ekman layer (13.31) | `Omega`, `omega`, `omega_x` |
| β | df/dy = 2Ω cos θ₀/R in $f=f_0+\beta y$ (13.10); once the haline contraction coefficient (§13.3); ch10 a scheme parameter, ch12 the log-linear coefficient | β = df/dy only; β_S for salt; `beta` |
| N | buoyancy frequency from potential density (ch12: a count of members) | N; counts are n, `nx`, `n_modes`; `N`, `N2` |
| H, h, η | layer depth, ocean depth, WKB scale of N, lid separation, scale height c²/g; h total depth over an uneven bottom; η surface displacement (ch12: Kolmogorov length) | H; H_N; H_s; H_e equivalent depth; h; η; `H`, `He`, `h`, `eta` |
| c | sound speed (§13.2–13.3); long-wave speed √(gH) (§13.8 on); modal speed c_n; complex phase speed c = c_r + ic_i (§13.16–13.17); c_x, **c**_g | c_s; c; c_n; "complex c" said each time; `c`, `c_n` |
| k, l, m, K | eastward, northward, vertical wavenumbers; K horizontal magnitude (the 3-D wavevector in §13.14); K the perturbation kinetic energy (§13.17); K₀, K₁, K₂ of the triad; K(z) our eddy-viscosity profile | k, l, m, K; KE; K_v(z); `k`, `l`, `m`, `K` |
| ζ | relative vorticity ∂v/∂x − ∂u/∂y (ch07: interface displacement; ch12: z/L_M) | ζ, ζ_g, ζ̄; `zeta` |
| θ | latitude (§13.4); angle of the wavevector with the horizontal (§13.14); potential temperature | θ latitude `lat_rad`; θ_K; θ_p |
| δ | Ekman e-folding thickness (earlier: boundary-layer thickness, Kronecker delta) | δ; D_E = πδ; `ekman_depth(convention=)` |
| Λ, λ | Rossby radius (three kinds); ch12: integral scales; λ wavelength (ch12: Taylor microscale) | Λ with the wave speed named; Λ_E; `rossby_radius*`; `lam` |
| ψ | vertical mode ψ_n(z) (§13.9); stream function with u = −∂ψ/∂y (§13.5, §13.16) | ψ_n; ψ with its sign stated; `psi` |
| E, Ro, R | Ekman number ν/(abs(f)L²) (earlier: energy); Rossby number; R Earth's radius (earlier: gas constant, correlation) | `ekman_number`, `rossby_number`, `EARTH_RADIUS_MEAN` |
| α | thermal expansion (§13.3 and our temperature form); Eady wavenumber α = NK/abs(f) (13.139); enstrophy flux (§13.18) | α_T (`alpha`, keyword-only); α (`alphaH`); α_Z (`alpha_ens`) |
| U, V, τ | velocity scale, geostrophic interior velocity, mean current, U₀z/H; **V = u + iv complex velocity** (§13.6–13.7); τ wind stress | U named each time; `U_g`, `V_g`; `tau_x`, `tau_y` |
| S(K), i | energy spectrum, one-sided, no ½ (ch12: two-sided); i = √−1 and the label in T_i (slip #5) | `enstrophy_spectrum`, `barotropic_spectrum` state their normalisation |

**⚠️ ch13: 14 printed slips and 18 traps** ("slip #1…#14", `ch13.book_slips()`, each row with `kind`): #1 the sign of ∇p
in (13.2); #2 the friction force per unit mass needs 1/ρ; #3 the sign of (13.17) is $2\Omega u=-\frac1\rho\frac{\partial p}{\partial y}$;
#4 the decay condition (13.26) is for z → −∞; #6 the first root of $\tan\frac{NH}{c_n}=\frac{c_nN}g$ (13.69) is ≪ 1; #7 the
ω³ term of the cubic is negligible for ω ≪ f; #8 N depends on depth in §13.14; **#14 westward flow over a step keeps a
permanent shift** (true as printed only for a ridge of finite width; the user has not yet ruled on it); #10 and #11 are
`kind = "loose"` (inconsistent numbers, not false statements); #9 four wrong cross-references; #13 a spelling (the
attribution is not asserted). `ekman_surface_printed` is a **trap**, not a slip: true as printed for f > 0.

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
| η ⚠️ (Stokes) | similarity variable y/(2√(νt)) of u = U erfc η | – | **η also = surface function (4.90)**; **ch08's (8.25) η = y/√(νt) has no factor 2** (see the ch08 η trap) | ch04 | `stokes_first_problem(y, t, U, nu)` (now in `core.laminar`, re-exported by ch04) |
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
| **ch08 — laminar flow** | | | | | |
| y, h ⚠️ (channel) | distance from the fixed wall; gap between the plates | m | walls y = 0 (fixed) and y = h (moving at U) — not ±h | ch01 → | `y`, `h` in `channel_flow` |
| U ⚠️ | wall speed (Couette, slider, Stokes' problems); stream speed (sphere); half the jump of the vortex sheet (Ex. 8.5) | m/s | u = ±U far from the sheet; **U₀, U_h = lower and upper gap walls** | ch01 → | `U`, `U_0`, `U_h` |
| dp/dx, dp/dz ⚠️ | streamwise pressure gradient (channel; pipe along z) | Pa/m | **book sign: favourable < 0**; `G=` alias = −dp/dx (ch04) | ch04 → | `dpdx`, `dpdz`, `G` |
| Q, V ⚠️ | flow rate (per unit width in a channel [m²/s]; volume flow in a pipe [m³/s]); mean velocity Q/h or Q/πa² | m²/s or m³/s; m/s | channel V = Q/h (printed middle form lacks 1/h); pipe u_max = 2V | ch08 | `channel_flow_rate` → (Q, V); `pipe_flow_rate` → (Q, V, u_max) |
| y* | height of the velocity extremum in Couette–Poiseuille flow, h/2 − μU/(h dp/dx) | m | interior only if 0 < y* < h; max (dp/dx < 0, toward the moving wall), min (dp/dx > 0, toward the fixed wall) | ch08 | `couette_poiseuille_state(...)["y_umax"]`, `["y_umin"]` |
| τ ⚠️, τ₀ | signed shear stress μ du/dy (channel), μ ∂u_z/∂R (pipe); pipe wall stress (a/2)dp/dz | Pa | τ₀ < 0 for forward flow (traction sign); the book draws \|τ\| | ch01 → | `channel_shear_stress`, `pipe_shear_stress`, `pipe_wall_stress` |
| f ⚠️ (friction) | Darcy friction factor 8τ₀/(ρV²) = 64/Re (laminar) | – | Re on the diameter and mean speed; **f = Coriolis parameter, traction, level set elsewhere** | ch08 → Ch. 12 | `pipe_friction_factor(Re)` |
| a ⚠️ | pipe radius (§8.2); sphere radius (§8.6) | m | Re_a = ρUa/μ uses it; **d = 2a** in the pipe Re | ch03 → | `a` |
| R₁, R₂, Ω₁, Ω₂ | inner/outer cylinder radii and rotation rates (circular Couette) | m, rad/s | R₂ = ∞ and R₁ = 0 are exact branches; Rayleigh-stable iff Ω₂/Ω₁ > (R₁/R₂)² (co-rotating) | ch08 → Ch. 11 | `circular_couette(R, R1, R2, Omega1, Omega2)` |
| A, B ⚠️ (Couette) | u_φ = AR + B/R constants: A = (Ω₂R₂² − Ω₁R₁²)/(R₂² − R₁²), B = −(Ω₂ − Ω₁)R₁²R₂²/(R₂² − R₁²) | 1/s, m²/s | ω_z = 2A; σ_Rφ = −2μB/R²; **A, B, C, D reused as integration constants in every ch08 derivation (D also = drag)** | ch08 → Ch. 11 | `circular_couette(..., return_coeffs=True)` |
| σ_Rφ ⚠️ | shear stress of circular Couette flow μR d(u_φ/R)/dR | Pa | on the fluid at R₁ the face normal is −e_R: power in = −2πR₁σ_Rφu_φ > 0 | ch08 | `circular_couette_shear_stress`, `circular_couette_power` |
| ε ⚠️ (fineness) | gap-to-length ratio h/L | – | ε ≪ 1; **ε = dissipation (ch04), ε² density ratio (ch07), kernel smoothing (ch05), roughness (ch01)** | ch08 → Ch. 9 | `lubrication_scales(...)["eps"]` |
| Re_L, ε²Re_L | passage Reynolds number ρUL/μ; reduced Reynolds number (weight of inertia in a thin gap) | – | lubrication needs ε²Re_L ≪ 1, not Re_L ≪ 1 (engine film 4350 vs 4.35e-3) | ch08 → Ch. 9 | `lubrication_scales(...)["Re_L"]`, `["eps2_Re_L"]` |
| Λ | bearing number μUL/(P_a h²) = viscous/pressure force ratio with p* = p/P_a | – | 49.35 for the engine film; Λ = 1 ⇔ the viscous scale μUL/h² | ch08 | `lubrication_scales(...)["Lambda"]` |
| P_a | absolute atmospheric pressure used as the lubrication pressure scale | Pa | 101 325 (`P_ATM`) | ch01 → | `p_a=P_ATM` |
| h(x, t) ⚠️, q | local gap (or film thickness); gap flux ∫₀^h u dy = −(h³/12μ)p_x + (U₀ + U_h)h/2 | m; m²/s | Reynolds equation h_t + q_x = 0; thin film q = −(ρg/3μ)h³h_x | ch08 | `lubrication_flux`, `thin_film_flux`, `reynolds_pressure_1d` |
| h₀, α ⚠️ (slider), L, C₁, W | slider inlet gap, taper h = h₀(1 + αx/L), pad length, pad-frame flux constant, load per unit width | m, –, m, m²/s, N/m | α > 0 widens toward x = L; W < 0 if αU < 0; optimum 1 + α = 2.18870; **α = thermal expansion, ray angle elsewhere; C₁ = Ex. 8.4's proportionality constant** | ch08 | `slider_bearing(x, h0, alpha, L, U, mu, p_e, model)`, `slider_bearing_load`, `slider_bearing_state` |
| h_m | gap at the slider's pressure peak, 2(1 + α)h₀/(2 + α) | m | recirculation where h > 1.5h_m ⇔ gap ratio > 2 | ch08 | (inside `slider_bearing_state`) |
| p_e | ambient (end) pressure of the slider | Pa | p(0) = p(L) = p_e | ch08 | `p_e` |
| φ ⚠️ (Hele-Shaw) | gap-averaged velocity potential, ū = ∇φ = −(h²/12μ)∇p | m²/s | the gap coordinate is **z** in Example 8.2; **φ = cylindrical azimuth elsewhere in ch08** | ch08 | `hele_shaw_potential`, `hele_shaw_mean_velocity` |
| β ⚠️, x_N, η_N | thin-film coefficient ρg/(3μ) (Huppert's β includes Δρ); front position; front constant 1.411245 | 1/(m s); m; – | x_N = η_N(βA³)^{1/5}t^{1/5}; A = area (volume per width) | ch08 → Ch. 13 | `viscous_current_similarity(x, t, area, …)`, `viscous_current_eta_N()`, `thin_film_state` |
| h_min | precursor film thickness in the spreading solver | m | keeps the diffusivity ρgh³/3μ > 0 ahead of the front; 1e-6 default | ch08 | `thin_film_spread(..., h_min=1e-6)` |
| η ⚠️ (similarity) | y/√(νt) (8.25) | – | **ch04 and the book's figures: y/(2√(νt))**; ch07: surface elevation; ch04: level set | ch08 → Ch. 9 | `similarity_variable(y, t, nu, half=False)` |
| F(η) | similarity profile u/U = F(η) | – | Stokes' first problem F = erfc(η/2) | ch08 | `similarity_ode_solve(case)` |
| δ₉₉ ⚠️, δ(t) ⚠️ | 99 % diffusion thickness 2 erfc⁻¹(0.01)√(νt) = 3.643√(νt); the time-dependent length scale of the ansatz (8.32) | m | **δ also Stokes-layer depth (§8.5), interface thickness (ch05)** | ch08 | `diffusion_thickness(t, nu, level=0.01)`, `transition_width` |
| γ ⚠️ (ansatz), n, m ⚠️, A, D | dependent variable of γ = At^{−n}F(ξ/δ(t)); decay and spreading exponents with δ = Dt^m | varies | matched powers of t + one conserved quantity fix n, m (Stokes n = 0, m = ½; sheet ½, ½; bead ⅕, ⅕); **γ = sheet strength (ch05), c_p/c_v (ch01), Stokes coefficient (ch07); n = electrons in Millikan; m = vertical wavenumber (ch07)** | ch08 → Ch. 9 | `similarity_reduce_sympy(case)`, `similarity_collapse_error(case, n, m, times)` |
| ξ ⚠️ | generic spatial coordinate of the ansatz (8.32) | m | **ξ = axial source coordinate / x − x_s in ch06** | ch08 | (sympy symbol) |
| ω ⚠️ (§8.5) | oscillation angular frequency of the plate | rad/s | **vorticity elsewhere in ch08** (see the trap) | ch07 → | `omega` in `stokes_second_problem`, `stokes_layer` |
| δ_e ⚠️, δ_book | Stokes-layer e-folding depth √(2ν/ω); the book's depth 4√(ν/ω) = 2√2 δ_e | m | amplitude e^{−y/δ_e}, phase lag y/δ_e; 5.91 % at δ_book; crest speed √(2νω) | ch08 → Ch. 13 (Ekman √(2ν/f)) | `stokes_layer(nu, omega)["delta_e"]`, `["delta_book"]` |
| k ⚠️ (Stokes layer) | complex decay rate ±(1 + i)√(ω/2ν) | 1/m | only the decaying root kept; **k = wavenumber (ch07), conductivity (ch01)** | ch08 | `stokes_second_sympy()` |
| Re ⚠️ (sphere), Re_a | 2aU/ν (diameter; (8.52)–(8.53)); ρUa/μ = Re/2 (radius; literature, far-field ratio) | – | C_D = 24/Re; Oseen 3/16 (diameter) = 3/8 (radius); settling "valid" flag Re < 0.1 (our choice) | ch04 → | `stokes_drag_coefficient`, `oseen_drag_coefficient`, `proudman_pearson_drag_coefficient`, `settling_state` |
| r, θ ⚠️ (sphere) | spherical radius; **polar angle from the downstream +x axis** | m, rad | rear stagnation θ = 0, front π; fields NaN for r < a; `frame="body"\|"fluid"` | ch06 → | `stokes_sphere_*`, `oseen_*` |
| ψ ⚠️ (Stokes, §8.6) | axisymmetric stream function about the stream axis, u_r = (1/r²sin θ)∂ψ/∂θ, u_θ = −(1/r sin θ)∂ψ/∂r (6.83) | m³/s | uniform stream ½Ur²sin²θ; fluid frame = body frame minus it | ch06 → | `stokes_sphere_streamfunction`, `oseen_streamfunction` |
| E², E⁴ | Stokes operator ∂²/∂r² + (sin θ/r²)∂/∂θ((1/sin θ)∂/∂θ) and its square | 1/m², 1/m⁴ | ω_φ = −E²ψ/(r sin θ); (E²)²ψ ≠ ∇⁴ψ | ch06 → | `E2(psi, r, theta)`, `E4_residual` (sympy) |
| D ⚠️ (drag), C_D | drag force on the sphere 6πμaU; drag coefficient D/(½ρU²πa²) | N; – | ⅓ pressure + ⅔ friction; **D = integration constant, diffusivity, diameter elsewhere** | ch04 → | `stokes_drag(mu, a, U, parts=)`, `stokes_drag_running` |
| U_t, ρ′ | terminal (settling) velocity 2(ρ′ − ρ)ga²/(9μ); particle density | m/s; kg/m³ | g = G0 = 9.80665 in `core.creeping`; warns when Re is not small | ch08 → Ch. 13 | `terminal_velocity(a, rho_p, rho, mu, g)`, `radius_from_terminal_velocity` |
| E ⚠️ (field), n ⚠️, e | electric field; number of electron charges on the drop; elementary charge | V/m; –; C | ne = 6πμa(U_up + U_t)/E; **E = constant of Ex. 8.7, energy (ch07), E² operator** | ch08 | `millikan_charge`, `synthetic_millikan(seed=0)`, `E_CHARGE` |
| Re_a r/a | far-field inertia/viscous ratio of the Stokes solution → ½Re_a r/a (axis and θ = π/2) | – | crossover r ≈ 2a/Re_a; the book's "r/a ~ 1/Re" is an order of magnitude | ch08 | `inertia_viscous_ratio(r, theta, U, a, nu)` |
| **ch09 — boundary layers, bluff bodies, jets** | | | | | |
| x, y ⚠️ (layer) | distance along the wall from the leading edge (or stagnation point); distance normal to the wall | m | u = ∂ψ/∂y, v = −∂ψ/∂x; jets: x along the jet, y across (free jet symmetric about y = 0) | ch09 → | `x`, `y` |
| U, U∞, U_e(x) ⚠️ | free-stream speed; speed at the edge of the layer (outer inviscid flow at the wall) | m/s | U_e = axⁿ (wedge), U₁/(1 + x/L) (diffuser), 2U sin φ (cylinder); dU_e/dx < 0 decelerating | ch09 → | `U`, `Ue`, `outer_flow(kind)` → `OuterFlow` |
| dp/dx ⚠️ (layer) | pressure gradient imposed by the outer flow, −(1/ρ)dp/dx = U_e dU_e/dx (9.11) | Pa/m | **> 0 adverse**, < 0 favourable (ch08 book sign) | ch08 → | `OuterFlow.dpdx(x, rho)`, `wall_curvature(dpdx, mu)` |
| Re, Re_x, Re_L ⚠️ | overall U∞L/ν (9.6); local U_e x/ν; plate UL/ν | – | six Re in ch09 (P200); cylinder/sphere on the **diameter**; jet Re_x = xu₀/ν; FS Re_x = ax^{n+1}/ν | ch01 → | `Re`, `Re_x`, `Re_L` |
| δ̄ | order-of-magnitude thickness L Re^{−1/2} (9.4) | m | a scale, not a measured thickness | ch09 | `boundary_layer_scales(...)["delta"]` |
| x*, y*, u*, v*, p* ⚠️ | scaled variables (9.6): x/L, y Re^{1/2}/L, u/U, v Re^{1/2}/U, (p − p∞)/ρU² | – | the star is **not** δ*'s star; ε = Re^{−1/2} plays ch08's ε = h/L | ch08 → | `to_bl_variables`, `from_bl_variables` |
| δ₉₉ ⚠️ | height where u = 0.99U_e | m | Blasius 4.910√(νx/U) (page 4.93); **ch08 δ₉₉ = 3.643√(νt) is a diffusion thickness** | ch08 → | `delta_level(y, u, Ue, level=0.99)`, `blasius_delta99(printed=)` |
| δ* ⚠️ | displacement thickness ∫(1 − u/U_e)dy (9.16) | m | lift of the outer streamlines; v∞ = U dδ*/dx; Blasius 1.7208√(νx/U) | ch09 → Ch. 12, 14 | `displacement_thickness(tail=)`, `blasius_delta_star` |
| θ ⚠️ (thickness) | momentum thickness ∫(u/U_e)(1 − u/U_e)dy (9.17) | m | a **length** here, an angle everywhere else; Blasius 0.6641√(νx/U) = 2f″(0)δ; θ₀ in (9.50) = initial θ | ch09 → Ch. 12, 14 | `momentum_thickness`, `blasius_theta`, `thwaites(...)["theta"]` |
| H | shape factor δ*/θ | – | Blasius 2.591, Hiemenz 2.216; H(λ) from the closure | ch09 → Ch. 12 | `thicknesses(...)["H"]`, `thwaites_H(lam, closure)` |
| τ₀, C_f, F_D, C_D ⚠️ | wall shear μ(∂u/∂y)₀; skin friction τ₀/(½ρU²); drag per width on one face; drag coefficient F_D/(½ρU²L) | Pa, –, N/m, – | C_f = 0.664/√Re_x, C_D = 1.328/√Re_L **one side** (`sides=1`); τ₀ ∝ x^{−1/2} integrable | ch08 → | `blasius_wall_shear`, `blasius_skin_friction`, `blasius_drag(sides=)`, `blasius_drag_coefficient(sides=)` |
| ψ, f(η) ⚠️ (similarity) | stream function ψ = Uδ(x)f(η) (9.19) (Blasius), √(νxU_e)f (9.34) (FS), u₀δf (9.59) (free jet), [νCx^{1/2}]^{1/2}f (9.82) (wall jet); f′ = u/U_e | m²/s; – | f(0) = f′(0) = 0 at a wall (9.28); jets f(0) = 0 on the axis (free) | ch04 → | `blasius_profile`, `falkner_skan`, `free_jet_profile`, `wall_jet_profile` |
| η ⚠️ (ch09) | y/δ(x) with the flow's own δ (see the trap) | – | Blasius y√(U/νx) = ch08's y/√(νt) with t = x/U | ch08 → | `eta` |
| f″(0) | dimensionless wall shear of a similarity profile | – | Blasius 0.3320573; Hiemenz 1.232588; 0 at the FS fold; wall jet f∞³/72 | ch09 | `blasius_constants()["fpp0"]`, `falkner_skan(m)["fpp0"]` |
| n ⚠️ (FS), m (code), β | Falkner–Skan exponent in U_e = axⁿ; code name; Hartree β = 2m/(m + 1) | – | n > 0 favourable, n < 0 adverse; fold n = −0.090429 (β = −0.198838); **ch08 n, m = similarity exponents** | ch09 → Ch. 11 | `falkner_skan(m, branch=)`, `falkner_skan_separation()` |
| a ⚠️ (wedge) | coefficient of U_e = axⁿ | m^{1−n}/s | **a = radius (pipe, sphere), row spacing (Kármán) elsewhere** | ch09 | `falkner_skan_fields(x, y, m, a, nu)` |
| I_δ, I_θ | ∫(1 − f′)dη, ∫f′(1 − f′)dη on a similarity profile | – | δ* = I_δδ, θ = I_θδ; H = I_δ/I_θ | ch09 | `falkner_skan_state(m)["I_delta"]`, `["I_theta"]` |
| λ ⚠️ (Thwaites) | Holstein–Bohlen parameter (θ²/ν)dU_e/dx (9.44) | – | λ < 0 adverse; stagnation 0.075 (Thwaites) / 0.0855 (exact); separation −0.090 (book) / −0.068148 (exact FS); **λ = eigenvalue in C10** (trap) | ch09 | `holstein_bohlen(theta, Ue_x, nu)`, `LAMBDA_SEP_BOOK`, `LAMBDA_SEP_FS` |
| l(λ), L(λ) ⚠️ | shear correlation τ₀θ/(μU_e) (9.45); L = 2l − 2(2 + H)λ (9.48) ≈ 0.45 − 6λ (9.49) | – | **capital L = body length elsewhere; l ≠ 1**; exact FS closure vs `"white"` fit | ch09 | `thwaites_l(lam, closure)`, `thwaites_L`, `thwaites_closure_table` |
| u_yy,wall | wall curvature (dp/dx)/μ | 1/(m s) | sign = sign of dp/dx (R17); inflection when > 0 | ch09 → Ch. 11 | `wall_curvature(dpdx, mu)`, `profile_inflection(y, u)` |
| x_sep, φ_sep ⚠️ | separation station (τ₀ = 0 or λ = criterion); separation angle from the **forward** stagnation point | m; deg | `phi_deg` in every BB function; 82° laminar, 125° turbulent (rounded) | ch09 → Ch. 14 | `separation_point(x, tau0)`, `thwaites_cylinder_separation(lam_sep)` |
| w, σ ⚠️ (marching) | von Mises unknown w = u² at fixed ψ; mapped coordinate ψ = ψ_max σ² | m²/s², – | **σ = surface tension, core radius, D14 coefficient elsewhere** | ch09 | `march_boundary_layer(..., ny, order)` |
| C_p, C_b | pressure coefficient (p − p∞)/(½ρU²); base (wake) pressure coefficient | – | ideal cylinder 1 − 4 sin²φ; C_b ≈ −1.2 (subcritical), −0.6 (supercritical), illustrative | ch06 → | `cp_ideal_cylinder`, `separated_cp(phi_deg, phi_sep_deg, cp_base)` |
| C_D,p | pressure (form) drag coefficient ½∮C_p cos φ dφ per unit span on the diameter | – | model sin φ_s(1 − (4/3)sin²φ_s − C_b); 0.884 (82°, −1.2), 0.578 (125°, −0.6) | ch09 → Ch. 14 | `pressure_drag_from_cp`, `separated_pressure_drag` |
| d, St, f ⚠️ (shedding) | cylinder diameter; Strouhal fd/U; shedding frequency | m, –, Hz | St ≈ 0.2 with the **cyclic** f (trap); `omega_rad` = 2πf; **f = Coriolis parameter, friction factor, similarity profile elsewhere** | ch04 → | `shedding_frequency(U, d, St)`, `strouhal_of_re` |
| a, b, Γ ⚠️ (street) | streamwise vortex spacing, distance between the rows, circulation magnitude | m, m, m²/s | stable b/a = arccosh(√2)/π = 0.28055; street speed (Γ/2a)tanh(πb/a); `offset` 0.5 staggered, 0 facing | ch05 → | `karman_street_ratio`, `karman_street_velocity(a, b, Gamma)`, `karman_street_spectrum(b_over_a, offset, k)` |
| γ, σ ⚠️ (D14) | ½ − sech²(πb/a); sinh(πb/a)/cosh²(πb/a) | – | eigenvalues (πΓ/2a²)(±γ ± iσ); see the trap | ch09 | `karman_street_growth_closed`, `ch09.karman_street_sympy` |
| λ ⚠️ (street eigenvalue) | eigenvalue of the linearised street, growth Re λ | 1/s | 0.48 Γ/a² at b/a = 0.15; π/4 for facing rows; code docstring calls it σ | ch09 → Ch. 11 | `karman_street_spectrum`, `karman_street_growth` |
| J ⚠️ | momentum flux per unit span ρ∫u²dy (9.58) | N/m | conserved in the free jet; ∝ x^{−1/4} in the wall jet | ch09 → Ch. 12 | `jet_momentum_flux(y, u, rho)`, `free_jet(x, y, J, rho, nu)` |
| u₀(x) ⚠️ | centreline speed (free jet, ∝ x^{−1/3}); velocity **scale** Cx^{−1/2} (wall jet, not the peak: peak ≈ 0.079 of it at f∞ = 1) | m/s | u₀ = (J²/(C²ρ²νx))^{1/3} (9.62) | ch09 | `free_jet_centreline`, `wall_jet(...)` |
| C ⚠️ (jets) | free jet ∫f′²dη = 4√6/3 (9.72) (`C_fj`); wall jet dimensional u₀ = Cx^{−1/2} | –; m^{3/2}/s | **C also Blasius proportionality, integration constants**; wall-jet C and f∞ are one physical constant (gauge f → λf(λη)) | ch09 | `free_jet_constants()["C"]`, `wall_jet_constants(rho, nu, Psi=, mdot=, x=)` |
| ṁ | mass flux per unit span ρ∫u dy | kg/(m s) | free jet (36Jρ²νx)^{1/3} ∝ x^{1/3} (9.73); wall jet ∝ x^{1/4} (9.84) — entrainment | ch09 → Ch. 12, 13 | `free_jet_mass_flux`, `wall_jet_mass_flux` |
| h₉₉ | jet half-width where u = 0.01u₀ | m | 7.3319[Cρν²x²/J]^{1/3} (page 5.6152 = 4 % point) | ch09 | `free_jet_halfwidth(x, J, rho, nu, level=0.01, printed=)` |
| Ψ ⚠️, f∞, K₁ | wall-jet invariant ∫u(∫_y^∞u²dy′)dy (9.80); f(∞) of the wall jet; ∫f′(∫f′²)dη = 1/40 at f∞ = 1 | m⁵/s³, –, – | Ψ = C²νf∞⁴/40; f″(0) = f∞³/72; **ψ = stream function** | ch09 | `wall_jet_invariant(y, u)`, `wall_jet_K1()`, `wall_jet_integrals(f_inf)` |
| u_e, u(z), R ⚠️ (teacup) | swirl speed above the bottom layer; swirl in the layer; radius | m/s, m/s, m | net inward force ρ(u_e² − u²)/R ≥ 0, zero in the core | ch09 → Ch. 13 | `secondary_flow_radial_force(u_inviscid, u_layer, R, rho)` |
| **ch10 — computational fluid dynamics** | | | | | |
| T, T_i^n | transported scalar (temperature, concentration) at node i and time level n (Fig. 10.1) | any (K, kg/m³) | node-based: T_i^n ≈ T(x_i, t_n), x_i = iΔx, t_n = nΔt | ch01 → | `T`, `T0` |
| u (1-D), D ⚠️ | constant advection speed; diffusivity of (10.1) | m/s; m²/s | u may be negative (upwind side from sign(u)); D ≥ 0; **not** the strain-rate tensor **D**[u] | ch10 | `u`, `D` |
| g ⚠️, q | Dirichlet value T(0, t) = g; Neumann slope ∂T/∂x(L, t) = q (10.2) | unit of T; unit of T per m | q = 0 insulated end (ghost node); g ≠ body force | ch10 | `g=`, `q=` |
| Δx, Δt, h, n ⚠️ | grid spacing; time step; element length (FE); number of cells/elements | m; s; m; – | h also "the step" of a stencil test; superscript n = time level | ch01 → | `dx`, `dt`, `h`, `n`, `n_el` |
| α ⚠️, β ⚠️ (FTCS) | FTCS numbers α = uΔt/(2Δx), β = DΔt/Δx² (10.11) | – | the ½ in α comes from the centred difference; β = ch01's r | ch10 | `alpha`, `beta`, `ftcs_coefficients` |
| C | Courant number uΔt/Δx = 2α | – | signed in `cfl_number`; stability uses \|C\| | ch10 → Ch. 13, 15 | `C`, `cfl_number(u, dt, dx)` |
| R ⚠️, R_cell | global Péclet uL/D (10.87); cell Péclet uΔx/D (10.31) | – | wiggles iff R_cell > 2 (centred); R < 0 puts the layer at x = 0; R = radius in ch08/ch09 | ch10 → Ch. 13 | `R`, `R_cell`, `cell_peclet(u, dx, D)` |
| δ ⚠️ (layer) | convection–diffusion layer thickness, δ/L = O(1/R) (10.89) | m | T = e⁻¹ one δ from the wall; ≠ Kronecker δ_AB; ≠ ch09 similarity length | ch10 | `cd_layer_thickness(R, L, level)` |
| E ⚠️ (truncation) | truncation error of a scheme (10.16)–(10.17) | unit of T per s | "exact minus scheme" over Δt; ≠ internal energy (ch01), ≠ x-flux **E** of (10.100) | ch10 | `truncation_error_sympy`, `truncation_terms` |
| e, K_e, a, b ⚠️ | error array T − T_exact (10.14); constant and rates in ‖e‖ ≤ K_eΔx^aΔt^b (10.15) | unit of T; –; – | K_e ≠ stiffness **K**; a, b ≠ FE α_t, β_t | ch10 | `error_norm(num, exact, kind)` |
| ξ, ξ̂ⁿ | round-off disturbance (10.18); its Fourier amplitude (book gⁿ(k)) | unit of T | obeys the same linear scheme (10.19) | ch10 | `propagate_error` |
| k ⚠️, θ ⚠️ | Fourier wavenumber in the book's e^{iπkx} (half-waves per unit length); phase per cell θ = kπΔx | 1/m; rad | θ = π is the zigzag; `convention="standard"`: e^{ikx}, θ = kΔx; θ ≠ Θ-scheme fraction, ≠ ch09 momentum thickness | ch10 → Ch. 11 | `theta`, `fourier_mode(x, k, convention)` |
| G ⚠️ | amplification factor per step g^{n+1}/gⁿ (10.23)–(10.24) | – (complex) | stable iff \|G\| ≤ 1 for every θ; ≠ G_A (weak-form row), ≠ ch02 velocity gradient G, ≠ ch08 `G=` = −dp/dx | ch10 → Ch. 11 | `amplification_factor(theta, alpha, beta, scheme)` |
| c_num/u | numerical phase speed over the true one, −arg G/(Cθ) | – | < 1: short waves lag (dispersion) | ch10 | `phase_error(scheme, C, theta)` |
| D_num | numerical diffusivity: steady upwind uΔx/2 = 0.5R_cell D (10.94); unsteady \|u\|Δx(1 − C)/2 | m²/s | negative for C > 1 (anti-diffusion) | ch10 → Ch. 12, 13 | `numerical_diffusivity(u, dx, D, scheme, C)` |
| r | root of the discrete recurrence: centred (1 + R_cell/2)/(1 − R_cell/2), upwind 1 + R_cell, downwind 1/(1 − R_cell) | – | r < 0 ⇒ node-to-node alternation | ch10 | `discrete_root(R_cell, scheme)` |
| w, 𝒮, V, H¹ | test function; trial space (Dirichlet built in, (10.32)); test space (w(0) = 0, (10.33)); finite-slope-energy space | – | essential BC in 𝒮, natural BC in the weak form | ch10 | `weak_residual(T_fn, w_fn, …)` |
| a(w, v) | bilinear form (10.42) u∫v_xw dx + D∫v_xw_x dx | (unit of w)(unit of v)/s | not symmetric when u ≠ 0 | ch10 | `bilinear_form(w, v, x, u, D)` |
| N_A, N₀, d_A, c_A | hat (shape) function of node A (10.59); extra function at the Dirichlet node; nodal unknowns T^h(x_A); test coefficients (10.45) | –; –; unit of T; – | N_A(x_B) = δ_AB; partition of unity | ch10 | `hat(x, x_nodes, A)`, `hat_basis`, `interpolate` |
| **M**, **K**, **F** ⚠️ | FE mass, stiffness, load (10.55)–(10.58) | M_AB [m]; K_AB [m/s]; F_A [unit of T · m/s] (1-D, dimensionless hats) | **M** ≠ Ma; **F** ≠ y-flux; K not symmetric for u ≠ 0 | ch10 | `assemble_1d(...)` → M, K, F |
| m^e, k^e, f^e, ξ (FE) ⚠️ | element blocks (10.74)–(10.77); parent coordinate ξ ∈ [−1, 1] (10.64) | – | m^e = (h/6)[[2,1],[1,2]]; local a = 1, 2 ↔ A = e − 1, e (R2); ξ ≠ disturbance ξ | ch10 | `element_matrices_linear(h, u, D)`, `connectivity(n_el)`, `map_to_parent` |
| U, E, F (conservation form) ⚠️ | state (ρ, ρu, ρv) and x-, y-fluxes of U_t + E_x + F_y = 0 (10.100) | per unit volume | ≠ FE **F**; E ≠ truncation error | ch10 → Ch. 15 | `maccormack_step(U, flux_E, flux_F, …)`, `ns_fluxes` |
| c ⚠️, Ma, p = c²ρ | sound speed; Mach number U/c; isothermal state (10.99) | m/s; –; Pa | non-dimensional p = ρ/Ma²; **book prints M** in (10.151)–(10.155); c₁–c₅ (10.109) are coefficients, not speeds | ch01 → | `c`, `Ma`, `pressure_isothermal(rho, c)` |
| σ (MacCormack) ⚠️ | safety factor of the time-step rules (10.110), (10.155) | – | ours 0.8 (book value private); ≠ surface tension (ch07), ≠ D14 coefficient (ch09) | ch10 | `sigma=` |
| Re_Δ | grid (mesh) Reynolds number in (10.110) | – | zero speeds left out (DEVIATION) | ch10 | inside `maccormack_dt` |
| ρ′ | density perturbation ρ − 1 stored instead of ρ (practice, N90) | – (non-dim) | round-off control; ≠ ch04 perturbation | ch10 | `perturbation=True` |
| A₁, A₂, [A₁, A₂] | split operators (10.112); commutator A₁A₂ − A₂A₁ | – | splitting error (Δt²/2)[A₁, A₂] per step | ch10 → Ch. 13 | `split_linear_system`, `marchuk_yanenko` |
| Θ, α_Θ, β_Θ ⚠️ | Glowinski Θ-scheme fraction and weights (10.129)–(10.133) | – | α_Θ + β_Θ = 1, β_Θ = Θ/(1 − Θ); three substeps of length ΘΔt, (1 − 2Θ)Δt, ΘΔt; second order only at Θ = 1 − 1/√2 (the book prints five digits; we compute it) | ch10 | `theta_scheme_linear(..., theta_split=1 - 1/np.sqrt(2), alpha_split=None)` |
| u*, u^{n+1/2}, p^{n+1} | predicted (divergent) velocity; projected pressure (10.116)–(10.118) | non-dim | ∇²p^{n+1} = ∇·u^{n+1/2}/Δt; correction −Δt∇p curl-free | ch10 → Ch. 11–13 | `predictor`, `solve_pressure`, `project` |
| u_{i+1/2, j}, p_{i, j} | staggered face velocity; cell-centre pressure (Fig. 10.4) | non-dim | `u[j, i+1]` ↔ u_{i+1/2, j}; shapes p[ny, nx], u[ny, nx+1], v[ny+1, nx] | ch10 → Ch. 13 | `MacGrid`, `face_coordinates` |
| ∇_d·, ∇²_d | discrete divergence (10.123) and discrete Laplacian of p (10.124) | 1/s; – | Poisson matrix singular (constant null vector); pin the mean | ch10 | `divergence`, `pressure_poisson_matrix(g, pin)` |
| ψ ⚠️ (cavity), ψ_min | stream function of the computed field (corner-based); its minimum (primary eddy) | non-dim | ψ = 0 on the walls; u = ∂ψ/∂y; ψ_min < 0 for the clockwise eddy (lid in +x) | ch04 → | `streamfunction(u, v, g)`, `primary_vortex_centre(psi, g)` |
| **A**, **B**, **M**_p, β_h ⚠️ | velocity block, divergence (pressure) block, pressure mass matrix of the saddle system (10.137); discrete inf–sup constant | – | β_h² = λ_min(**B**ᵀ**A**⁻¹**B**, **M**_p) on zero-mean pressures; P2–P1 ≈ 0.366; β_h ≠ FTCS β | ch10 → Ch. 11 | `infsup_constant(n, pair)` |
| φ_a, ψ_b, ζ ⚠️ | P2 velocity shapes (10.185); P1 pressure shapes (10.187); barycentric ζ = 1 − ξ − η | – | ψ_b ≠ stream function; ζ ≠ vorticity/displacement | ch10 | `p2_shape(xi, eta)`, `p1_shape` |
| J ⚠️ (FE) | Jacobian x_ξy_η − x_ηy_ξ (10.198) | – | = 2 × area for a straight triangle; ≠ ch09 jet momentum flux J | ch10 | `jacobian(xe, ye, xi, eta)` |
| α_t, β_t | time-derivative weights (10.163): backward Euler (1, 0), trapezoidal (2, 1) | – | ≠ FTCS α, β; ≠ Θ-scheme | ch10 | `time_derivative(..., alpha_t, beta_t)`, `march_unsteady(alpha_t=2, beta_t=1)` |
| u′, p′ ⚠️ (Newton) | Newton corrections (10.164) | non-dim | p′ ≠ ch04/ch07 perturbation pressure | ch10 | inside `newton_solve` |
| St, f_s, τ̄ | Strouhal number f_s d/U = 1/τ̄ (non-dimensional period τ̄) | –; Hz; – | **cyclic** frequency (book n); confined (W = 5d) St 0.2054 is qualitative | ch09 → | `strouhal_from_period`, `dominant_frequency` |
| p (order), r (ratio), f₀, GCI | observed order; refinement ratio; Richardson extrapolate; grid convergence index (F_s = 1.25) | – | p ≠ pressure; NaN when rᵖ − 1 = 0 | ch10 → all | `grid_convergence_index(f1, f2, f3, r)`, `richardson_three` |
| **ch11 — instability** | | | | | |
| σ ⚠️, σ_r, σ_i | complex growth rate of a normal mode e^{σt} (11.1); its real part (growth) and imaginary part (frequency) | 1/s (or κ/d², ν/d², U₀/L when scaled) | σ_r > 0 unstable; **σ = −ikc**; ≠ surface tension (σ_s here), ≠ vortex core radius (ch03), ≠ MacCormack safety factor (ch10) | ch11 → Ch. 12, 13 | `sigma`, `ST.sigma_from_c(K, c)`, `stability_class`, `marginal_type` |
| c, c_r, c_i | complex phase speed of e^{ik(x − ct)}; phase speed; growth is kc_i | m/s (or U₀) | c_i > 0 unstable; unstable c lies in Howard's semicircle on [U_min, U_max] | ch07 (real c) → ch11 (complex) → Ch. 13 | `c`, `ST.c_from_sigma`, `in_howard_semicircle` |
| k, m, K ⚠️ | streamwise (or axial) wavenumber; spanwise wavenumber; K = (k, m, 0), ∣K∣ horizontal wavenumber magnitude (Bénard K in units of 1/d) | 1/m (or 1/L, 1/d) | k > 0; Squire: k̄ = √(k² + m²) | ch07 → ch11 | `k`, `m`, `K` |
| σ_s | surface tension in the KH relation (Ex. 11.1) | N/m | adds σ_s k³ restoring; ch01/ch07 wrote σ | ch11 | `surface_tension` |
| ΔU, U₁, U₂, ρ₁, ρ₂ | velocity jump U₁ − U₂; upper (1) and lower (2) stream speeds and densities | m/s; kg/m³ | ρ₂ > ρ₁ bottom-heavy; growth independent of the sign of ΔU, c_r not | ch11 | `kh_phase_speed(k, U1, U2, rho1, rho2, …)`, `kh_min_shear` |
| Ra ⚠️ | Rayleigh number gαΓd⁴/(κν) (11.21), Γ = −dT̄/dz | – | > 0 heated from below (§11.4); **signed with dT̄/dz in §11.5** (negative when heated from below); Ra_c 1707.76 rigid–rigid, 1100.65 rigid–free, 27π⁴/4 free–free | ch11 → Ch. 13 | `rayleigh_number(alpha, dT, d, kappa, nu, g=G0)`, `thermal_rayleigh_signed` |
| Γ (Bénard) ⚠️ | conduction temperature gradient −dT̄/dz = ΔT/d (11.21), (11.24) | K/m | positive heated from below; **opposite in sign to Kundu's Ch. 1 Γ ≡ dT/dz**; equal to the meteorological lapse rate | ch11 §11.4 only | `gamma_conventions(dT, d)`; functions take `dT` or `dTdz` |
| d | layer depth (Bénard, double diffusion); gap width R₂ − R₁ (Taylor) | m | the length scale of Ra, Rs, Ta | ch11 | `d` |
| W, T̂, D | scaled vertical-velocity and temperature amplitudes of (11.36)–(11.37); D = d/dz | – | W ≡ (Γd²/κ)ŵ; rigid wall W = DW = T̂ = 0; free wall W = D²W = T̂ = 0 | ch11 | inside `benard_growth_rate`, `benard_eigenfunction` |
| Pr ⚠️ | Prandtl number ν/κ | – | does not move the Bénard margin (only σ); in the Lorenz system it is the literature's σ | ch04 → ch11 | `Pr` |
| a² | π² + K² (free–free Bénard), the eigenvalue of −(D² − K²) on the first sine mode | – | Ra = a⁶/K² (11.44) | ch11 | inside `benard_free_free_Ra`, `benard_free_free_sigma` |
| S, β ⚠️, κ_s, τ | salinity; haline contraction coefficient; salt diffusivity; τ = κ_s/κ ≈ 0.01 | –; 1/(salinity unit); m²/s; – | ρ = ρ₀[1 − α(T − T₀) + β(S − S₀)]; τ ≠ stress, ≠ torsion | ch11 → Ch. 13 | `dSdz`, `beta_S`, `kappa_s`, `linear_eos` |
| Rs, Rs′ | salinity Rayleigh numbers gβd⁴(dS/dz)/(νκ_s) and gβd⁴(dS/dz)/(νκ) | – | fingers when Rs − Ra > 27π⁴/4 (11.46) with §11.5's signed Ra | ch11 | `salinity_rayleigh`, `salinity_rayleigh_prime`, `double_diffusive_margin` |
| R_ρ | density ratio αT_z/(βS_z) | – | shown in the E4 status when the column is density-stable | ch11 → Ch. 13 | `salt_finger_regime(...)["R_rho"]` |
| Ω₁, Ω₂, R₁, R₂, μ ⚠️ | angular speeds and radii of the inner (1) and outer (2) cylinder; μ = Ω₂/Ω₁ | rad/s; m; – | **`mu` is the speed ratio, not viscosity**; book α = μ − 1; Rayleigh's line μ = (R₁/R₂)² | ch08 → ch11 | `mu`, `couette_rayleigh_line`, `taylor_critical(mu)` |
| A, B (Couette) | circular Couette constants U_φ = AR + B/R (ch08 (8.9)) | 1/s; m²/s | −4AΩ₁ > 0 where Γ² falls outward | ch08 → ch11 | `core.laminar.circular_couette`, inside `taylor_number` |
| Ta | Taylor number −4AΩ₁d⁴/ν² (11.52) | – | > 0 only beyond Rayleigh's line; narrow gap, inner only: 2(Ω₁R₁d/ν)²(d/R₁); Ta_c ≈ 1708/(½(1 + μ)) (11.54) | ch11 | `taylor_number`, `taylor_number_narrow_inner`, `taylor_critical_approx` |
| Γ (ring) ⚠️ | circulation 2πrU_θ of a fluid ring | m²/s | unstable where Γ² decreases outward (Rayleigh); → Ch. 13 with absolute angular momentum | ch05 → ch11 → Ch. 13 | `ring_interchange_energy`, `rayleigh_circulation_criterion` |
| N² | buoyancy frequency squared −(g/ρ₀)dρ̄/dz, ch07 **(7.127)** | 1/s² | > 0 stable; may be negative in `gradient_richardson` | ch01 → ch07 → ch11 | `N2` (callable or array) |
| Ri | gradient Richardson number N²/(dU/dz)² (11.66) | – | Ri > ¼ **everywhere** ⇒ stable (11.67); Ri < ¼ somewhere only allows instability; ±∞ at U′ = 0 with the sign of N² | ch04 (named) → ch11 → Ch. 13 | `gradient_richardson`, `miles_howard_stable` |
| J ⚠️, R ⚠️ (tanh layer) | bulk Richardson number (Ri at the centre) of U = tanh z, N² = J sech²(Rz); R = shear/density thickness ratio | – | neutral curve J = k(1 − k) at R = 1; ≠ ch09 J, ≠ ch10 J | ch11 | `tg_growth(k, J, R=1.0)`, `tg_tanh_neutral_J`, `richardson_profiles` |
| ψ̂, φ ⚠️, F, Q | stream-function amplitude; φ = ψ̂/(U − c)^{1/2} (Miles–Howard) or the OS/Rayleigh amplitude; F = ψ̂/(U − c); Q = ∣F′∣² + k²∣F∣² ≥ 0 | – | û = φ′, v̂ = −ikφ in §11.8; φ = velocity potential in §11.3 | ch11 | `phi` (eigenvectors of `os_mode`, `rayleigh_eigs`) |
| U, U′, U″, U_I, y_c | basic parallel flow and its derivatives; speed at the inflection point; critical level where U = c_r | U₀; U₀/L; U₀/L²; –; L | Rayleigh: U″ changes sign; Fjørtoft: U″(U − U_I) < 0 somewhere | ch11 → Ch. 13 (U″ − β) | `U`, `Up`, `Upp` callables; `U_I`; `critical_layer` |
| Re ⚠️, Re_c, k_c | Reynolds number U₀L/ν of the chosen flow; its critical value and wavenumber | – | **L per flow**: half-width (Poiseuille), δ* (Blasius), L (tanh, sech²); Squire Re̅ = kRe/k̄ | ch04 → ch11 | `Re`; `poiseuille_critical`, `blasius_critical`, `bickley_critical`, `table_11_1` |
| ω, F (frequency) ⚠️ | wave frequency kc_r; Blasius frequency parameter ων/U∞² | U₀/L; – | F ≈ 2.3e-4 at Re_c; ≠ FE load, ≠ Howard's F | ch11 | `blasius_neutral_curve(in_frequency=True)` |
| E, P, Λ | disturbance kinetic energy per wavelength; production −∫⟨uv⟩U′dy; viscous dissipation (11.88) | per unit ρU₀²L | dE/dt = P − Λ = 2kc_iE; P > Λ inside the neutral curve | ch11 → Ch. 12 | `ST.disturbance_energy_budget` → `production`, `dissipation`, `dEdt`, `ratio`, `residual` |
| y_max, s, δ (path) ⚠️ | truncation box of an unbounded profile; scale of the tan map; depth of the complex path below/above the real axis | L | y_max = max(profile box, 12/k); s = min(0.5, max(0.035, 0.025/k)); δ inside U's strip of analyticity; ≠ boundary-layer δ | ch11 → every unbounded eigenproblem | `ST.decay_box`, `decay_map_scale`, `decay_box_map_scale`; `delta=` in `rayleigh_eigs_contour` |
| N (degree) ⚠️ | Chebyshev polynomial degree (N + 1 nodes, descending) | – | ≠ buoyancy frequency; errors fall exponentially then rise with round-off (N ≳ 60 for Bénard) | ch11 | `N=` in every solver |
| X, Y, Z, r ⚠️, b ⚠️ | Lorenz amplitudes (roll speed, horizontal temperature contrast, distortion of the mean profile); r = Ra/Ra_c(k); b = 4π²/(π² + k²) | – | time in units of d²/((π² + k²)κ); pitchfork at r = 1, Hopf at r_H = Pr(Pr + b + 3)/(Pr − b − 1) | ch11 → Ch. 12, 13 | `lorenz_rhs(t, s, Pr, r, b)`, `lorenz_r`, `lorenz_b`, `lorenz_hopf_r` |
| δ₀, λ (Lyapunov) ⚠️ | initial separation of two runs; largest Lyapunov exponent (slope of ln∣δ∣) | –; 1/time | λ ≈ 0.9 at r = 28 (**qualitative**); ≠ wavelength, ≠ ch09 Thwaites λ, ≠ ch10 eigenvalue | ch11 | `lorenz_separation`, `lorenz_largest_lyapunov`, `lorenz_predictability_time` |
| A (logistic), δ_F | control parameter of x_{n+1} = Ax_n(1 − x_n); Feigenbaum ratio | – | period doubling at 3, 1 + √6, …; δ_F → 4.6692 | ch11 | `logistic_map`, `period_doubling_points`, `feigenbaum_estimate` |
| **ch12 — turbulence** | | | | | |
| ũ_i, U_i, u_i ⚠️; T̃, T̄, T′ | total, mean and fluctuating velocity (12.24); the same for temperature | m/s; K | **lower case = fluctuation**; $\overline{u_i}=0$; T̄, T′ potential temperature in buoyancy terms | ch12 → Ch. 13 | `samples` (members on axis 0), `mean`, `fluct`; `TS.reynolds_decompose` |
| $\overline{(\ )}$, N (members) ⚠️ | ensemble average over N realizations (12.1); time average over Δt (12.2); volume average (12.3) | – | ensemble unless the name says `time_`/`volume_`; N ≠ buoyancy frequency (`N2`) | ch12 → | `ensemble_average`, `time_average(t, u, window)`, `volume_average`; `n_members` |
| R_ij ⚠️, r_ij, R₁₁(τ), r₁₁(τ) | correlation tensor (12.12); correlation coefficient (12.14)–(12.15); autocorrelation in lag form (12.17) and its coefficient | m²/s²; – | −1 ≤ r ≤ 1; R₁₁ even in τ; R_ij(−τ) = R_ji(τ); ≠ ch02–ch03 rotation tensor | ch12 → Ch. 13 | `correlation`, `correlation_coefficient`, `autocorrelation`, `cross_correlation` |
| τ (lag) ⚠️, t_c, τ_c | time lag; correlation time (first zero of r); memory time of the Ornstein–Uhlenbeck test signal | s | τ ≥ 0 returned; ≠ stress τ̄, τ₀ | ch12 → | `lag`, `correlation_time`, `tau_c` |
| Λ_t, Λ_f, Λ_g ⚠️ | integral time scale ∫₀^∞r₁₁dτ (12.18); longitudinal and transverse integral length scales (12.39) | s; m | Λ_g = Λ_f/2 (3-D isotropic); `upto="first_zero"` is not (12.18) when r has a negative lobe — use `"all"`; ≠ ch11 dissipation Λ | ch12 → Ch. 13 | `integral_scale(lag, r, upto=)`, `Lambda_t` |
| λ_t, λ_f, λ_g, λ_T ⚠️ | Taylor microscales: −2/r″(0) in time (12.19), longitudinal and transverse (12.39), generic λ_T in (12.52) | s; m | λ_g = λ_f/√2; does not exist for an OU signal (cusp) | ch12 | `taylor_microscale`, `lambda_f`, `lambda_g` |
| S_e(ω), S₁₁(k₁), E(K), S_T ⚠️ | frequency spectrum (12.20); one-dimensional longitudinal wavenumber spectrum (12.45); three-dimensional energy spectrum (book S(K)); temperature spectrum (12.113) | m²/s; m³/s²; m³/s²; K² m | **two-sided**, 1/2π forward, ∫S = variance; E on K ≥ 0 with ∫E dK = ē; `one_sided=True` doubles | ch12 → Ch. 13 | `spectrum_from_correlation`, `periodogram`, `inertial_spectrum_1d/3d`, `model_spectrum`, `scalar_spectrum` |
| ω, k₁, K, U₀ | angular frequency; streamwise wavenumber; wavenumber magnitude; probe or sweeping speed of the frozen-field swap k₁ = ω/U₀ | rad/s; rad/m; rad/m; m/s | S(k) = U₀S(ω); frozen-field error ∝ u_rms/U₀ | ch12 | `omega`, `k1`, `K`, `U0`; `taylor_frozen`, `frequency_to_wavenumber_spectrum` |
| $\overline{u_iu_j}$, −ρ₀$\overline{u_iu_j}$, τ̄_ij | velocity covariance; Reynolds stress; mean stress tensor $-P\delta_{ij}+2\mu\bar S_{ij}-\rho_0\overline{u_iu_j}$ of (12.30) | m²/s²; Pa; Pa | $\overline{uv}$ < 0 where dU/dy > 0; **`uv_plus` = `minus_uv_plus` = −$\overline{uv}$/u_*²** | ch12 → Ch. 13, 14 | `velocity_covariance`, `reynolds_stress`, `mean_stress_tensor`, `uv_plus` |
| ē ⚠️, Ē, $\overline{u^2}$ | turbulent kinetic energy ½$\overline{u_i^2}$ per unit mass (the "k" of k–ε); mean-flow kinetic energy ½U_i²; variance of **one** component in §12.6 | m²/s² | ē = (3/2)$\overline{u^2}$ when isotropic; ≠ ch01 internal energy e | ch12 → Ch. 13 | `e`, `E_mean`, `u2`; `turbulent_kinetic_energy` |
| ε̄, ε̄_T | dissipation rate of ē, 2ν$\overline{S'_{ij}S'_{ij}}$ (12.42); dissipation of ½$\overline{T'^2}$ | m²/s³; K²/s | ε̄ ~ (ΔU)³/L (12.49), independent of ν | ch04 ε → ch12 → Ch. 13 | `eps`, `eps_T`; `dissipation_rate`, `dissipation_isotropic`, `dissipation_outer_scaling` |
| P (production), g α$\overline{wT'}$ | shear production −$\overline{u_iu_j}$∂U_i/∂x_j of (12.47); buoyant production | m²/s³ | production > 0 takes energy from the mean flow; buoyant term > 0 for an upward heat flux | ch11 (11.88) → ch12 → Ch. 13 | `shear_production`, `buoyant_production`, `tke_budget`, `mean_energy_budget` |
| f(r), g(r) ⚠️ | longitudinal and transverse correlation coefficients (12.38) | – | g = f + (r/2)f′ in 3-D (12.41); g = d(rf)/dr in 2-D | ch12 | `f`, `g_corr`; `transverse_from_longitudinal(dim=)` |
| η ⚠️, u_K, τ_η, η_T | Kolmogorov length (ν³/ε̄)^{1/4}, velocity (νε̄)^{1/4} (12.50), time (ν/ε̄)^{1/2}; Batchelor scale η(κ/ν)^{1/2} | m; m/s; s; m | ηu_K/ν = 1; η/L ~ Re_L^{−3/4} (12.51); ≠ surface elevation, ≠ similarity variable | ch12 → Ch. 13 | `kolmogorov_scales(nu, eps)` → (η, u_K, τ_η); `batchelor_scale` |
| ΔU, L, Re_L, R_λ | outer velocity and length scales; ΔUL/ν; λ$\sqrt{\overline{u^2}}$/ν | m/s; m; –; – | name which λ (λ_f or λ_g: a factor √2) | ch12 | `dU`, `L`, `Re_L`; `taylor_reynolds_number` |
| C, C₁ ⚠️ | Kolmogorov constants of E(K) = Cε̄^{2/3}K^{−5/3} and of $S_{11}=C_1\bar\varepsilon^{2/3}k_1^{-5/3}$ (12.54) | – | C ≈ 1.5; one-sided C₁ = (18/55)C; `C1` argument always one-sided | ch12 | `KOLMOGOROV_C`, `kolmogorov_constants` |
| J_s, Ṁ_s, U_CL, δ(x), ξ, F, G, Ψ | momentum flux per span (12.62); slot-fluid mass flux; centreline speed; layer width; similarity variable y/δ (= y/x for the jet); velocity and stress profile functions; stress amplitude | N/m; kg/(m s); m/s; m; –; –; –; m²/s² | J_s invariant; δ ∝ x, U_CL ∝ x^{−1/2}; ξ½ = half-width in ξ; J ≠ ch11 bulk Richardson number, Ψ ≠ ch09 wall-jet invariant | ch09 J → ch12 | `Js`, `Ms`, `xi`, `xi_half`, `C5`, `C3`; `plane_jet_*` |
| v_e, V̇ | entrainment velocity ½dV̇/dx = −V(+∞); volume flux per span | m/s; m²/s | V̇ ∝ x^{1/2}; V does not vanish at the jet edge (slip #18) | ch12 → Ch. 13 plumes | `plane_jet_entrainment_velocity`, `plane_jet_volume_flux` |
| u_*, l_ν, y⁺, U⁺, Re_τ = δ⁺ | friction velocity √(τ₀/ρ) (12.81); viscous length ν/u_*; wall units; friction Reynolds number δu_*/ν | m/s; m; –; –; – | y from the wall; sublayer y⁺ < 5, buffer 5–30, log layer to y/δ ≈ 0.15 | ch12 → Ch. 13, 14 | `u_star`, `friction_velocity`, `viscous_length`, `wall_units`, `friction_reynolds_number`, `layer_name` |
| κ ⚠️, B, A, Π ⚠️ | von Kármán constant and additive constant of (12.88); constant of the defect form (12.89); Coles' wake strength | – | **required keywords**; a fitted pair belongs to its window; `LOG_LAW_CONSTANTS` holds cited presets; ≠ thermal diffusivity, ≠ Π groups | ch12 → Ch. 13, 14 | `kappa`, `B`, `A`, `Pi`; `log_law`, `fit_log_law`, `composite_profile` |
| h ⚠️, δ, a, d | full channel height; half-height, pipe radius or layer thickness; pipe radius; diameter | m | $dP/dx=-2\tau_0/h$ (12.90), $-4\tau_0/d$ (12.91) | ch12 | `h`, `delta`, `a`, `d` |
| y₀ (z₀), C_D | roughness length where the rough-wall law (12.93) vanishes; neutral drag coefficient [κ/ln(z/z₀)]² | m; – | not the size of the roughness elements; C_D depends on the reference height | ch12 → **Ch. 13** | `y0`, `z0`; `rough_wall_log_law`, `drag_coefficient_neutral`, `friction_velocity_from_wind` |
| ν_T, κ_T, κ_m, l_T, u_T, A⁺ | eddy viscosity (12.94); eddy diffusivities of heat and of a scalar (12.95)–(12.96); mixing length; turbulent velocity scale; van Driest damping constant | m²/s; m²/s; m²/s; m; m/s; – | ν_T ~ l_Tu_T (12.98); l_T = κy(1 − e^{−y⁺/A⁺}); properties of the flow, not of the fluid | ch12 → **Ch. 13** | `nu_T`, `kappa_T`, `l_T`, `A_plus`; `eddy_viscosity_stress`, `mixing_length_*`, `gradient_diffusion_flux` |
| C_μ, C_ε1, C_ε2, σ_e, σ_ε | constants of the k–ε model (12.103)–(12.105) | – | ν_T = C_μē²/ε̄ (12.104); decay n = 1/(C_ε2 − 1); implied κ = 0.433 | ch12 | `K_EPSILON_CONSTANTS`, `k_epsilon_*` |
| Rf, Rf_cr, Ri ⚠️, Pr_T | flux Richardson number (12.107); its observed critical value ≈ ¼; gradient Richardson number (12.108); turbulent Prandtl number ν_T/κ_T | – | Rf < 0 unstable; Ri = Pr_T·Rf (12.109); Rf_cr an observation, ch11's Ri > ¼ a theorem; ±inf/NaN at zero shear | ch04 Ri → ch11 → ch12 → **Ch. 13** | `flux_richardson`, `gradient_richardson_thermal(…, *, Gamma_a)`, `turbulence_regime`, `turbulent_prandtl` |
| wT, H, α | kinematic heat flux $\overline{wT'}$; surface heat flux ρC_p$\overline{wT'}$; thermal expansion coefficient (1/T for a perfect gas) | K m/s; W/m²; 1/K | **positive upward** | ch12 → **Ch. 13** | `wT`, `H`, `alpha`; `turbulent_heat_flux`, `monin_obukhov_from_fluxes` |
| L_M, ζ, φ_m, ψ_m, β | Monin–Obukhov length (12.110); stability parameter z/L_M; dimensionless shear (κz/u_*)dU/dz; its integrated correction; coefficient 5 of the log-linear profile | m; –; –; –; – | L_M > 0 stable, < 0 unstable, inf neutral; Rf = ζ (12.111); log-linear form valid for 1 + βζ > 0; Businger–Dyer coefficients unread first-hand | ch12 → **Ch. 13** | `L_M`, `zeta`, `beta`; `monin_obukhov_length`, `dimensionless_shear`, `surface_layer_wind`, `surface_layer_state` |
| Γ, Γ_a, Γ_met, Γ_d ⚠️ | Kundu lapse rate dT/dz and its adiabatic value −g/C_p; meteorological −dT/dz and the dry-adiabatic rate +g/C_p | K/m | stable ⇔ dT/dz > Γ_a ⇔ Γ_met < Γ_d; the library string "Γ < Γa" means Γ_met < Γ_d | ch01 → ch11 → ch12 → **Ch. 13** | `Gamma_a` (required), `convention="kundu"\|"met"`; `ch12.adiabatic_lapse_rate()` |
| X_α, r_α(τ), D_T | displacement of a fluid particle in direction α (no sum); Lagrangian velocity autocorrelation; eddy diffusivity ½d$\overline{X_\alpha^2}$/dt (12.127) | m; –; m²/s | Λ_t here is the **Lagrangian** integral scale; D_T grows as $\overline{u_\alpha^2}$t, saturates at $\overline{u_\alpha^2}$Λ_t for t ≫ Λ_t | ch12 → Ch. 13 | `taylor_dispersion*`, `eddy_diffusivity_*`, `langevin_particles`, `dispersion_regime` |
| σ ⚠️ (Gaussian width), R_n, L (step) | standard deviation of a diffusing cloud, σ² = 2νt per coordinate (12.126); end-to-end distance of an n-step walk, (R_n)_rms = L√n (12.125) | m | ch03's core radius used 4νt; ≠ ch11 growth rate, ≠ surface tension | ch12 | `sigma_z`, `diffusivity_from_variance`, `random_walk(n_steps, n_walkers, L)` |
| **ch13 — geophysical fluid dynamics** | | | | | |
| f ⚠️, f₀, θ (latitude) ⚠️ | Coriolis parameter 2Ω sin θ (13.8); its value at the central latitude; latitude | 1/s; 1/s; rad | **signed** (north +, south −); abs(f) for scales, sign(f) for directions; `ValueError` at f = 0 | ch04 → ch13 | `f`, `f0`, `lat_rad` (`coriolis_parameter`, `coriolis_parameter_deg`, `hemisphere`) |
| Ω ⚠️, T_i | Earth's rotation rate; inertial period 2π/abs(f) | rad/s; s | sidereal `OMEGA_EARTH` = 7.292115e-5 (solar-day value 0.27 % smaller); the i of T_i is a label (slip #5) | ch04 → ch13 | `Omega`, `OMEGA_EARTH`, `OMEGA_SOLAR_DAY`, `SIDEREAL_DAY`, `inertial_period` |
| β ⚠️, R | northward gradient of f, 2Ω cos θ₀/R, in $f=f_0+\beta y$ (13.10); Earth's mean radius | 1/(m s); m | β ≥ 0 in both hemispheres; ≠ ch10 scheme β, ch12 log-linear β | ch13 | `beta`, `beta_parameter`, `beta_plane`, `EARTH_RADIUS_MEAN` |
| ν_H, ν_v, F_x, F_y, F_z | horizontal and vertical eddy viscosities; friction force per unit mass (13.6) | m²/s; m/s² | F = (1/ρ)∂τ_ij/∂x_j (slip #2); ν_v ≪ ν_H typically; illustrative values labelled ours | ch13 | `nu_H`, `nu_v`, `eddy_friction`, `anisotropic_eddy_stress` |
| Ro ⚠️, E | Rossby number U/(abs(f)L) (13.13); Ekman number ν/(abs(f)L²) (13.18) | – | ch04's Ro is U/(2ΩL) (factor sin θ); E ≠ energy | ch04 → ch13 | `rossby_number(U, f, L)`, `ekman_number(nu, f, L)` |
| p, ρ (perturbations) ⚠️, ρ₀ | pressure and density departures from the resting state, primes dropped from §13.4 on; reference density | Pa; kg/m³ | geostrophic u = −(1/ρ₀f)∂p/∂y, v = (1/ρ₀f)∂p/∂x | ch13 | `dpdx`, `dpdy`, `rho`, `rho0` |
| V ⚠️, V_g, τ (stress) | complex horizontal velocity u + iv (§13.6–13.7); its geostrophic part; wind stress τ_x + iτ_y at the sea surface | m/s; m/s; Pa | multiplying by i turns 90° to the left; potential flow's complex velocity is u − iv (ch06) | ch13 | `as_complex=True`, `U_g`, `V_g`, `tau_x`, `tau_y` |
| δ ⚠️, D_E | Ekman e-folding thickness √(2ν_v/abs(f)) (13.29); Ekman depth πδ | m | `convention="efold"` (default) or `"pi"`; ≠ boundary-layer δ | ch08 (Stokes layer) → ch13 | `ekman_depth`, `eddy_viscosity_from_depth` |
| M_x, M_y, w_E | Ekman volume transport per unit width (τ_y, −τ_x)/(ρf) (13.30); pumping velocity at the base of the layer | m²/s; m/s | 90° to the right of the stress for f > 0, left for f < 0; w_E > 0 upward (cyclonic stress curl, f > 0) | ch13 | `ekman_transport`, `ekman_pumping`, `ekman_pumping_from_curl`, `sverdrup_transport` |
| η ⚠️, H, h | surface displacement; undisturbed depth; total depth H + η (over an uneven bottom in §13.13) | m | z = 0 at the flat bottom in §13.8; η ≠ Kolmogorov length (ch12) | ch07 → ch13 | `eta`, `H`, `h` |
| c ⚠️, c_n, H_e | long-wave speed √(gH); speed of vertical mode n; equivalent depth c_n²/g (13.62) | m/s; m/s; m | c_n = NH/(nπ) for uniform N (13.71); c also sound speed in §13.2 and complex phase speed in §13.16–13.17 | ch07 → ch13 | `long_wave_speed`, `baroclinic_mode_speed`, `equivalent_depth`, `Modes.c`, `Modes.He` |
| ψ_n(z) ⚠️, w_n, ρ_n, p_n | vertical structure of mode n (13.52), (13.56); modal amplitudes | –; 1/s; kg/m²; m²/s² | ψ_n(0) = 1; orthogonal with **weight 1** for either lid; with a rigid lid index 0 is the first baroclinic mode | ch13 | `Modes.psi`, `w_structure`, `rho_structure`, `modal_amplitudes`, `orthogonality_matrix(kind=)` |
| k, l, m, K ⚠️ | eastward, northward, vertical wavenumbers; horizontal magnitude √(k² + l²) | rad/m | signed; Rossby waves with ω > 0 have k < 0; K also kinetic energy (§13.17) and K₀, K₁, K₂ (§13.18) | ch07 → ch13 | `k`, `l`, `m`, `K` |
| ω (three roots) | frequencies of rotating shallow water, (13.76): two Poincaré roots and one Rossby root | rad/s | fast roots of opposite sign with abs(ω) > abs(f); NaN where the discriminant is negative | ch13 | `shallow_water_omega`, `shallow_water_discriminant`, `poincare_omega`, `rossby_omega`, `kelvin_omega` |
| Λ ⚠️, Λ_E | Rossby radius of deformation c/abs(f); Eady radius NH/abs(f) | m | name the wave speed (external, internal with π, two-layer); Λ_E has **no π**; ≠ ch12 integral scales | ch07 (√(g′H)/f) → ch13 | `rossby_radius`, `rossby_radius_internal(with_pi=)`, `rossby_radius_two_layer` |
| ζ ⚠️, q | relative vorticity ∂v/∂x − ∂u/∂y; shallow-water potential vorticity (ζ + f)/h (13.94) | 1/s; 1/(m s) | q conserved following the motion; linear form ζ − fη/H | ch05 → ch13 | `zeta`, `potential_vorticity`, `sw_potential_vorticity`, `step_vorticity` |
| ψ (stream function) ⚠️ | geostrophic or barotropic stream function, p/(fρ₀) on an f-plane | m²/s | **u = −∂ψ/∂y, v = ∂ψ/∂x** (ch04, ch06, ch11: u = ∂ψ/∂y) | ch13 | `geostrophic_streamfunction`, `barotropic_velocity` |
| θ_K ⚠️, N | angle of the wavevector with the horizontal, tan θ_K = m/k (the book's θ in §13.14); buoyancy frequency | rad; 1/s | internal waves exist for abs(f) < ω < N; N from potential density | ch07 → ch13 | `inertia_gravity_omega`, `inertia_gravity_band`, `N`, `N2` |
| c_x, **c**_g | zonal phase speed ω/k (13.119); group velocity ∇_K ω | m/s | Rossby: c_x < 0 always (no mean flow); for l = 0, c_gx < 0 when abs(k)Λ < 1 and > 0 for shorter waves | ch07 → ch13 | `rossby_phase_speed`, `rossby_group_velocity`, `poincare_group_velocity`, `inertia_gravity_group_velocity` |
| U(y), β − U″ | zonal basic current; northward gradient of absolute vorticity (13.124) | m/s; 1/(m s) | must change sign for barotropic instability (necessary) | ch11 → ch13 | `absolute_vorticity_gradient`, `rayleigh_kuo_criterion`, `rayleigh_kuo_eigs`, `ST.rayleigh_eigs_contour(beta=)` |
| U₀, α ⚠️, αH, c = c_r + ic_i, σ | Eady wind difference between the lids; scaled wavenumber NK/abs(f) (13.139); its product with the lid separation; complex phase speed (13.141); growth rate abs(k)c_i | m/s; 1/m; –; m/s; 1/s | growth for αH < 2.39936; σ_max = 0.30982 abs(f)(dU/dz)/N; α ≠ thermal expansion, ≠ enstrophy flux | ch11 (σ = −ikc) → ch13 | `U0`, `eady_alpha`, `alphaH`, `eady_phase_speed`, `eady_growth_rate`, `eady_max_growth_rate` |
| S(K) ⚠️, K²S(K), α_Z, L_β | energy spectrum of §13.18; enstrophy spectrum; enstrophy flux; Rhines length √(u_rms/β) | m³/s²; m/s²; 1/s³; m | **one-sided, no ½**, mean(u²) = ∫₀^∞ S dK (ch12: two-sided) | ch12 → ch13 | `enstrophy_spectrum`, `two_d_cascade_spectrum(K, K0, eps, alpha_ens)`, `rhines_length`, `barotropic_spectrum` |
| Γ, Γ_a, Γ_met, Γ_d ⚠️ | Kundu lapse rate dT/dz and its adiabatic value; meteorological −dT/dz and the dry-adiabatic rate | K/m | computed in Kundu's sign, the meteorological form always shown beside it; `Gamma_a` required keyword | ch01 → ch12 → ch13 | `lapse_rate_table(dT_dz, *, Gamma_a)` |

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
| ch08 | channel: x along the plates, y across, walls y = 0 (fixed) and y = h (moving); pipe and circular Couette: cylindrical (R, φ, z), capital R; lubrication: x along, y across the gap h(x, t) (Hele-Shaw: z across, (x, y) in the plane); §8.4 y normal to the plate, η = y/√(νt); Ex. 8.6 plane polar (r, θ azimuth); §8.6 spherical (r, θ, φ) with θ from the downstream +x axis, body frame (sphere at rest) or fluid frame | fixed lower wall; pipe axis; slider inlet x = 0 with gap h₀; plate at y = 0; sphere centre | p absolute or gauge (only gradients matter) except lubrication p* = p/P_a (absolute atmospheric); dp/dx book sign (favourable < 0), `G` = −dp/dx alias; signed τ; power into the fluid −2πR₁σ_Rφu_φ; §8.6 p − p∞ | lubrication: L along, h = εL across, U, P_a (or μUL/h²); similarity: √(νt); Stokes layer δ_e = √(2ν/ω); low Re: pressure μU/L; Re per section (pipe diameter, Re_L, sphere diameter 2aU/ν, Re_a radius) | dimensional SI throughout; the scaled lubrication and low-Re equations live in the sympy engines (`lubrication_nondim_sympy`, `low_re_scaling_sympy`); g = G0 in `core.lubrication` and `core.creeping` |
| ch09 | x along the wall from the leading edge (or stagnation point), y normal to it (2-D); jets: x along the jet, y across (free jet symmetric about y = 0, wall jet y = 0 the wall); body angles φ from the **forward** stagnation point; Kármán street rows at y = ±b/2 with spacing a (upper row −Γ in the picture, +Γ in D14 step 1 — an open explainer item); teacup cylindrical (R, z) with z up from the floor | leading edge x = 0 (Blasius δ → 0 there, (9.22)); Thwaites start x[0] with θ₀ (stagnation start: the finite limit); marching inlet x₀ > 0 with a supplied or local Falkner–Skan profile; jet slot x = 0 (virtual origin) | −(1/ρ)dp/dx = U_eU_e′; **dp/dx > 0 adverse**; u = ψ_y, v = −ψ_x; signed τ₀ (> 0 attached, 0 at separation); C_p = (p − p∞)/(½ρU²); C_D of (9.33) for one face; μu_yy(wall) = dp/dx | overall L, U (δ̄ = L Re^{−1/2}); similarity lengths √(νx/U), √(νx/U_e), (Cρν²x²/J)^{1/3}, ∝ x^{3/4}; cylinder/sphere Re on the diameter; St with the cyclic frequency | dimensional SI; similarity profiles and the Thwaites closure non-dimensional; scaled (9.7)–(9.8) in `bl_nondim_sympy`; ρ default 1.2 (air) in BL helpers, 1000 in the teacup; g only in `ball_swing_deflection` (9.81) |
| ch10 | 1-D x ∈ [0, L] (FD, FEM1); 2-D (x, y) with the `[j, i]` layout, y up (cavity lid at y = 1 moving in +x; block and cylinder channels with the stream in +x); **three grids**: node-based (FD, MCK), element meshes (FEM1 [x_{e−1}, x_e]; FEM2 vertices then mid-edge nodes, Fig. 10.17 numbering), staggered C-grid (MAC: p centres, u x-faces, v y-faces); Fourier angle θ = kπΔx (book e^{iπkx_i}) | x = 0 Dirichlet end (g), x = L Neumann end (q); cavity corner (0, 0), walls on node lines (MCK) or cell faces (MAC); block/cylinder centred in the channel (ours: H = 4 block sides, 8 ahead, 20 behind; cylinder W = 5d) | CFL and upwind side with sign(u), \|u\|; truncation error E = "exact minus scheme"/Δt on the left of (10.16); weak continuity with a minus sign (B, Bᵀ symmetric); pressure pinned by its mean (defined up to a constant); ψ = 0 on the walls, ψ_min < 0 clockwise eddy; drag positive downstream | §10.2–10.3 dimensional (u m/s, D m²/s, L m); §10.4 NS non-dimensional (L, U, L/U, ρU², Re (10.81)); cavity p = ρ/Ma², Re = ρ₀UL/μ; block side and cylinder diameter d as lengths; force coefficients per span on ½ρU²d | §10.2–10.3 dimensional SI; §10.4–10.5 solvers non-dimensional (ρ = 1, c = 1/Ma, μ = 1/Re); our run parameters Ma 0.08 (cavity), 0.06 (block), σ = 0.8 (book values private) |
| ch11 | KH: (x, z), z up, interface at z = 0, upper stream 1, lower stream 2 (depth h below in Ex. 11.1); Bénard and double diffusion: z up across the layer, **z ∈ [−½, ½]** in units of d; Taylor: cylindrical (R, φ, z), x = (R − R₁)/d ∈ [0, 1]; stratified shear: (x, z), z up, layer centred at z = 0 or between walls; parallel viscous flows: x streamwise, y across (Poiseuille y ∈ [−1, 1], Blasius y ≥ 0 in δ*, free layers y ∈ (−∞, ∞) truncated at ±y_max); Chebyshev nodes **descending** (`grid.y[0]` is the top/right end); Lorenz phase space (X, Y, Z) | interface z = 0; layer mid-plane z = 0; inner cylinder x = 0; wall y = 0 (Blasius); critical level y_c where U = c_r | normal mode e^{ikx + σt} = e^{ik(x − ct)}, σ = −ikc; Γ = −dT̄/dz in (11.21) (code takes `dT` = T_bottom − T_top); §11.5 Ra signed with dT̄/dz; ψ: §11.7 u = ∂ψ/∂z, §11.8 u = ∂ψ/∂y, §11.14 u = −∂ψ/∂z; Ri sign follows N² at U′ = 0; production P = −∫⟨uv⟩U′dy positive when it feeds the wave | Bénard d, d²/κ; Taylor d, d²/ν; parallel flows L and U₀ per flow (half-width and centreline speed; δ* and U∞; L and U₀ of tanh/sech²); Lorenz time d²/((π² + k²)κ); default g = G0 = 9.80665 | KH and the salt-finger examples dimensional SI; every eigen-solver non-dimensional with the scales stated in its docstring |
| ch12 | signals: time (or space) on the last axis, **ensemble members on axis 0**; shear flows: x streamwise, y across (jet symmetric about y = 0, ξ = y/x; wall at y = 0, channel walls at y = 0 and y = **h = full height**, δ = h/2); surface layer: z up from the ground; isotropic fields: periodic boxes with r the separation vector; dispersion: X_α from the release point | wall y = 0; ground z = 0 with the wind vanishing at z₀; jet origin at the slot (virtual origin dropped unless `x0=`); release at X = 0, t = 0 | Reynolds stress −ρ₀$\overline{u_iu_j}$ (`uv_plus` keys hold −$\overline{uv}$/u_*²); P = deviation from hydrostatic in §12.8–12.10; heat flux positive upward; L_M > 0 stable; Kundu Γ ≡ dT/dz in code with `Gamma_a` required, meteorological Γ_met shown alongside; T̄, T′ potential temperature; two-sided spectra with ∫S = variance | outer ΔU, L (Re_L); Kolmogorov η, u_K, τ_η; wall units u_*, ν/u_* (Re_τ = δ⁺); surface layer u_*, L_M; dispersion u_rms, Λ_t (Lagrangian); default g = G0 = 9.80665; κ = 0.41 in §12.9–12.10, 0.4 in the surface-layer block | dimensional SI for statistics, jets, the surface layer and dispersion; wall functions in wall units with dimensional wrappers (`wall_units`, `from_wall_units`, `composite_profile`); channel budgets in u_*⁴/ν; empirical constants are required keywords |
| ch13 | local tangent plane: **x east, y north, z up**, (u, v, w); latitude θ in radians; f-plane or β-plane (y measured north from the central latitude); grids `[j, i]` = (y, x) on a C-grid (η centres, u west faces, v south faces, vorticity at corners); vertical-mode nodes from z[0] = −H to z[-1] = 0; Eady: z ∈ [0, H] between two lids | **z = 0 changes by section** (trap T4): sea surface (§13.6, §13.9, §13.14), solid surface (§13.7), flat bottom (§13.8), lower lid (§13.17); coast at y = 0 with the sea in y ≥ 0 (Kelvin wave); step at x = 0 (adjustment, flow over a step) | p, ρ perturbations (primes dropped); **f signed**; Coriolis term −fv, +fu on the left; wind stress at the surface τ = ρν_v ∂u/∂z; **ψ with u = −∂ψ/∂y**; eddy fluxes of the Eady wave positive northward (y) and upward (z); lapse rate Kundu dT/dz with the meteorological form shown | f-plane 1/abs(f) and Λ = c/abs(f); Ekman δ; Eady Λ_E = NH/abs(f) and NH/(abs(f)U₀); numbers from `ch13.illustrative_inputs()` (35° N/S, 60° N for Ekman layers, 12° N for Rossby waves) | dimensional SI throughout; Eady functions also in αH and σNH/(fU₀); model runs dimensional; the Rayleigh–Kuo solver non-dimensional (jet width and speed) |
