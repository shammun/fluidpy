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

## Coordinate and sign conventions per chapter
| Chapter | Axes (which is "up") | Origin / reference level | Stress / pressure sign | Reference scales (L, U, T) | Dimensional or non-dimensional code |
|---|---|---|---|---|---|
| ch01 | z up (§1.7, §1.10); Couette y from fixed wall (0) to moving plate (h) | p0 at z = 0; θ reference p_ref = 1000 hPa; parcel rest height z_o | pressure absolute and isotropic, acts along the inward normal; τ_xy = +μ ∂u/∂y signed; lapse rate Kundu dT/dz; first law q in / w on | none fixed (Π groups carry their own; E2 uses the clock h²/ν) | dimensional SI throughout; only Π groups are dimensionless |
| ch02 | right-handed x₁x₂x₃ (no preferred "up"); rotated frame shares the origin; angles counterclockwise about e₃; grid arrays `[k, j, i]` = (z, y, x), components and `h` in (x, y, z) | origin of both frames; boxes/loops centred at x0 | τ_ij tensile positive, +e_i face → +e_j; traction f = n·τ (first index); ∇·τ on the second index; pressure τ = −pδ; passive C (x' = Cᵀx); book A:B = A_ij B_ji; R = G − Gᵀ ↔ ω = ∇×u, A = ½R ↔ ½∇×u; Stokes n_c into A, t counterclockwise about n; outward n on closed surfaces | none (pure mathematics) | dimensional where physical (Pa, 1/s, m); most results unit-agnostic |
| ch03 | right-handed Cartesian (no preferred "up"); plane polar (r, θ from +x), cylindrical (R, φ, z), spherical (r, θ from +z, φ); angles and rotation counterclockwise positive (shear spins clockwise: ω₃ = −γ); field callables `u(x, t)` with coordinates on axis 0; Galilean frame O′ at constant U (x = x′ + Ut + x′_o, u′ = u − U); rotating frame u = Ω × x + u′ at the coinciding instant | cylinder centre at the origin at t = 0 for every observer (E3); Ex. 3.1 port at the origin; vortices centred at the origin; RTT shapes with explicit reference geometry (E7 interval [1, 3] + ȧt, ḃt; ellipse a = 2 + ȧt, b = 1 + ḃt) | R = G − Gᵀ (no ½), ω = ∇×u, spin ½ω; γ = 2S₁₂; RTT outward n, signed b·n; Leibniz lower term subtracted | none (E6 draws in r/σ; its real-vortex modes use metres) | dimensional SI throughout |
| ch04 | right-handed Cartesian, **z up**, g = −g e_z, Φ = gz (4.18); cylindrical (R, φ, z) for rotating flows and Ex. 4.5; noninertial frame O′ translating at U(t) and rotating at Ω(t) with basis e′_i; NH Ω_z > 0; latitude φ | CVs with explicit geometry (E1 boxes around the wake, bore, jet, rocket, balloon); hydrostatic base state p_s(z), ρ_s(z) for Boussinesq; free surface η = 0 | τ = −pδ + σ, tensile positive; traction f_j = n_iτ_ij and Cauchy's divergence on the **first** index; outward n, signed (u − b)·n; drag on the body +x, on the fluid −F_D; acceleration terms +2Ω × u′, +Ω × (Ω × x′) (4.43) vs forces −2Ω × u′, −Ω × (Ω × x′) (4.45); Stokes assumption μ_v = 0; 2-D ψ: u = ∂ψ/∂y; primes: rotating frame (§4.7), perturbation (§4.9), dummy (4.67) | `core.similarity.Scales` holds one set per use with its time scale (1/Ω (4.100) or l/U (4.109)) and pressure scale (ρU², μU/l or ρgl); Ro = U/(2Ωl) forward pointer | dimensional SI in functions; non-dimensional via `Scales` and the `nondimensional_*` coefficient routines; g default 9.81 in `ch04`, 9.80665 in `core` |
