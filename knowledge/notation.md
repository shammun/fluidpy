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
| ω ⚠️ | vector of R = vorticity ∇×u | 1/s | ω_k = −½ ε_ijk R_ij; **ch07 uses ω for angular frequency** | ch02 → | `omega`, `vector_from_antisymmetric(rotation_tensor(G))` |
| ½ω₃ = A₂₁ | spin rate of a fluid element in a plane flow (E3 readout) | 1/s | never written "ω" | ch02 | `vector_from_antisymmetric(antisymmetric_part(G))` |
| Γ ⚠️ | shear/strain rate of a linear flow (u₁ = Γx₂; Ex. 2.4 S₁₂) | 1/s | **ch01 Γ = lapse rate; ch05 Γ = circulation** | ch02 | `Gamma` |
| Γ_circ ⚠️ | circulation ∮u·t ds | m²/s | + counterclockwise about n; the notebook writes Γ_circ to avoid the clash | ch02 → (ch05 calls it Γ) | `circulation` |
| K ⚠️ | strength of the irrotational vortex u_θ = K/r | m²/s | Γ_circ = 2πK round the core; not Taylor's K | ch02 | `K` in `irrotational_vortex_field` |
| m ⚠️ | 2-D point-source strength (outflux per unit depth) | m²/s | not a molecule mass | ch02 | `m` in `point_source_field` |
| V, A ⚠️ | volume of a region; area of an open surface (or of a loop) | m³; m² | A here is an area, not a tensor or matrix | ch02 → | `volume`, `Loop.area` |
| n_c | in-surface normal to the rim C, pointing **into A** | – | t = n_c × n | ch02 → | `boundary_tangent(n_c, n)` |
| t ⚠️ | unit tangent of a curve (not time here) | – | counterclockwise about n | ch02 → | `Loop.tangents` |
| s ⚠️, ds | arc length along a curve (not entropy here) | m | | ch02 → | `Loop.ds` |
| u_i,j | comma notation for ∂u_i/∂x_j (2.36) | 1/s | a comma index transforms as a vector index | ch02 → | `comma_to_partial` |
| singular_at | point where a test field is not differentiable | m | Stokes/Gauss report `hypothesis_ok = False` | ch02 → | `VectorField.singular_at` |

## Coordinate and sign conventions per chapter
| Chapter | Axes (which is "up") | Origin / reference level | Stress / pressure sign | Reference scales (L, U, T) | Dimensional or non-dimensional code |
|---|---|---|---|---|---|
| ch01 | z up (§1.7, §1.10); Couette y from fixed wall (0) to moving plate (h) | p0 at z = 0; θ reference p_ref = 1000 hPa; parcel rest height z_o | pressure absolute and isotropic, acts along the inward normal; τ_xy = +μ ∂u/∂y signed; lapse rate Kundu dT/dz; first law q in / w on | none fixed (Π groups carry their own; E2 uses the clock h²/ν) | dimensional SI throughout; only Π groups are dimensionless |
| ch02 | right-handed x₁x₂x₃ (no preferred "up"); rotated frame shares the origin; angles counterclockwise about e₃; grid arrays `[k, j, i]` = (z, y, x), components and `h` in (x, y, z) | origin of both frames; boxes/loops centred at x0 | τ_ij tensile positive, +e_i face → +e_j; traction f = n·τ (first index); ∇·τ on the second index; pressure τ = −pδ; passive C (x' = Cᵀx); book A:B = A_ij B_ji; R = G − Gᵀ ↔ ω = ∇×u, A = ½R ↔ ½∇×u; Stokes n_c into A, t counterclockwise about n; outward n on closed surfaces | none (pure mathematics) | dimensional where physical (Pa, 1/s, m); most results unit-agnostic |
