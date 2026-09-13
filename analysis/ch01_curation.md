# Chapter 1 — Introduction: curation
(from `analysis/ch01.md` — 106 inventory rows, 40 numbered equations, 33 derivations; `book.yaml` ch01; knowledge files
are empty templates because this is the first chapter, so there is nothing to RECAP from earlier chapters and no primer
to reuse: every prerequisite below needs a primer. 2026-09-12, concept-curator.)

**Decisions in one paragraph.** Every NEW row is CORE unless it is a figure, a value, a remark, a special case or an
alternative form of a CORE item of this chapter that can be shown inside its parent block, in which case it is a NOTE
with the parent named. The analyst's prose-SKIPs were re-tiered: a new definition is a new idea (unclear → CORE). Only
the exercises and the bibliography are SKIP. The three pre-book rows of §1.2 are RECAP (source: pre-book school physics).
CORE IDs follow the teaching order. Four derivations the book never writes but the lesson needs were **added** (D34
kinetic pressure, D35 continuum noise law, D36 N² from θ, D37 Archimedes from hydrostatics). With them D18 (parcel
buoyancy) no longer relies on a pre-book fact, and the notebook shows why sign N² = sign dθ/dz. The analyst's
**suspected typos** (§9 of the analysis) are taught as written there: Clausius–Duhem with the actual δq, and (1.35) as
"has the same sign as".

**Totals:** CORE 76 · RECAP 3 · NOTE 25 · SKIP 2 (= 106 rows) · derivations 37 (★ 19 · ★★ 15 · ★★★ 3) · explainers 5
+ 1 backup.

## 1. Teaching order (CORE IDs grouped by book section; one sentence each: "once you see X, Y follows")

**§1.1 Fluid Mechanics**
- C01 What fluid mechanics is and its three routes (analysis, computation, experiment): once you see the map of the subject, you can see where every later idea in this chapter will be used.

**§1.2 Units of Measurement** (no new idea: R01–R03 recap SI units, prefixes and kelvin with `pint`, so every later number carries its unit)

**§1.3 Solids, Liquids, and Gases**
- C02 Fluid vs solid: once you see that a fluid keeps deforming under *any* shear stress while a solid settles at a fixed strain, the words "fluid" and "flow" have an exact meaning.
- C03 Normal stress: compression vs tension: once shear is understood, the other half of stress follows: fluids resist being squeezed but liquids break (cavitate) when pulled.
- C04 Liquid vs gas (molecular spacing, free surface): once you see that molecules in a gas are about 10× farther apart than in a liquid, it follows that gases fill their container and liquids form a free surface.

**§1.4 Continuum Hypothesis**
- C05 Pressure as the average of molecular impacts: once you see pressure as momentum delivered by many molecules, you can ask how many molecules an average needs, which is the continuum question.
- C06 The continuum hypothesis (density at a point as an average over a window of box sizes): once averages settle to a plateau as the box grows, the idea of a "value at a point" follows.
- C07 Knudsen number Kn = l/L: once the continuum needs boxes much larger than molecular scales, the single ratio l/L follows as the test of whether continuum mechanics applies.

**§1.5 Molecular Transport Phenomena**
- C08 Transport by random molecular motion (diffusion down gradients): once molecules carry their properties as they wander, it follows that any uneven profile gets smoothed out.
- C09 Mass fraction Y and partial density ρY: once we need to follow one species in a mixture, the quantity that diffuses has to be named.
- C10 Fick's law (1.1): once Y diffuses, the flux points down its gradient with strength ρκ_m, which is the law.
- C11 Fourier's law (1.2): the same picture applied to molecular energy gives heat flux = −k∇T.
- C12 Newton's law of viscosity (1.3): the same picture applied to momentum gives shear stress = μ du/dy. A stress is a momentum flux.
- C13 How μ depends on temperature (gases up, liquids down): once viscosity comes from molecules carrying momentum (gases) or from sticking to neighbours (liquids), the opposite temperature trends follow.
- C14 Kinematic viscosity ν = μ/ρ (1.4): once momentum diffuses, its diffusivity has units m²/s like κ_m, and ν rather than μ sets how fast motion spreads.

**§1.6 Surface Tension**
- C15 Surface tension σ (force per length or energy per area): once an interface is pictured as a stretched membrane, curved interfaces must carry a pressure difference.
- C16 Laplace pressure jump (1.5): once tension pulls along a curved patch, the balance gives Δp = σ(1/R₁ + 1/R₂); a sphere gives 2σ/R.

**§1.7 Fluid Statics**
- C17 Static, absolute and gauge pressure (atm, bar): once pressure is a number at a point, we have to say what it is measured from.
- C18 Pressure at a point is the same in every direction (1.6): once a shrinking wedge's weight vanishes faster than its face forces, pressure can only be a scalar.
- C19 Pascal's law, no horizontal gradient at rest (1.7): once pressure is a scalar and gravity is vertical, the horizontal forces on a cube can only cancel.
- C20 Hydrostatic law dp/dz = −ρg (1.8): the vertical balance on the same cube gives the law that everything in §1.10 builds on.
- C21 Uniform density: p = p₀ − ρgz (1.9), with buoyancy: integrate (1.8) once and you get pressure growing linearly with depth, and the net pressure force on a body turns out to be Archimedes' buoyancy.
- C22 Capillary rise (Ex. 1.1): once you have Laplace (1.5) and hydrostatics (1.9), the height a meniscus lifts water is a two-line balance.

**§1.8 Classical Thermodynamics**
- C23 Thermodynamic system, equilibrium and relaxation time: once molecules collide ~10¹⁰ times a second, a small parcel settles to equilibrium almost instantly compared with flow times.
- C24 Fluid particle: once a parcel is big enough to average (C06) and relaxes fast (C23), we can give it p, ρ, T and follow it.
- C25 First law δq + δw = Δe (1.10): once a fluid particle has an energy, heat in plus work done on it must change that energy.
- C26 Path functions (q, w) vs state functions (e); reversible processes: once you draw two paths between the same states, q and w differ but Δe does not.
- C27 Reversible first law de = dq − p dv (1.11): once work is a slow push on a boundary, it is −p dv per unit mass.
- C28 Equations of state (1.12): once two properties fix the state, p = p(v, T) and e = e(p, T) are surfaces.
- C29 Enthalpy h = e + pv (1.13): once constant-pressure processes are common, grouping e with pv gives the natural energy.
- C30 C_p = (∂h/∂T)_p (1.14): once h is a state function, its slope with T at fixed p is a property.
- C31 C_v = (∂e/∂T)_v (1.15): the same idea for e at fixed v.
- C32 Heat per degree along isochoric and isobaric paths: once C_v and C_p are defined as slopes of state functions, you can show they equal heat per degree only along those two reversible paths.
- C33 Entropy (1.16): once some integral along reversible paths turns out not to depend on the path, it defines a new state function s.
- C34 Clausius–Duhem inequality: once entropy exists, irreversible processes (stirring, free expansion) can only raise it beyond ∫δq/T. This is why μ and k must be positive.
- C35 Gibbs relations (1.18): once T ds = dq and de = dq − p dv hold, eliminating dq gives relations between state functions that hold for *any* process.
- C36 Speed of sound c² = (∂p/∂ρ)_s (1.19): once you have an equation of state, how stiff it is at constant entropy sets the speed of pressure signals, and infinite stiffness is the incompressible limit.
- C37 Thermal expansion coefficient α (1.20): the same kind of partial derivative, taken at constant pressure, measures how density changes when a fluid is heated.

**§1.9 Perfect Gas**
- C38 Molecular gas law pV = nk_BT (1.21): once pressure is molecular impacts (C05) and temperature is molecular kinetic energy, this law follows for non-interacting molecules.
- C39 Universal and specific gas constants (k_B, Avogadro, R_u, M_w, R): once we count molecules in kilomoles, the constants turn molecule counts into masses.
- C40 Continuum perfect-gas law p = ρRT (1.22): once n/V × m is the continuum density, (1.21) becomes the law used in the rest of the book.
- C41 A perfect gas has e = e(T), h = h(T): once Gibbs (C35) is combined with p = ρRT, internal energy cannot depend on volume.
- C42 R = C_p − C_v (1.23): once e and h depend on T only, h = e + RT differentiates to this.
- C43 Ratio of specific heats γ (1.24): once C_p and C_v are both known, their ratio is the one number that sets adiabatic behaviour.
- C44 Adiabatic vs isentropic processes: once you have entropy, "no heat" and "no heat and no friction" turn out to be different claims.
- C45 Isentropic law p/ρ^γ = const (1.25): once ds = 0 in the Gibbs relations, separating variables gives this power law.
- C46 Isentropic T and ρ ratios (1.26): once p/ρ^γ is fixed, p = ρRT gives how T and ρ follow p.
- C47 c = √(γRT) (1.27): differentiating p = Kρ^γ in (1.19) gives the sound speed of air from its temperature.
- C48 α = 1/T for a perfect gas (1.28): once ρ = p/RT, heating at fixed p gives this simple α, used later in the Boussinesq approximation.

**§1.10 Stability of Stratified Fluid Media**
- C49 In a static fluid, hydrostatics plus an equation of state tie p, ρ and T together: once one profile is known the other two follow, as in the standard atmosphere built from T(z).
- C50 The displaced-parcel equation: once a parcel keeps its own (adiabatic) density while its surroundings change, Newton's law with buoyancy gives ζ'' + N²ζ = 0.
- C51 Brunt–Väisälä frequency N² (1.29): once the parcel equation is linear, its coefficient is the frequency of vertical oscillation.
- C52 Stable / neutral / unstable from the sign of N²: once N² can be negative, cos becomes cosh and the parcel runs away.
- C53 Lapse rate Γ ≡ dT/dz: once temperature profiles matter, their slope needs a name, and we follow the book, so Γ is negative where temperature falls with height (the standard troposphere has Γ = −6.5 K/km).
- C54 Adiabatic lapse rate Γ_a = −gαT/C_p (1.30): once a parcel rises isentropically through hydrostatic surroundings, its temperature changes at this rate, Γ_a = −g/C_p ≈ −9.8 K/km for air (negative: the parcel cools as it rises).
- C55 Potential temperature θ (1.31): once adiabatic cooling is predictable, bringing a parcel to a reference pressure removes it and leaves θ.
- C56 dθ/dz and lapse rates (1.32): once you differentiate the log of (1.31), dθ/dz is proportional to Γ − Γ_a.
- C57 Stability from the sign of dθ/dz (and N² = (g/θ)dθ/dz): once (1.32) holds, a gas column is stable exactly where θ increases upward.
- C58 Potential density ρ_θ (1.33): the same reference-pressure trick applied to density.
- C59 ρ_θ and θ carry the same information (1.34): once θρ_θ is constant, their log-gradients are equal and opposite.
- C60 Seawater: salinity S and ρ(T, p, S): once the fluid is the ocean, density needs a third variable, and parcels keep their S.
- C61 Ocean stability criterion (1.35): once a parcel's density changes only through compression (−ρg/c²), stability is the sign of dρ/dz + ρg/c².
- C62 Isothermal atmosphere p = p₀e^{−gz/RT}: once T is constant, (1.8) with p = ρRT is separable and pressure decays exponentially.
- C63 Scale height H = RT/g: once the decay is exponential, its e-folding height is the natural vertical length of the atmosphere (≈ 8 km).

**§1.11 Dimensional Analysis**
- C64 Dimensional homogeneity: once laws cannot depend on the units we choose, every term in a correct equation must have the same dimensions.
- C65 Step 1, choose the variables (pipe pressure drop, (1.38)): once the answer can only depend on the variables listed, choosing that list is the physics step.
- C66 Step 2, dimensions [q] and base dimensions M, L, T, θ: once quantities are written in base dimensions, each one becomes a vector of exponents.
- C67 Step 2, the dimensional matrix (1.39): stacking those vectors as columns gives one matrix that holds all the information.
- C68 Step 3, rank from minors: once there is a matrix, its rank r counts the independent dimensions really present.
- C69 Buckingham's Π theorem (1.37), with n − r: once dimensionless groups are exactly the null space of the matrix, rank–nullity gives n − r of them. (We teach it after C67–C68, not in book order, because its proof needs the matrix and its rank.)
- C70 Step 5, exponent algebra (with inspection as a shortcut): once you pick r repeating variables with a nonzero determinant, each group's exponents come from an r × r linear solve.
- C71 Combining groups (products and powers are groups too, only n − r are independent): once groups are vectors, a different basis gives the same law.
- C72 Steps 6–7, the dimensionless pipe law (1.40) and using physics to simplify: once the groups are known, the unknown function has 3 arguments instead of 6, and Δp ∝ Δx removes one more.
- C73 Ex. 1.2, scale height by dimensional analysis: once a problem has a single group, that group must be a constant, and H ∝ R_uT/gM_w follows without solving anything.
- C74 Ex. 1.3, Pythagoras by dimensional analysis: once area = C²φ(β) and similar triangles share φ, A² + B² = C² follows.
- C75 Ex. 1.4, G. I. Taylor's blast energy E = Kρ D⁵/t²: once only E, ρ, D, t matter, a photo sequence gives the yield (D ∝ t^{2/5}).
- C76 Ex. 1.5, Rayleigh scattering ∝ λ⁻⁴: once physics fixes the powers of d and V, the wavelength dependence follows, and with it the blue sky.

## 2. Tiers
| ID | Item | § | Tier | Why this tier | Treatment (parent for NOTE, pointer text for SKIP, source chapter for RECAP) |
|---|---|---|---|---|---|
| C01 | What fluid mechanics is; the three routes (analysis, computation, experiment) [#1] | 1.1 | CORE | new in this chapter: the subject's definition and method (the Prandtl/boundary-layer anecdote is history and gets one pointer line inside) | plain words + figure: the book's chapter map drawn from `book.yaml` with the route to Ch. 13 (GFD) highlighted; code builds the graph |
| R01 | SI base quantities and derived units N, Pa, J, W, Hz [#2] | 1.2 | RECAP | SEEN (pre-book school physics) | source: pre-book; recap with `core.units` (`Q_`, `.to()`), unit table printed by code |
| R02 | SI prefixes [#3] | 1.2 | RECAP | SEEN (pre-book) | source: pre-book; one line of `pint` prefix conversions (hPa, mN/m, nm) used later in the chapter |
| R03 | Celsius–kelvin relation [#4] | 1.2 | RECAP | SEEN (pre-book) | source: pre-book; `core.units.celsius_to_kelvin`, warning on pint offset units (`delta_degC`) |
| C02 | Fluid vs solid: a fluid deforms without limit under any shear stress [#5] | 1.3 | CORE | new definition | primer on stress/strain → γ(t) for solid (τ/G) and Newtonian fluid (τt/μ) → from-scratch loop vs `ch01.shear_deformation_history` → animation of the two elements (Fig. 1.1 idea) |
| N01 | Plastic (yield) and viscoelastic (partial memory) materials [#6] | 1.3 | NOTE | intermediate cases of the C02 deformation picture | parent C02: overlay a Bingham (no flow below τ_y) and a Maxwell (τ/G + τt/μ) curve on the C02 γ(t) figure, one code line each |
| C03 | Normal stress: compression and tension; liquids cavitate under tension [#7] | 1.3 | CORE | new idea (the other half of stress) | traction split into normal and shear parts (dot product) → arrow figure on a surface element; tiny example; code `ch01.traction_components` |
| C04 | Liquid vs gas; free surface [#8] | 1.3 | CORE | new definition | number density from ρ and molecular mass → mean spacing n^{-1/3} (water 0.31 nm, air 3.4 nm) → side-by-side molecule scatter at true relative spacing; code `ch01.mean_molecular_spacing` |
| N02 | Fig. 1.1: solid element deflects to a fixed shape, fluid keeps deforming [#9] | 1.3 | NOTE | figure of C02 | parent C02: it is the C02 animation (our drawing), with load removed half-way |
| C05 | Pressure as the statistical average of molecular collision force per area [#10] | 1.4 | CORE | new definition | derivation D34 (kinetic pressure) → Maxwellian velocity histogram → from-scratch wall-momentum sum vs `ch01.molecular_pressure` → figure: pressure estimate converging as molecule count grows |
| C06 | Continuum hypothesis: properties at a point as averages over δV [#11] | 1.4 | CORE | new, load-bearing | D35 noise ∝ N^{-1/2} → from-scratch particle count in boxes vs `ch01.sample_density` → log–log figure (noise → plateau → macroscopic variation) + zoom animation + explainer `continuum_averaging_volume` |
| C07 | Knudsen number Kn = l/L [#12] | 1.4 | CORE | new dimensionless group | tiny example (a 1 µm droplet in air) → `ch01.knudsen_number`, `mean_free_path_jennings` → plotly slider (altitude → l) of Kn vs body size with regime bands; explainer `continuum_averaging_volume` |
| N03 | Mean free path of air at room conditions; water's is much smaller [#13] | 1.4 | NOTE | value used in C07 | parent C07: computed l ≈ 67 nm (Jennings) marked on the Kn figure; order-of-magnitude label |
| C08 | Random molecular motion diffuses species, heat and momentum down gradients [#14] | 1.5 | CORE | new physical picture | random walkers spreading (figure) → model PDE ∂f/∂t = D∂²f/∂y² (primer: finite differences) → from-scratch FTCS loop vs `core.diffusion.ftcs_diffusion_1d` → profile relaxation animation; explainer `viscosity_momentum_diffusion` |
| C09 | Mass fraction Y and partial density ρY [#15] | 1.5 | CORE | new quantity | tiny example: dry air N₂/O₂/Ar mole → mass fractions → stacked bar figure (mole vs mass fraction); code `ch01.mass_fractions` |
| C10 | Fick's law J_m = −ρκ_m∇Y (1.1) [#16] | 1.5 | CORE | numbered law | gradient primer → `ch01.fick_mass_flux` → profile Y(y) with flux arrows; mode "species" of explainer `viscosity_momentum_diffusion` |
| N04 | Fig. 1.2: Y increases with y, flux across AB points down [#17] | 1.5 | NOTE | figure of C10 | parent C10: the C10 arrow figure (our drawing; caption typo "Y" → "y" mentioned) |
| C11 | Fourier's law q = −k∇T (1.2) [#18] | 1.5 | CORE | numbered law | `ch01.fourier_heat_flux` → heat flux through a window pane (tiny example) → T profile + flux arrow figure; mode "heat" of explainer `viscosity_momentum_diffusion` |
| C12 | Newton's law of friction τ = μ du/dy (1.3); dynamic viscosity μ [#19] | 1.5 | CORE | numbered law, load-bearing | shear stress as momentum flux → tiny example (Couette gap) → from-scratch central differences vs `ch01.shear_stress_profile` → u(y) and τ(y) figure; explainer `viscosity_momentum_diffusion` |
| N05 | Fig. 1.3: τ on AB, u(y) relaxes toward the dashed profile [#20] | 1.5 | NOTE | figure of C12 | parent C12: momentum-relaxation animation from C08's FTCS with D = ν, τ(y) panel |
| C13 | μ of a gas rises (≈ T^{1/2}), μ of a liquid falls with T [#21] | 1.5 | CORE | new physical behaviour | `ch01.viscosity_power_law`, `sutherland_viscosity`, `water_viscosity` → plotly slider on the exponent n vs Sutherland (air) + water curve on a second panel |
| C14 | Kinematic viscosity ν = μ/ρ (1.4) [#22] | 1.5 | CORE | numbered definition | tiny example: air ν ≈ 15× water ν although μ_air ≈ μ_water/55 → `ch01.kinematic_viscosity` → bar figure μ vs ν (log) for air, water, oil, honey + diffusion time h²/ν; explainer `viscosity_momentum_diffusion` |
| N06 | The transport laws are linear and contain first derivatives only [#23] | 1.5 | NOTE | property shared by C10–C12 | parent C14: a three-row table (flux, gradient, diffusivity κ_m, k/ρC_p, ν in m²/s) with one line showing doubling the gradient doubles the flux |
| C15 | Surface tension σ: force per length = energy per area [#24] | 1.6 | CORE | new quantity | wire-frame film tiny example (F = 2σL) → `ch01.surface_tension_water` (IAPWS) → σ(T) figure with 20 °C marker |
| N07 | Spherical interface: σ(2πR) = Δp πR², Δp = 2σ/R [#25] | 1.6 | NOTE | special case R₁ = R₂ of (1.5) | parent C16: derived first (D01), then `laplace_pressure_jump(sigma, R)` line; Δp vs R curve on the C16 figure |
| C16 | Laplace pressure jump Δp = σ(1/R₁ + 1/R₂) (1.5) [#26] | 1.6 | CORE | numbered law | D01 then D02 → tiny example (1 mm droplet) → `ch01.laplace_pressure_jump` → Δp vs R log–log for droplets/bubbles + 3-D plotly surface (N09) |
| N08 | Capillary tube, capillarity (Ch. 4 treats it as a boundary condition) [#27] | 1.6 | NOTE | vocabulary introduced with C15 | parent C15: two sentences and a forward pointer to C22 (Ex. 1.1) and Ch. 4 |
| N09 | Fig. 1.4: hemispherical drop forces; patch with two principal radii [#28] | 1.6 | NOTE | figure of C16 | parent C16: plotly 3-D sphere cap and saddle with principal circles (our drawing) |
| C17 | Static pressure; absolute vs gauge; standard atmosphere; bar [#29] | 1.7 | CORE | new definitions | tiny example (tyre gauge 2 bar → absolute) → `core.statics.gauge_pressure` → number-line figure (vacuum, atm, gauge readings) |
| C18 | Pressure at a point is the same in all directions (1.6) [#30] | 1.7 | CORE | numbered result | D03 → `ch01.wedge_pressure_difference` → plotly slider shrinking the wedge (dz) with p₂ − p₁ → 0 and the force triangle closing |
| N10 | Fig. 1.5: triangular element and pressure forces [#31] | 1.7 | NOTE | figure of C18 | parent C18: the slider figure's left panel (our drawing) |
| C19 | Pascal's law ∂p/∂x = ∂p/∂y = 0 at rest (1.7) [#32] | 1.7 | CORE | numbered law | D04 → `np.gradient` of a 3-D hydrostatic field from `hydrostatic_pressure_uniform` shows zero x, y components → heatmap of p on a vertical slice with horizontal isobars |
| C20 | Hydrostatic law dp/dz = −ρg (1.8) [#33] | 1.7 | CORE | numbered law, load-bearing | D05 → from-scratch Euler integration vs `core.statics.integrate_hydrostatic` (ρ(z) given) → p(z) figure for a layered tank; cube diagram (N11) |
| C21 | Uniform density p = p₀ − ρgz (1.9); pressure rise ρgh; buoyancy from pressure [#34] | 1.7 | CORE | numbered result | D06 then D37 (Archimedes, added) → tiny example 10 m of water ≈ 1 atm → `core.statics.hydrostatic_pressure_uniform`, `buoyancy_force` → plotly slider on the upper-layer density of a two-layer tank with a block showing face forces |
| N11 | Fig. 1.6: fluid cube, top/bottom pressure difference balances weight [#35] | 1.7 | NOTE | figure of C20 | parent C20: cube with p and p + dp arrows (our drawing) |
| C22 | Capillary rise h = 2σ sin α/(ρgR) (Ex. 1.1) [#36] | 1.7 | CORE | worked result | D07 → tiny example (1 mm tube, water) → `ch01.capillary_rise` (+ contact-angle conversion warning) → plotly slider on α: h vs R, and p(z) along E–F |
| N12 | Fig. 1.7: meniscus, pressure along E–F, force balance on ABCD [#37] | 1.7 | NOTE | figure of C22 | parent C22: p(z) panel of the C22 slider figure (our drawing) |
| C23 | Thermodynamic system, equilibrium, relaxation time [#38] | 1.8 | CORE | new definitions | tiny example: collision time l/c̄ ≈ 10⁻¹⁰ s vs flow time L/U → `ch01.mean_molecular_speed`, `collision_time` → log time-scale bar figure + exponential approach to equilibrium |
| C24 | Fluid particle: fixed molecules, big enough to average, relaxes fast [#39] | 1.8 | CORE | new, load-bearing (Ch. 3 material volume) | combines C06 and C23: a 10 µm air parcel has ~2.5×10¹⁰ molecules, noise ~10⁻⁵ → figure: size–time "fluid-particle window" with the C06 noise curve and C23 time scale; explainer `continuum_averaging_volume` |
| C25 | First law δq + δw = Δe (1.10) [#40] | 1.8 | CORE | numbered law | sign convention (work on the system positive) → tiny example (stirred insulated cup) → `core.thermo.process_heat_work` → energy-bar figure (q, w, Δe); explainer `heat_work_paths` |
| C26 | Heat and work are path functions, e is a state function; reversible process [#41] | 1.8 | CORE | new, classic confusion | two paths on a p–v diagram between the same states → from-scratch trapezoid −∫p dv vs `process_heat_work` → figure with shaded work areas + piston animation; explainer `heat_work_paths` |
| C27 | Reversible first law de = dq − p dv (1.11) [#43] | 1.8 | CORE | numbered law | D08 → isothermal path tiny example (w = −RT ln(v₂/v₁)) → `process_heat_work` → p–v figure with −p dv strip highlighted; explainer `heat_work_paths` |
| N13 | Specific volume v = 1/ρ [#42] | 1.8 | NOTE | alternative form of density (C06) | parent C27: one line `specific_volume(1.225)`; v axis of the C27 figure |
| C28 | Equations of state p = p(v, T), e = e(p, T) (1.12); two properties fix the state [#44] | 1.8 | CORE | numbered statement | `core.thermo.perfect_gas_state` (give any two of p, ρ, T) → plotly 3-D p(v, T) surface with a state dot; seawater needs S (pointer to C60) |
| C29 | Enthalpy h = e + pv (1.13) [#45] | 1.8 | CORE | numbered definition | tiny example for air at 300 K → `core.thermo.enthalpy` → stacked bars e and pv vs T |
| C30 | C_p = (∂h/∂T)_p (1.14) [#46] | 1.8 | CORE | numbered definition | primer "partial derivative with a variable held fixed" → from-scratch central difference vs `core.thermo.specific_heat_cp` → h(T) curves at two pressures with tangent slopes |
| C31 | C_v = (∂e/∂T)_v (1.15) [#47] | 1.8 | CORE | numbered definition | `core.thermo.specific_heat_cv` on e(T, v) → e(T) at two volumes with tangent; C_v vs C_p on the same axes |
| C32 | Heat per degree = C_v (constant v) or C_p (constant p) for reversible p dv work; stirring shows why this is not the definition [#48] | 1.8 | CORE | new unnumbered result | D09 → isochoric and isobaric heating via `process_heat_work` → q vs T figure with slopes C_v and C_p, plus a stirred path breaking q = C_vΔT |
| C33 | Second law (i): entropy s₂ − s₁ = ∫dq_rev/T (1.16) [#49] | 1.8 | CORE | numbered theorem | two reversible paths give the same ∫dq/T → `core.thermo.entropy_change_reversible` → running ∫dq/T along both paths converging to one value; explainer `heat_work_paths` |
| N15 | T ds = dq for a reversible process (1.17) [#52] | 1.8 | NOTE | differential form of (1.16) | parent C33: equation + one line differentiating the cumulative entropy array |
| C34 | Second law (ii): Clausius–Duhem inequality, s₂ − s₁ ≥ ∫δq/T [#50] | 1.8 | CORE | numbered theorem (typo taught with actual δq) | free expansion and stirring: q = 0 but Δs = R ln 2 > 0 via `perfect_gas_entropy_change` → bar figure Δs vs ∫δq/T for reversible and irreversible processes; mode "irreversible" of `heat_work_paths` |
| N14 | Second law (iii): μ > 0 and k > 0 [#51] | 1.8 | NOTE | consequence of C34 | parent C34: one line: negative μ would un-mix momentum and lower entropy; the argument checks in `newton_shear_stress` raise |
| C35 | Gibbs relations T ds = de + p dv = dh − v dp (1.18) [#53] | 1.8 | CORE | numbered result, load-bearing | D10 → sympy residual cell + `core.thermo.perfect_gas_entropy_change` → T–s diagram showing isotherms and isentropes; explainer `heat_work_paths` |
| C36 | Speed of sound c² = (∂p/∂ρ)_s (1.19); incompressible limit c → ∞ [#54] | 1.8 | CORE | numbered definition | `core.thermo.sound_speed_from_eos` on p = Kρ^γ and a stiff Tait-like liquid EOS → figure c vs bulk stiffness (log) with air and water marked; proof deferred to Ch. 15 (pointer) |
| C37 | Thermal expansion coefficient α = −(1/ρ)(∂ρ/∂T)_p (1.20) [#55] | 1.8 | CORE | numbered definition | `core.thermo.thermal_expansion_coefficient` on air and on a cited water density fit → ρ(T) figure (water max at 4 °C: α changes sign) |
| C38 | Molecular perfect-gas law pV = nk_BT (1.21) [#56] | 1.9 | CORE | numbered law | links C05: from-scratch molecules in a box vs `core.thermo.molecular_gas_pressure` → figure p vs T for three molecule counts |
| C39 | k_B, Avogadro per kmol, R_u = k_BA_o, M_w, R = R_u/M_w [#57] | 1.9 | CORE | new constants | CODATA values in `core.thermo` → tiny example R_air = 8314/28.96 → bar figure of R for He, air, H₂O vapour, CO₂ |
| C40 | Continuum perfect-gas law p = ρRT (1.22) [#58] | 1.9 | CORE | numbered law, load-bearing | D11 → from-scratch chain k_B → R vs `core.thermo.perfect_gas_pressure`, `perfect_gas_density` → USSA sea-level check → isobars ρ(T) figure |
| C41 | A perfect gas has e = e(T), h = h(T) (and conversely) [#62] | 1.9 | CORE | new theorem (book leaves the proof to an exercise) | D13 (★★★, sympy) → figure: e(v) at fixed T flat for a perfect gas vs sloped for van der Waals (`core.thermo.van_der_waals_internal_energy`) |
| C42 | R = C_p − C_v (1.23) [#59] | 1.9 | CORE | numbered result | D12 → `core.thermo.cv_from_cp` → h(T) and e(T) lines separated by RT (figure) |
| C43 | γ = C_p/C_v (1.24) [#60] | 1.9 | CORE | numbered definition | `core.thermo.gamma_from_cp` → bar figure γ for monatomic, diatomic, triatomic gases (5/3, 7/5, ≈ 1.3) |
| N16 | Air values of γ and C_p; C_p and C_v rise with T [#61] | 1.9 | NOTE | values of C43 | parent C43: air marked on the C43 bars with the benchmark note (constant-C assumption flagged for C45–C48) |
| C44 | Adiabatic (no heat) vs isentropic (adiabatic + frictionless) [#63] | 1.9 | CORE | new definitions | p–v figure: isotherm, isentrope, and an adiabatic but stirred (irreversible) path, with Δs readouts from `perfect_gas_entropy_change`; explainer `heat_work_paths` |
| C45 | Isentropic perfect gas p/ρ^γ = const (1.25) [#64] | 1.9 | CORE | numbered result | D14 → from-scratch ds = 0 ODE integration vs `core.thermo.isentropic_pressure` → p–ρ log–log (slope γ) vs isothermal (slope 1); explainer `heat_work_paths` |
| C46 | Isentropic ratios T/T₀ = (p/p₀)^{(γ−1)/γ}, ρ/ρ₀ = (p/p₀)^{1/γ} (1.26) [#65] | 1.9 | CORE | numbered result, load-bearing | D15 → tiny example (bicycle pump, 1 → 2 bar) → `core.thermo.isentropic_ratios` → figure of both ratios vs p/p₀ for γ = 1.4 and 5/3 |
| C47 | Speed of sound in a perfect gas c = √(γRT) (1.27) [#66] | 1.9 | CORE | numbered result (the special case needs C40–C45 taught first, so it cannot sit inside C36) | D16 → 340 m/s at 288 K → `core.thermo.perfect_gas_sound_speed` → c(z) through the USSA troposphere/stratosphere figure |
| C48 | α = 1/T for a perfect gas (1.28) [#67] | 1.9 | CORE | numbered result (needs C40 first, so not inside C37) | D17 → `perfect_gas_expansion_coefficient` vs the numerical C37 function → α(T) curve with water's α for contrast |
| C49 | Static medium: (1.8) with (1.12) tie p, ρ, T; one profile fixes the others [#68] | 1.10 | CORE | new idea | `core.statics.atmosphere_from_temperature` with the USSA lapse −6.5 K/km → three-panel p, ρ, T profile figure vs USSA points |
| C50 | Displaced-parcel equation ζ'' − (g/ρ)(dρ/dz − dρ_a/dz)ζ = 0 [#69] | 1.10 | CORE | new result (book leaves it to an exercise) | D18 → from-scratch Euler–Cromer loop vs `core.stratification.parcel_displacement` and `parcel_ode` → animation of three parcels (stable/neutral/unstable); explainer `parcel_stability` |
| C51 | Brunt–Väisälä frequency N² (1.29) [#70] | 1.10 | CORE | numbered definition, load-bearing | tiny example: ocean thermocline N ≈ 10⁻² s⁻¹ (period ≈ 10 min) → `core.stratification.brunt_vaisala_sq` → ζ(t) figure for several N² with period 2π/N marked; explainer `parcel_stability` |
| C52 | Stable (N² > 0), neutral (= 0), unstable (< 0) [#71] | 1.10 | CORE | new criterion | vectorised `classify_stability` → a ρ(z) profile coloured by regime; explainer `parcel_stability` status line |
| C53 | Lapse rate Γ ≡ dT/dz [#72] | 1.10 | CORE | new definition (the book's sign convention is used throughout) | T(z) figure with slope triangles: troposphere Γ = −6.5 K/km (negative, cooling upward), isothermal layer Γ = 0, inversion Γ > 0; code `np.gradient(T, z)` gives Γ directly. One-line side note only: some meteorology texts quote the magnitude as a positive number. |
| C54 | Adiabatic lapse rate Γ_a = −gαT/C_p (1.30) [#73] | 1.10 | CORE | numbered result | D19 (★★★, sympy) → `core.stratification.adiabatic_lapse_rate` → environment T(z) vs parcel adiabats figure; explainer `parcel_stability` |
| N17 | Value of Γ_a for Earth's atmosphere ≈ −9.8 K/km [#74] | 1.10 | NOTE | value of C54 | parent C54: computed Γ_a = −g/C_p printed in K/km (≈ −9.8 K/km); the benchmark test compares its magnitude with the published 9.8 K/km value |
| C55 | Potential temperature θ (1.31) [#76] | 1.10 | CORE | numbered definition, load-bearing (Ch. 13) | D20 → tiny example (air at 500 hPa, 250 K → θ ≈ 305 K) → from-scratch formula vs `core.stratification.potential_temperature` → T(z) and θ(z) panels; explainer `parcel_stability` |
| C56 | (T/θ)dθ/dz = dT/dz + g/C_p = Γ − Γ_a (1.32) [#78] | 1.10 | CORE | numbered result | D21 → `potential_temperature_gradient` vs finite differences of θ(z) → figure: dθ/dz and Γ − Γ_a overlaid |
| N19 | Log-derivative of (1.31) [#77] | 1.10 | NOTE | intermediate step of C56 | parent C56: step 2–3 of D21, with a one-line numeric check |
| C57 | Stability from the sign of dθ/dz; lab-scale T ≈ θ; N² = (g/θ)dθ/dz [#79] | 1.10 | CORE | new criterion (N² link added) | D36 (added) → synthetic boundary-layer profile → plotly slider on the inversion strength: T(z), θ(z), N²(z) coloured stable/unstable; explainer `parcel_stability` |
| N18 | Fig. 1.9a: T(z) with near-neutral layer, inversion, stable layer, neutral reference lines [#75] | 1.10 | NOTE | figure of C57 | parent C57: panel (a) of the C57 figure (our synthetic profile, not digitised) |
| C58 | Potential density ρ_θ (1.33) [#80] | 1.10 | CORE | numbered definition | D22 → `core.stratification.potential_density` → ρ(z) vs ρ_θ(z) figure for the C57 profile |
| C59 | −(1/ρ_θ)dρ_θ/dz = (1/θ)dθ/dz (1.34) [#82] | 1.10 | CORE | numbered result | D23 → finite-difference check of both sides → overlay figure of the two log-gradients |
| N20 | θρ_θ = p_o/R = const [#81] | 1.10 | NOTE | intermediate result of C59 | parent C59: step of D23 plus `np.allclose(theta*rho_theta, p_ref/R)` line |
| C60 | Seawater: salinity S, ρ = ρ(T, p, S), potential density at constant S [#83] | 1.10 | CORE | new quantity (important for the ocean) | tiny example: S = 35 g/kg → `ch01.seawater_density_linear` (cited linear EOS) → ρ(T, S) contour figure with typical water masses |
| C61 | Ocean stability criterion dρ_θ/dz ≅ dρ/dz + ρg/c² (1.35) [#85] | 1.10 | CORE | numbered criterion | D24 → tiny example ρg/c² ≈ 4.5×10⁻³ kg m⁻⁴ → `ocean_potential_density_gradient` → deep-ocean profile that looks stable in ρ but is near-neutral in ρ_θ (figure) |
| N21 | Parcel density gradient dρ_a/dz = −ρg/c² in the ocean [#84] | 1.10 | NOTE | first half of C61's derivation | parent C61: steps 1–4 of D24; `isentropic_density_gradient` line; dashed reference line on the C61 figure |
| C62 | Isothermal atmosphere p = p₀e^{−gz/RT} [#86] | 1.10 | CORE | new result | D25 → `core.statics.isothermal_pressure` vs `integrate_hydrostatic` → plotly slider on T: exponential p(z) vs USSA table |
| C63 | Scale height H = RT/g [#87] | 1.10 | CORE | new quantity | tiny example T = 250 K → H ≈ 7.3 km → `core.statics.scale_height` → semilog p(z) figure with e-folds marked at H, 2H, 3H |
| C64 | Dimensional homogeneity; laws can be written dimensionless [#88] | 1.11 | CORE | new principle, load-bearing | `core.units.dimensional_check` on (1.9) and on a wrong formula → figure: dimension-exponent bars per term (match vs mismatch) |
| C65 | Step 1: choose variables; pipe problem f(Δp, Δx, d, ε, U, ρ, μ) = 0 (1.38) [#91] | 1.11 | CORE | method step with a trap (a missing variable breaks the result) | variable dictionary in code → figure: synthetic laminar pipe data plotted with and without μ in the group list (no collapse without it) |
| C66 | Step 2: [q] notation, base dimensions M, L, T, θ [#92] | 1.11 | CORE | new notation | `core.dimensional.dimension_vector` via pint → heatmap of exponent vectors for the chapter's quantities (p, μ, ν, σ, C_p, R, k_B·θ) |
| C67 | Dimensional matrix of the pipe problem (1.39) [#93] | 1.11 | CORE | numbered object | `core.dimensional.dimensional_matrix` → annotated heatmap of (1.39); explainer `buckingham_pi_machine` |
| C68 | Step 3: rank as the largest nonzero minor [#94] | 1.11 | CORE | new maths tool | D26 → from-scratch cofactor expansion over `itertools.combinations` vs `core.dimensional.rank_by_minors` and `np.linalg.matrix_rank` → figure highlighting the zero and nonzero 3×3 minors; explainer `buckingham_pi_machine` |
| C69 | Buckingham's Π theorem φ(Π₁ … Π_{n−r}) = 0 (1.37) [#90] | 1.11 | CORE | numbered theorem, load-bearing | D28 (★★★, sympy nullspace) → `core.dimensional.pi_groups` → figure: groups keep their values when units change (cgs vs SI scatter on the identity line); explainer `buckingham_pi_machine` |
| N22 | The relation f(q₁ … q_n) = 0 (1.36) [#89] | 1.11 | NOTE | starting form of C69 | parent C69: the input dictionary of `pi_groups` |
| N23 | Step 4: number of groups is n − r [#95] | 1.11 | NOTE | part of C69's statement | parent C69: one printed line `n - r` for the pipe problem (7 − 3 = 4) |
| C70 | Step 5: exponent algebra with a repeating set [#96] | 1.11 | CORE | new method | D27 → from-scratch `np.linalg.solve` on the 3×3 system vs `core.dimensional.solve_exponents` → figure of the linear system (matrix · exponents = target) for Π₁; explainer `buckingham_pi_machine` |
| N24 | Step 5 by inspection: remove M, L, T one at a time [#97] | 1.11 | NOTE | alternative method for C70 | parent C70: the ratio chain Δp → Δp/ρ → Δp/ρU² shown as three printed lines checked with `dimension_vector` |
| C71 | Combining groups: products and powers are groups; only n − r independent [#98] | 1.11 | CORE | new idea | D33 → `core.dimensional.groups_independent` → exponent-matrix heatmaps of an independent and a dependent set |
| C72 | Step 6: pipe law Δp/ρU² = φ(Δx/d, ε/d, μ/ρUd) (1.40) [#99] | 1.11 | CORE | numbered result | synthetic laminar data (many d, U, μ) collapses onto one line in Π coordinates (figure); explainer preset in `buckingham_pi_machine` |
| N25 | Step 7: use physics (proportionality, linearity) to simplify [#100] | 1.11 | NOTE | follow-on step of C72 | parent C72: Δp ∝ Δx turns (1.40) into Π₁ = (Δx/d)φ₁(ε/d, Re); the collapse figure redrawn in the reduced variables; forward pointer to C73–C76 which use it |
| C73 | Ex. 1.2: scale height by dimensional analysis [#101] | 1.11 | CORE | worked result | D29 → `pi_groups` on (H, T_o, M_w, g, R_u) → H vs T by both routes (C63 line vs Π result) on one figure |
| C74 | Ex. 1.3: Pythagoras by dimensional analysis [#102] | 1.11 | CORE | worked result | D30 → `pi_groups` on (a, β, C) and `ch01.pythagoras_phi` → triangle split into similar triangles figure + numeric check of A² + B² = C² |
| C75 | Ex. 1.4: G. I. Taylor's blast energy E = KρD⁵/t² [#103] | 1.11 | CORE | worked result | D31 → `ch01.blast_energy`, `blast_radius` → log–log fit of synthetic D(t) (slope 2/5) + growing-hemisphere animation; preset in `buckingham_pi_machine` |
| C76 | Ex. 1.5: Rayleigh scattering S/I = V²φ₃(n_s)/(d²λ⁴) [#104] | 1.11 | CORE | worked result | D32 → `ch01.rayleigh_scattering_ratio` → plotly slider on particle size: S/I across the visible spectrum coloured by wavelength; blue/red ratio 5.85 |
| S01 | Exercises 1.1–1.30 [#105] | Ex | SKIP | exercises | pointer: "The chapter's 30 exercises are in the book and are not reproduced here. The derivations the book leaves to Exercises 1.10, 1.11, 1.13 and 1.14 are written out above (D13, D14, D18, D19)." |
| S02 | Literature cited and supplemental reading [#106] | Ref | SKIP | bibliography | pointer: "The book's reading list is at the end of the printed chapter; the public data we used (CODATA, USSA-1976, IAPWS, Jennings, Taylor) are cited in reference/ch01/SOURCES.md." |

## 3. Section coverage
| § | Title | CORE | RECAP | NOTE | SKIP |
|---|---|---|---|---|---|
| 1.1 | Fluid Mechanics | C01 | — | — | — |
| 1.2 | Units of Measurement | — | R01, R02, R03 | — | — |
| 1.3 | Solids, Liquids, and Gases | C02, C03, C04 | — | N01, N02 | — |
| 1.4 | Continuum Hypothesis | C05, C06, C07 | — | N03 | — |
| 1.5 | Molecular Transport Phenomena | C08, C09, C10, C11, C12, C13, C14 | — | N04, N05, N06 | — |
| 1.6 | Surface Tension | C15, C16 | — | N07, N08, N09 | — |
| 1.7 | Fluid Statics | C17, C18, C19, C20, C21, C22 | — | N10, N11, N12 | — |
| 1.8 | Classical Thermodynamics | C23, C24, C25, C26, C27, C28, C29, C30, C31, C32, C33, C34, C35, C36, C37 | — | N13, N14, N15 | — |
| 1.9 | Perfect Gas | C38, C39, C40, C41, C42, C43, C44, C45, C46, C47, C48 | — | N16 | — |
| 1.10 | Stability of Stratified Fluid Media | C49, C50, C51, C52, C53, C54, C55, C56, C57, C58, C59, C60, C61, C62, C63 | — | N17, N18, N19, N20, N21 | — |
| 1.11 | Dimensional Analysis | C64, C65, C66, C67, C68, C69, C70, C71, C72, C73, C74, C75, C76 | — | N22, N23, N24, N25 | S01, S02 (end of chapter) |

## 4. Prerequisites needing primers (concept or tool | needed by | why it is not CORE/RECAP)
Nothing is primed yet (`knowledge/primers.md` is empty), so every entry below needs a 📎 primer the first time it is used.

| Concept or tool | Needed by | Why it is not CORE/RECAP |
|---|---|---|
| **Physics vocabulary** | | |
| stress = force per area; normal and tangential (shear) components; strain angle γ and strain rate | C02, C03, C12 | pre-book mechanics, used without definition in §1.3 |
| elastic shear modulus G (solid: τ = Gγ) | C02 | not in the book; our model for Fig. 1.1 |
| vapour pressure (why a liquid under tension boils/cavitates) | C03 | mentioned in passing, taught in no chapter |
| Newton's second law F = ma and free-body diagrams | C16, C18, C20, D03, D18 | pre-book |
| weight = mg, gravitational acceleration g | C20, C21 | pre-book |
| temperature as average molecular kinetic energy; Boltzmann constant k_B (informal, before C38) | C05, C23 | used in §1.4 before §1.9 defines it |
| mole, kilomole, molecular weight, Avogadro's number (informal, before C39) | C04, C09 | needed for molecular spacing and mass fractions before §1.9 |
| perfect-gas law as a working model (forward look at C40) | C26, C27, C28, C33 | the thermodynamics demos need a concrete e and p before §1.9 |
| specific heat (informal) and thermal diffusivity κ = k/ρC_p | C08, C11 | C_p is defined only in §1.8 (C30) |
| boundary condition, no-slip at a wall, steady vs transient | C08, C12, explainer `viscosity_momentum_diffusion` | formal treatment in Ch. 4 |
| internal energy and kinetic energy per unit mass | C25 | pre-book physics; the book uses e without defining it |
| light as waves: wavelength and colour; intensity ∝ amplitude²; oscillating dipole | C76 | pre-book optics, needed for Ex. 1.5's physical arguments |
| similarity solution (why K comes from a full solution) | C75 | full treatment is beyond Ch. 1; one sentence |
| **Mathematics** | | |
| ordinary derivative as a slope; d/dy of a profile | C12, C53 | pre-book calculus |
| partial derivative ∂/∂x | C19, C20, D04, D05 | pre-book calculus, first used in (1.7) |
| partial derivative with a variable held fixed, ( )_p notation | C30, C31, C36, C37 | thermodynamic notation, not school calculus |
| total vs partial derivative | D12 | needed when e(T, v) becomes e(T) |
| gradient vector ∇ and "down-gradient" | C10, C11 | vector calculus comes formally in Ch. 2 |
| dot product and unit normal vector | C03, C10 | Ch. 2 formalises it |
| first-order Taylor expansion f(x + dx) ≈ f + f′dx | D04, D05, D18 | pre-book |
| limits and orders of smallness (dz² ≪ dz) | D03 | pre-book, but the key move of D03 |
| small-angle approximation sin x ≈ x | D02 | pre-book trigonometry |
| radius of curvature; principal radii; signed curvature | C16, D02 | differential geometry, not taught elsewhere |
| definite integral; integrating a constant | D06, D25, D37 | pre-book |
| surface integral of pressure over a closed body (net force) | D37, D01 | needed for Archimedes and projected-area arguments |
| natural log and exponential; e-folding | C62, C63, D14, D25 | pre-book |
| separation of variables for a first-order ODE | D14, D25 | pre-book ODEs |
| exponent rules (a^m a^n = a^{m+n}) | D15, D20, D22, D23 | pre-book algebra |
| logarithmic differentiation d ln f = df/f | D21, D23, D36 | not school-standard |
| differentials; exact (d) vs inexact (δ) differentials | C25, C26, C27 | the book's δ/d notation needs a primer |
| line integral along a path (∫p dv as area on a p–v diagram) | C26, C33 | calculus-of-paths idea |
| product rule for differentials d(pv) = p dv + v dp | D09, D10 | pre-book, easy to drop a term |
| chain rule (one and several variables) | D19, D24 | pre-book |
| equality of mixed partials, exact differentials, Maxwell relations | D13, D19 | the chapter never introduces them (analysis §9) |
| Gibbs free energy g = h − Ts (as a device) | D19 | not in Ch. 1, needed for the Maxwell relation |
| linear second-order ODE ζ″ + N²ζ = 0; cos, sinh/cosh solutions | C50, C51, D18 | pre-book ODEs |
| square root of a negative number → growth rate (imaginary N) | C51, C52 | reading N² < 0 |
| Poisson counting statistics; relative noise ∝ N^{-1/2} | C06, D35 | statistics tool |
| mean square speed; Gaussian (Maxwell) velocity components | C05, D34 | kinetic theory tool |
| power laws and log–log plots | C06, C13, C45, C75 | reading slopes as exponents |
| finite differences: central difference, FTCS scheme, stability limit DΔt/Δy² ≤ ½ | C08, C12, C30 | numerical tool (Ch. 10 formal) |
| trapezoid rule | C26, C33 | numerical integration |
| explicit ODE time stepping (Euler, Euler–Cromer) | C20, C50 | from-scratch versions |
| matrices, determinants, minors, cofactor expansion | C67, C68, D26 | linear algebra |
| linear independence and rank | C68, C71 | linear algebra |
| solving linear systems Ax = b | C70, D27, D29 | linear algebra |
| null space and rank–nullity theorem | C69, D28 | linear algebra |
| similar triangles and angle chasing | D30 | pre-book geometry |
| **Python** | | |
| numpy arrays, vectorised arithmetic, broadcasting | C02 onward | first Python cell |
| `np.linspace`, `np.logspace` | C06, C07 | log-spaced box sizes and bodies |
| `np.random.default_rng(seed)`: `normal`, `poisson` | C05, C06 | seeded sampling |
| `np.gradient` (and `edge_order=2`) | C12, C19 | finite differences on arrays |
| vectorised `np.where` / `np.select` | C52, C57 | regime labels |
| `matplotlib`: `subplots`, `plot`, `loglog`, `semilogy`, `fill_between`, `quiver`/`arrow`, `imshow`/`pcolormesh` | C02 onward | first figure cells |
| `fluidpy.core.anim.animate` + `show_animation` | C02, C08, C50 | first animation |
| `fluidpy.core.interact.slider_figure` (plotly) | C07 | first slider figure |
| `plotly.graph_objects.Surface` | C28, N09 | 3-D surfaces |
| `pint`: `Q_`, `.to()`, `.dimensionality`, offset units | R01, R03, C64, C66 | units in code |
| `scipy.integrate.solve_ivp` (dense output, events) | C20, C49, C50 | ODE solver |
| `scipy.integrate.cumulative_trapezoid` / `np.trapezoid` | C26, C33 | path integrals |
| `sympy`: `symbols`, `diff`, `simplify`, `Matrix`, `nullspace`, `Rational`, `linsolve` | D13, D19, D28, C67–C70 | symbolic checks |
| `itertools.combinations` | C68 | enumerating minors |
| `np.linalg.det`, `matrix_rank`, `solve` | C68, C70 | numerical linear algebra |
| `np.polyfit` on logarithms | C75 | power-law fit |
| functions passed as arguments / `lambda` | C30, C36, C37, C20 | `partial_derivative(f, …)`, `integrate_hydrostatic(rho_fn)` |
| dictionaries as variable lists | C65 | `pi_groups` input |
| f-strings with units; `assert np.allclose` | C02 onward | from-scratch checks |
| `show_viz(chapter, slug)` | first explainer (C06) | embedding explainers |

## 4b. Derivations (parsed by tools: ID first, CORE id in a column, ★★★ for hard, explainer slugs backticked in the LAST column)
| ID | Result (Eq.) | CORE | Difficulty | Steps | Tools used | Traps | Shown in |
|---|---|---|---|---|---|---|---|
| D01 | spherical jump p_i − p_o = 2σ/R (unnumbered) | C16 | ★ | 6 | free-body diagram (primer), surface integral of pressure → projected area (primer) | using the curved area 2πR² instead of the projected πR²; which side is higher (the concave side) | notebook |
| D02 | Laplace jump (1.5) (book: "a similar analysis"; all moves filled in) | C16 | ★★ | 9 | small-angle sin x ≈ x (primer), radius of curvature (primer) | counting each pair of edges once, not twice; unsigned radii for a saddle; forgetting the patch area ds₁ds₂ cancels | notebook |
| D03 | pressure is isotropic (1.6) | C18 | ★★ | 10 | force components and trigonometry (primer), Newton's second law (primer), limits and orders of smallness (primer) | keeping the weight term in the limit; dropping the y-balance on the end faces; forgetting there is no shear at rest (C02) | notebook |
| D04 | Pascal's law (1.7) | C19 | ★ | 5 | first-order Taylor (primer), partial derivative (primer) | a sign on the p + (∂p/∂x)dx face; dividing by the volume before cancelling | notebook |
| D05 | hydrostatic law (1.8) | C20 | ★ | 5 | first-order Taylor (primer), partial derivative (primer), C19 | sign of g with z upward; writing ∂p/∂z before C19 shows p depends on z only | notebook |
| D06 | uniform-density hydrostatics (1.9) and p − p₀ = ρgh | C21 | ★ | 5 | definite integral (primer) | h = −z, so the sign flips; p₀ is the value at z = 0, not the atmosphere by default | notebook |
| D07 | capillary rise h = 2σ sin α/(ρgR) (Ex. 1.1) | C22 | ★ | 7 | C16, C21, force balance (primer) | book's α is the complement of the usual contact angle; meniscus volume neglected; p_E below atmospheric | notebook |
| D08 | reversible first law (1.11) | C27 | ★ | 6 | differentials (primer), C25 | sign: work *on* the particle is −p dv; including frictional work | notebook · `heat_work_paths` |
| D09 | heat per degree equals C_v (constant v) and C_p (constant p) (unnumbered) | C32 | ★★ | 8 | product rule for differentials (primer), C27, C29 | dropping v dp in dh; treating Q as a state function; capital-Q notation for a per-mass quantity | notebook |
| D10 | Gibbs relations (1.18) | C35 | ★ | 5 | product rule for differentials (primer), N15 (T ds = dq), C27, C29 | concluding they hold only for reversible processes (they link state functions, so they hold for any process) | notebook · `heat_work_paths` |
| D11 | p = ρRT (1.22) from pV = nk_BT (1.21) | C40 | ★ | 6 | C38, C39, unit bookkeeping with kmol (primer) | mixing mol and kmol (factor 1000); R_u vs R | notebook |
| D12 | R = C_p − C_v (1.23) | C42 | ★★ | 7 | total vs partial derivative (primer), C41, C29 | differentiating pv as if p were constant; using it for a non-perfect gas | notebook |
| D13 | perfect gas ⇔ e = e(T), h = h(T) (claim; Exercise 1.10) | C41 | ★★★ | 14 | exact differentials and mixed partials (primer), Maxwell relations (primer), C35, C40, first-order PDEs (primer) | sign of the Maxwell relation; forgetting the converse direction needs integrating two PDEs; confusing (∂e/∂v)_T with C_v | notebook |
| D14 | isentropic law p/ρ^γ = const (1.25) (Exercise 1.11) | C45 | ★★ | 9 | C35, C42, C43, separation of variables (primer), logarithms (primer) | dv/v = −dρ/ρ sign; assuming constant C_p, C_v silently; dividing the two Gibbs forms in the wrong order | notebook · `heat_work_paths` |
| D15 | isentropic ratios (1.26) | C46 | ★ | 5 | exponent rules (primer), C40, C45 | 1 − 1/γ vs (γ − 1)/γ written inconsistently; inverting ρ₀/ρ | notebook |
| D16 | c = √(γRT) (1.27) | C47 | ★ | 5 | differentiating a power (primer), C36, C45, C40 | differentiating at constant T (Newton's error, gives √(RT)) instead of constant s | notebook |
| D17 | α = 1/T (1.28) | C48 | ★ | 4 | partial derivative held fixed (primer), C37, C40 | losing the minus sign so α comes out negative | notebook |
| D18 | parcel equation ζ″ + N²ζ = 0 and (1.29) (Exercise 1.13) | C50 | ★★ | 11 | Newton's second law (primer), D37 buoyancy (C21), first-order Taylor (primer), linear second-order ODE (primer) | sign of buoyancy; using the environment's gradient for the parcel; keeping O(ζ²) terms; ρ_p ≈ ρ(z_o) only in the inertia term | notebook · `parcel_stability` |
| D19 | adiabatic lapse rate (1.30) (Exercise 1.14, no perfect-gas relations) | C54 | ★★★ | 13 | Gibbs free energy (primer), Maxwell relations (primer), chain rule (primer), C35, C20, C37, C30 | Maxwell-relation sign; using dh = C_p dT (true only for a perfect gas); dropping the minus sign: with Γ ≡ dT/dz the result Γ_a = −gαT/C_p is negative | notebook · `parcel_stability` |
| D20 | potential temperature (1.31) | C55 | ★ | 4 | exponent rules (primer), C46 | p_o (reference pressure) confused with p₀ at z = 0 | notebook |
| D21 | dθ/dz and lapse rates (1.32) | C56 | ★★ | 9 | logarithmic differentiation (primer), C20, C40, C42, C43, C54 | (γ − 1)/γ = R/C_p identity not shown; Γ − Γ_a = dT/dz + g/C_p because Γ_a = −g/C_p is negative; forgetting T/θ > 0 so signs agree | notebook |
| D22 | potential density (1.33) | C58 | ★ | 4 | exponent rules (primer), C46 | exponent 1/γ vs (γ − 1)/γ | notebook |
| D23 | θρ_θ = p_o/R and (1.34) | C59 | ★★ | 7 | exponent rules (primer), logarithmic differentiation (primer), C40 | adding exponents (γ − 1)/γ + 1/γ = 1; the minus sign in (1.34) | notebook |
| D24 | ocean criterion (1.35) with dρ_a/dz = −ρg/c² | C61 | ★★ | 8 | chain rule (primer), C36, C20, inverse derivative (primer) | (∂ρ/∂p)_s = 1/c², not c²; the parcel follows its surroundings' pressure; (1.35) holds as "same sign", not equality | notebook |
| D25 | isothermal p = p₀e^{−gz/RT} and H = RT/g | C62 | ★ | 6 | separation of variables (primer), logarithms and e-folding (primer), C40 | treating ρ as constant; ln of a dimensional quantity (use p/p₀) | notebook |
| D26 | rank of the pipe matrix (1.39) is 3 | C68 | ★ | 6 | determinants and cofactor expansion (primer), rank (primer) | one zero minor does not mean rank < 3; rank ≤ number of rows | notebook · `buckingham_pi_machine` |
| D27 | Π₁ = Δp/ρU², Π₂ = Δx/d, Π₃ = ε/d, Π₄ = μ/ρUd | C70 | ★★ | 10 | linear systems (primer), C67, C68 | sign in the L-equation (−3c); choosing a singular repeating set; skipping the final dimension check | notebook · `buckingham_pi_machine` |
| D28 | Buckingham Π theorem (1.37) (book states it only) | C69 | ★★★ | 12 | null space and rank–nullity (primer), C64, C67, C68, scaling of base units (primer) | thinking the groups are unique; forgetting the relation must be dimensionally homogeneous and the variable list complete | notebook · `buckingham_pi_machine` |
| D29 | Ex. 1.2: H = const·R_uT_o/(gM_w), const = 1 | C73 | ★★ | 9 | linear systems (primer), C63, C69 | the kmol dimension dropped without saying why; a single group must be a constant (not zero) | notebook |
| D30 | Ex. 1.3: A² + B² = C² | C74 | ★★ | 8 | similar triangles and angle chasing (primer), C69 | needs φ(β) ≠ 0; each sub-triangle's hypotenuse is A or B, not C | notebook |
| D31 | Ex. 1.4: E = KρD⁵/t² and D ∝ t^{2/5} | C75 | ★ | 6 | exponent rules (primer), logarithms (primer), C69 | expecting dimensional analysis to give K; forgetting ambient pressure was dropped | notebook |
| D32 | Ex. 1.5: S/I = V²φ₃(n_s)/(d²λ⁴) | C76 | ★★ | 8 | power laws (primer), C69, rank (primer) | why r = 2 (T row = −3 × M row); d⁻² must enter as (λ/d)²; V² forces λ⁻⁴ | notebook |
| D33 | combined groups Δp d²ρ/μ² = Π₁/Π₄², ε/Δx = Π₃/Π₂ | C71 | ★ | 4 | exponent rules (primer), rank (primer) | counting a combined group as an extra independent group | notebook |
| D34 | kinetic pressure p = ⅓ (N/V) m⟨u²⟩ (added: the book states the idea only) | C05 | ★★ | 8 | Newton's second law (primer), Gaussian velocity components and mean square speed (primer) | factor 2 from the bounce; ⟨u_x²⟩ = ⅓⟨u²⟩; counting only molecules moving toward the wall | notebook |
| D35 | relative density noise ≈ N^{-1/2} ∝ (δV)^{-1/2} (added) | C06 | ★ | 5 | Poisson counting statistics (primer), power laws (primer) | noise ∝ L^{-3/2} in box side, not L^{-1/2}; mixing absolute and relative noise | notebook · `continuum_averaging_volume` |
| D36 | N² = (g/θ)dθ/dz for a perfect-gas atmosphere (added: links (1.29) and (1.32)) | C57 | ★★ | 8 | logarithmic differentiation (primer), C46, C51, C56, C42, C43 | using the environment's density gradient for the parcel; (γ − 1)/γ = R/C_p identity; sign of dp/dz | notebook · `parcel_stability` |
| D37 | buoyancy force = ρ_fluid gV from p = p₀ − ρgz (added: grounds D18) | C21 | ★ | 6 | surface integral of pressure over a closed body (primer), definite integral (primer), C21 | side-face forces do cancel; the net force is up; body's own density is irrelevant to the buoyancy | notebook |

## 5. Interactive explainers (4–5 + backup)

### E1 · continuum_averaging_volume
- **CORE:** C06, C07, C24 (+ C05 in the molecules view)
- **Confusion removed:** "density at a point" sounds meaningless when matter is made of molecules. The reader sees that a well-defined value exists only inside a *window* of box sizes: larger than the molecular spacing, smaller than the scale over which the flow itself varies.
- **Why interactive:** the idea is a sweep over eight decades of box size with random noise. Dragging the box and resampling shows noise that no static curve conveys: the same box can read 30 % high one moment and 20 % low the next, then settles as it grows.
- **Stage:** (1) a zoomable cross-section of gas or liquid with molecules as dots and the sampling cube drawn on it, with a background density gradient visible at large zoom; (2) measured ρ vs box side L on log axes: a faint expected band (±N^{-1/2}), bold samples so far, the plateau and the macroscopic drift; (3) Kn = l/L strip with continuum/slip/free-molecular bands.
- **Controls:** box side L (log slider) · medium (chips: air at sea level, air at 80 km, water) · macroscopic gradient strength · "resample" button · body size for Kn (optional).
- **Equations shown:** ρ = δm/δV (continuum definition, §1.4), relative noise ≈ N^{-1/2} (D35), Kn = l/L (§1.4); kinetic pressure (D34) in the Equations tab only.
- **Mirrors:** `fluidpy.ch01_introduction.sample_density`, `density_noise_expected`, `knudsen_number`, `mean_free_path_jennings`.
- **Derivations:** D35.
- **Depth features:** Explain tab (N = nL³ computed with your numbers → noise → verdict), synced Code tab, + linked views (3), presets (air / high altitude / water / microchannel), status verdict ("🎲 molecular noise ±x %" / "✅ continuum plateau" / "📈 box sees the flow's own variation"), inspector (click a sample: count, mass, δm/δV arithmetic), "Right now" notes with a highlighted table of real sizes (virus, aerosol, raindrop, pipe, cloud).
- **Follows reference:** `angular_frequency_explorer_1.html` (linked views, modes, highlighted real-world table, Right-now notes).
- **Aha:** a continuum value is not a point property of matter. It is the plateau of an average, and it exists because there is a wide range of box sizes between the molecules and the flow.

### E2 · viscosity_momentum_diffusion
- **CORE:** C08, C12, C14 (+ C10, C11 through modes)
- **Confusion removed:** viscosity pictured as "stickiness" or friction between layers, and μ confused with ν. Momentum *diffuses* from the moving plate, and ν = μ/ρ, not μ, sets how fast. Air (small μ) spreads motion faster than water.
- **Why interactive:** the phenomenon is a transient: a profile creeping into the gap and settling to the linear Couette state with uniform τ. Playing, scrubbing and switching fluids shows the time scale h²/ν at work, and a static figure cannot show the process.
- **Stage:** (1) the gap between plates with dyed vertical tracer lines shearing and a few molecules hopping between layers carrying momentum colour; (2) u(y, t) with a faint steady linear profile (ghost) and the bold current profile; (3) wall shear stress τ_w(t) = μ ∂u/∂y at both plates, both tending to μU/h, with a marker at t = h²/ν.
- **Controls:** fluid preset or μ (log) · ρ · gap h · plate speed U · mode (momentum / heat / species).
- **Equations shown:** (1.3), (1.4), (1.1), (1.2); model diffusion equation ∂f/∂t = D∂²f/∂y² (labelled "ours; derived in Ch. 4").
- **Mirrors:** `fluidpy.core.diffusion.ftcs_diffusion_1d`, `fluidpy.ch01_introduction.newton_shear_stress`, `kinematic_viscosity`, `thermal_diffusivity`.
- **Derivations:** none (the book derives nothing in §1.5; the Explain tab works τ, ν and h²/ν out with numbers).
- **Depth features:** Explain tab, synced Code tab, + transport (play/step/scrub, end-of-run summary card), linked views (3), modes (momentum / heat / species with D = ν, k/ρC_p, κ_m), presets (air, water, glycerine, honey), status verdict ("⏳ momentum still spreading: t/(h²/ν) = 0.12" / "✅ steady Couette: τ uniform").
- **Follows reference:** `forced_damped_vibrations.html` (system + graph on one clock, transient → steady state, "at the current time" values, regime-dependent interpretation).
- **Aha:** honey has about 10⁴ times the μ of air, but ν decides how fast motion spreads, and air's ν is about 15 times water's.

### E3 · heat_work_paths
- **CORE:** C25, C26, C27, C33, C34, C35, C44, C45
- **Confusion removed:** heat and work treated as things a gas "has"; δ vs d; why entropy is a property although heat is not; adiabatic vs isentropic.
- **Why interactive:** path dependence only clicks when the reader drags a path between the *same* two states and watches q and w change while Δe and Δs stay fixed. Term bars updating live beat any pair of static diagrams.
- **Stage:** (1) piston–cylinder with gas: the piston moves along the chosen path, a heater glows when δq > 0, an insulation jacket appears for adiabatic legs, a stirrer runs in irreversible mode; (2) p–v diagram with isotherms (faint), the chosen path bold "so far", and the work area shaded; (3) T–s diagram with the same path, with its heat area shaded.
- **Controls:** path (chips: isothermal · isochoric→isobaric · isobaric→isochoric · isentropic + isochoric · custom corner) · end state v₂/v₁ · T₂/T₁ · gas (γ = 5/3 or 7/5) · mode (reversible / irreversible: stirring or free expansion).
- **Equations shown:** (1.10), (1.11), (1.13), (1.16), (1.17), (1.18), (1.25), and the Clausius–Duhem inequality (with actual δq).
- **Mirrors:** `fluidpy.core.thermo.process_heat_work`, `entropy_change_reversible`, `perfect_gas_entropy_change`, `isentropic_pressure`.
- **Derivations:** D08, D10, D14.
- **Depth features:** Explain tab (w = −∫p dv leg by leg with your numbers, q = Δe − w, Δs two ways), synced Code tab, + linked views (3), transport along the path, term bars (q, w, Δe, and Δs vs ∫δq/T, which sum or compare live), presets (the five classic paths), modes (reversible / irreversible), status ("🔁 same Δe, different q" / "⚠️ irreversible: Δs > ∫δq/T").
- **Follows reference:** `amplitude_phase_second_order_II_3.html` (several windows linked by one state, numbered live derivation in the explanation), with term bars after `fid_formula_lab.html`.
- **Aha:** along every path between the same two states, Δe and Δs are identical while q and w are not. That is what "state function" means.

### E4 · parcel_stability
- **CORE:** C50, C51, C52, C53, C54, C55, C57 (C49 as the background builder)
- **Confusion removed:** "colder air aloft means unstable" (wrong: the standard troposphere cools with height and is stable). Stability compares the environment's lapse rate with the parcel's own adiabatic cooling, which θ does in one step. N² < 0 means exponential runaway, not a faster oscillation.
- **Why interactive:** the reader has to *push a parcel* and watch it oscillate, stay put or escape while the environment profile is reshaped. The link between a profile's slope, the parcel adiabat and ζ(t) is a three-way dependence that sliders make obvious.
- **Stage:** (1) a vertical column coloured by θ (or ρ) with a parcel that moves on the clock, drawn with its own adiabat through the release height; (2) the profile panel: environment T(z) bold, parcel dry adiabat dashed (or ρ and ρ_a in ocean mode), θ(z) toggle; (3) ζ(t) with the linear solution (cos / linear / cosh) as ghost and the nonlinear `parcel_ode` path bold, period 2π/N or e-folding time marked.
- **Sign convention (user decision):** the book's Γ ≡ dT/dz everywhere: slider, readouts, equations, Explain and Derivation tabs. The slider runs from about −15 to +10 K/km, the dry adiabat sits at Γ_a ≈ −9.8 K/km, and "stable" means Γ > Γ_a. At most one Explain-tab hint line: "some meteorology texts quote the magnitude, 9.8 K/km, as a positive number".
- **Controls:** environment lapse rate Γ = dT/dz in K/km, negative when T falls with height (or dρ/dz in ocean/lab mode) · inversion toggle/strength · release displacement ζ₀ · mode (atmosphere / ocean / lab tank) · show θ.
- **Equations shown:** parcel equation (§1.10), (1.29), (1.30), (1.31), (1.32), (1.35) in ocean mode, N² = (g/θ)dθ/dz.
- **Mirrors:** `fluidpy.core.stratification.brunt_vaisala_sq`, `parcel_displacement`, `parcel_ode`, `adiabatic_lapse_rate`, `potential_temperature`, `classify_stability`.
- **Derivations:** D18, D19 (★★★), D36.
- **Depth features:** Explain tab (N² with your lapse rate → N → period or e-folding → θ gradient), synced Code tab, + linked views (3), transport (release/play/scrub, end-of-run card), presets (dry adiabatic = neutral, standard atmosphere −6.5 K/km, nocturnal inversion, superadiabatic surface layer, ocean thermocline), status verdict ("🌊 stable: period 11 min" / "⚖️ neutral" / "🚀 unstable: e-folds in 3 min"), modes (atmosphere / ocean / lab tank), "Right now" notes.
- **Follows reference:** `forced_damped_vibrations.html` (system + x(t) graph + regime-dependent interpretation, "reading the current setting").
- **Aha:** air that gets colder with height can still be stable. What matters is whether it cools more slowly than a rising parcel cools by expanding, and that is exactly whether θ increases upward.

### E5 · buckingham_pi_machine
- **CORE:** C65, C66, C67, C68, C69, C70, C72 (C73, C75, C76 as presets)
- **Confusion removed:** the Π theorem as a magic recipe. The reader sees it as linear algebra: groups are null-space vectors of the dimensional matrix, their number is n − r, a different repeating set gives a different but equivalent basis, and forgetting a variable changes everything.
- **Why interactive:** toggling variables on and off and swapping the repeating set changes the matrix, its rank and the groups at once. The reader experiments with the method instead of watching it done once.
- **Stage:** (1) variable chips with their dimension vectors (on = column in the matrix); (2) the dimensional matrix with the chosen r × r repeating minor highlighted and its determinant shown (red when singular); (3) the resulting groups as monomials with exponents, plus a units-change check: switch SI → cgs → imperial and every Π value stays the same.
- **Controls:** problem preset (pipe Δp · pendulum period · sphere drag · Ex. 1.2 scale height · Ex. 1.4 blast · Ex. 1.5 Rayleigh) · variable toggles · repeating-set selector · unit system.
- **Equations shown:** (1.36), (1.37), (1.38), (1.39), (1.40), n − r, the exponent system of Π₁.
- **Mirrors:** `fluidpy.core.dimensional.dimensional_matrix`, `rank_by_minors`, `solve_exponents`, `pi_groups`, `groups_independent`.
- **Derivations:** D26, D27, D28 (★★★).
- **Depth features:** Explain tab (matrix → minors → r → n − r → each exponent solve with numbers → final dimension check), synced Code tab, + presets (the classic problems), inspector (click a matrix cell or a group: its exponent arithmetic), live status ("✅ 4 independent groups" / "⚠️ singular repeating set" / "⚠️ a variable is dimensionless on its own"), linked views (3).
- **Follows reference:** `amplitude_phase_second_order_II_3.html` (numbered live derivation), with presets and "formula with numbers plugged in" after `stride_padding_playground.html`.
- **Aha:** dimensionless groups are the null space of a small matrix. Buckingham's n − r is just rank–nullity, and every group's value survives a change of units.

### B1 · pressure_in_still_fluid (backup)
- **CORE:** C17, C19, C20, C21, C22 (C16 through the capillary mode)
- **Confusion removed:** gauge vs absolute pressure; "pressure pushes down only"; why a submerged block feels a net upward force; why water climbs a thin tube.
- **Why interactive:** drag a fluid cube or block through a layered tank and watch the face-pressure arrows: horizontal ones cancel, vertical ones differ by the weight, and in a different fluid the difference is buoyancy. Shrinking the tube radius makes the meniscus lift water higher.
- **Stage:** tank with 1–2 layers and a draggable cube/block with pressure arrows · p(z) (absolute/gauge) · mode "capillary tube" with meniscus, p along E–F and h vs R.
- **Controls:** upper-layer density · layer thickness · probe depth · block density · mode (tank / capillary) with R and α.
- **Equations shown:** (1.5), (1.7), (1.8), (1.9), Ex. 1.1, buoyancy (D37).
- **Mirrors:** `fluidpy.core.statics.hydrostatic_pressure_uniform`, `integrate_hydrostatic`, `gauge_pressure`, `buoyancy_force`; `fluidpy.ch01_introduction.capillary_rise`, `laplace_pressure_jump`.
- **Derivations:** D04, D05, D06, D37, D07 (would be added to those rows' "Shown in" only if this backup is built).
- **Depth features:** Explain, Code, + linked views, presets (water, oil over water, mercury manometer, glass capillary), inspector (click a face: p × area), term bars (top force, bottom force, weight, buoyancy), modes.
- **Follows reference:** `angular_frequency_explorer_1.html` (modes, presets, Right-now notes).
- **Aha:** buoyancy is nothing but the pressure difference between a body's bottom and top faces. The fluid's hydrostatic law does all the work.

## 6. Python animations and interactive figures (CORE ID → what, why, player/figure kind)
| Kind | CORE | What moves / what the slider controls | Why | Player / figure |
|---|---|---|---|---|
| Animation | C02 (+ N01, N02) | a solid and a fluid element under the same shear force; the load is removed half-way; γ(t) curves (with Bingham and Maxwell overlays) drawn so far | the definition is about what happens *over time* and after unloading | `show_animation(player="video")`, ≤ 60 frames |
| Animation | C06 | zoom-in: the sampling cube grows from 1 nm to 1 m over a gas with a density gradient; ρ estimate trace | the noise → plateau → drift sequence is a process of scale | `player="frames"` (step through decades) |
| Animation | C08, C12 (N05) | u(y, t) diffusing from a suddenly moved plate with τ(y) panel (FTCS, D = ν); the same run with D = κ for heat | diffusion is a transient; the settling to a linear profile must be seen | `player="video"` |
| Animation | C26 | piston following two paths between the same states with running q, w, Δe counters | path dependence seen as motion | `player="frames"` (leg by leg) |
| Animation | C50 | three parcels released at once in stable, neutral and unstable columns; ζ(t) traces | cos vs straight line vs cosh is a behaviour in time | `player="video"` |
| Animation | C75 | hemispherical blast front growing as t^{2/5} with a log–log D(t) trace | shows why a few photographs with times are enough | `player="video"`, 40 frames |
| Plotly slider | C07 (N03) | altitude (0–100 km): mean free path l and Kn vs body size with regime bands | the continuum's failure depends on two scales at once | `slider_figure`, 30 steps |
| Plotly slider | C13 | exponent n of μ ∝ T^n vs Sutherland for air (water curve fixed) | shows how rough the T^{1/2} law is | `slider_figure`, 20 steps |
| Plotly slider | C18 (N10) | wedge size dz: force triangle closing and p₂ − p₁ → 0 | the continuum limit in D03 made visible | `slider_figure`, 25 steps |
| Plotly slider | C21 | upper-layer density of a two-layer tank: p(z) kink and block face forces | how density sets the slope of p(z) | `slider_figure`, 25 steps |
| Plotly slider | C22 (N12) | wetting angle α: h vs tube radius R and p along E–F | 1/R law and the role of α | `slider_figure`, 20 steps |
| Plotly slider | C57 (N18) | inversion strength: T(z), θ(z) and N²(z) coloured by regime | reading stability from θ rather than T | `slider_figure`, 20 steps |
| Plotly slider | C62, C63 | temperature T: isothermal p(z) and H vs the USSA table | how good the isothermal model is | `slider_figure`, 20 steps |
| Plotly slider | C76 | particle volume V: S/I across 380–750 nm coloured by wavelength | λ⁻⁴ and why small particles scatter blue | `slider_figure`, 20 steps |
| Plotly 3-D | C16 (N09), C28 | sphere cap / saddle with principal circles (dropdown); p(v, T) surface with a state dot | curvature and "two properties fix the state" are 3-D ideas | `go.Surface`, height 520 |
| Live widget (kernel only, plus the slider figure for the page) | C51 | N² and ζ₀ with the nonlinear `parcel_ode` | free exploration beyond the precomputed slider | `live(...)` |

## 7. From-scratch moments
Every section with computable CORE items has at least one; each is followed by `assert np.allclose(mine, library)`.

| § | CORE | Hand-written version | Tested function it must agree with |
|---|---|---|---|
| 1.3 | C02 | loop over times: γ = τ/G (solid) and γ += τΔt/μ (fluid), unloaded half-way | `ch01.shear_deformation_history` |
| 1.4 | C05 | sum 2mu_x over molecules hitting a wall in Δt, divided by AΔt | `ch01.molecular_pressure` (and n k_B T) |
| 1.4 | C06 | count seeded random positions inside a cube, ×m/L³, repeat for the std | `ch01.sample_density` (same seed) |
| 1.5 | C08 | five-line FTCS loop f[1:-1] += r(f[2:] − 2f[1:-1] + f[:-2]) | `core.diffusion.ftcs_diffusion_1d` |
| 1.5 | C12 | central differences (u[i+1] − u[i−1])/(2Δy) with one-sided second-order ends, ×μ | `ch01.shear_stress_profile` |
| 1.6 | C16 | σ(1/R₁ + 1/R₂) with the sphere and flat special cases | `ch01.laplace_pressure_jump` |
| 1.7 | C20 | explicit Euler march p[k+1] = p[k] − ρ(z_k)gΔz on a fine grid | `core.statics.integrate_hydrostatic` (to O(Δz)) |
| 1.8 | C26 | trapezoid sum of −p dv along each leg, Δe = c_vΔT, q = Δe − w | `core.thermo.process_heat_work` |
| 1.8 | C30 | (h(T + ΔT, p) − h(T − ΔT, p))/(2ΔT) | `core.thermo.specific_heat_cp` |
| 1.9 | C40 | ρ = p/(k_B A_o/M_w · T) from the constants chain | `core.thermo.perfect_gas_density` |
| 1.9 | C45 | integrate dp/dρ = γp/ρ with small steps from (ρ₀, p₀) | `core.thermo.isentropic_pressure` |
| 1.10 | C50 | Euler–Cromer loop for ζ″ = −N²ζ | `core.stratification.parcel_displacement` |
| 1.10 | C55 | θ = T(p_ref/p)^{(γ−1)/γ} | `core.stratification.potential_temperature` |
| 1.11 | C68 | cofactor expansion of every 3×3 minor of (1.39) with `itertools.combinations` | `core.dimensional.rank_by_minors` and `np.linalg.matrix_rank` |
| 1.11 | C70 | `np.linalg.solve` on the 3×3 exponent system for Π₁ | `core.dimensional.solve_exponents` |

(§1.1 has no computable CORE item; §1.2 is recap only, but its R01 cell still converts units with `pint`.)

## 8. Notes for the implementer
Functions the notebook figures and explainers need that analysis §4 did not plan (all SI, docstrings per contract; mark
our own models "not in the book" and cite any external coefficient):
1. **Parity namespace.** `tools/shot.py` exposes `ch01.<fn>` (the chapter module) and `core` (the package, which does
   **not** import its submodules). Re-export every function an explainer mirrors from `fluidpy/ch01_introduction.py`
   (e.g. `from fluidpy.core.stratification import brunt_vaisala_sq, parcel_displacement, …`), so parity rows read
   `ch01.brunt_vaisala_sq(...)`.
2. **§1.3:** `ch01.traction_components(force_vec, normal_vec) -> (normal, shear)` (C03);
   `ch01.mean_molecular_spacing(rho, M_w)` and `number_density(rho, M_w)` (C04); extend `shear_deformation_history`
   with `kind="bingham"` (yield stress `tau_y`) and `kind="maxwell"` (N01).
3. **§1.4:** `ch01.density_noise_expected(L, number_density)` = (nL³)^{-1/2}, deterministic, for E1 parity and D35
   (the seeded `sample_density` cannot be matched by JS draws); `sample_density` needs the optional macroscopic gradient
   the E1 stage shows.
4. **§1.5:** `ch01.mass_fractions(mole_fractions, molar_masses)` with dry-air composition (cite USSA-1976);
   `ch01.thermal_diffusivity(k, rho, cp)`; `ch01.diffusion_time(L, D)`; `ch01.water_viscosity(T)` (cited fit, e.g.
   IAPWS 2008 or Vogel; benchmark row needed); `ch01.wall_shear_history(u_hist, dy, mu)` for E2's τ_w(t); optionally
   `ch01.couette_startup_profile(y, t, U, h, nu, nterms)` (series solution) as an analytic V1 check of the FTCS run.
   A Fourier-series primer is then needed only in `tests/`, not in the notebook.
5. **§1.7:** `core.statics.buoyancy_force(rho_fluid, volume, g)` and `net_pressure_force_on_box(p_fn, box)`
   (face-by-face integration, D37); `ch01.capillary_rise_deg` as planned.
6. **§1.8:** `ch01.mean_molecular_speed(T, m)` = √(8k_BT/πm), `ch01.collision_time(l, T, m)` (C23);
   `core.thermo.process_path(kind, state1, state2, n)` generating isothermal / isochoric / isobaric / isentropic legs
   (C26, E3); an irreversible-stirring path helper `stirred_isochoric_path` (C32, C34, E3 irreversible mode);
   `core.thermo.tait_pressure(rho, rho0, K0, n)` stiff-liquid EOS for C36 (cite); `ch01.water_density(T)` cited fit
   (e.g. Kell 1975 or IAPWS) for C37's 4 °C maximum.
7. **§1.9:** `core.thermo.van_der_waals_pressure(T, v, a, b, R)` and `van_der_waals_internal_energy(T, v, a, cv)` for
   C41 and the D13 sympy check (cite a, b for CO₂); `core.thermo` table of M_w for He, H₂O, CO₂ (C39) and γ for mono-,
   di-, triatomic gases (C43).
8. **§1.10:** `core.stratification.parcel_temperature(T0, z, cp, g)` (dry adiabat through a point);
   `brunt_vaisala_sq_from_theta(theta, dtheta_dz, g)` and `brunt_vaisala_sq_from_lapse(T, dT_dz, cp, g)` (D36, E4);
   `stability_timescale(N2) -> (kind, period_or_efold)` (E4 status); `ch01.synthetic_boundary_layer_profile(z,
   inversion_strength)` (our profile for N18/C57); `ch01.seawater_density_linear(T, S, rho0, alpha_T, beta_S, T0, S0)`
   with cited coefficients (e.g. Vallis 2017, TEOS-10 linearisation) for C60 and E4 ocean mode.
9. **§1.11:** `core.dimensional.group_value(group, values)` and `rescale_units(values, system)` (C69 units-invariance
   figure, E5 unit switch); preset dictionaries `PIPE`, `PENDULUM`, `SPHERE_DRAG`, `SCALE_HEIGHT`, `BLAST`, `RAYLEIGH`;
   `ch01.poiseuille_pressure_drop(mu, U, dx, d)` as a data generator for the C65/C72 collapse figures (label it
   "result derived in Ch. 8"); `ch01.pythagoras_phi(beta)` = ¼ sin 2β (C74); `ch01.wavelength_to_rgb(lam)` for the C76
   spectrum colours (visual only). E5's JS needs exact rational arithmetic (a small Fraction class) for the null-space
   solve. It is a promotion candidate for `viz_lib.js` once a second chapter needs it (Ch. 4 §4.11).
10. **§1.1:** C01's chapter map reads `book.yaml` through `fluidpy.core.project` and draws it with matplotlib. No
    physics function is needed, but put the drawing in `scripts/ch01_book_map.py` so the notebook cell stays short.
11. **Budget (< 5 min on Colab CPU):** cache the FTCS runs (C08, N05, E2 presets are JS-side) and the `sample_density`
    sweeps with `FAST` grids (N = 200 in FAST, 400 full; ≤ 40 slider steps); `parcel_ode` capped with a `solve_ivp`
    event for the unstable case; sympy checks (D13, D19, D28) are sub-second but must avoid `simplify` on large
    expressions (use `expand` + `cancel`).
12. **Teaching traps to carry into code names and docstrings:** lapse rates follow the book (user decision):
    `adiabatic_lapse_rate(...)` returns Γ_a = dT/dz in K/m, **negative** for air (≈ −9.8e-3), and any environment
    lapse-rate argument is named `dT_dz` (Γ ≡ dT/dz, negative when T falls with height). Do not add a positive
    "meteorological" variant. Tests compare against the published magnitude with `abs()`, and docstrings get at most one
    sentence noting that some texts quote the magnitude as positive;
    `capillary_rise` takes the book's α (complement of the contact angle); Clausius–Duhem is taught with the actual δq;
    (1.35) holds as "same sign" (analysis §9 typo 7); reference pressure `p_ref` ≠ surface pressure `p0`; kmol vs mol
    in every gas constant.
