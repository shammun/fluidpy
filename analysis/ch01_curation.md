# Chapter 1 — Introduction: curation
(from `analysis/ch01.md` — 106 inventory rows, 40 numbered equations, 33 derivations + 4 added; `book.yaml` ch01. First
chapter: nothing to RECAP from earlier chapters, no primer to reuse. First curation 2026-09-12; **re-curated for depth
2026-09-12** on the user's instruction, concept-curator.)

**What changed in the re-curation (user decision, binding).** The first curation made all 76 new ideas CORE at full
depth, which was too slow and gave every idea the same volume. Coverage stays exhaustive (all 106 rows keep their IDs
and appear once in the chapter map, §2), but each row now has a **depth**:

| Depth | What the reader gets | Tier word (parsed by `tools/nbkit.py`) |
|---|---|---|
| **A · full treatment** | picture → question → derivation step by step → worked number → code (tested function + from-scratch) → figure → how to read it | CORE |
| **B · stated and explained** | a paragraph, the equation, a number where it helps, written **inside its A parent's block**; no separate derivation (a demoted derivation's result is given, §4c) | NOTE (RECAP stays RECAP) |
| **C · named** | one sentence with a pointer to where the idea is used later | NOTE (RECAP stays RECAP, SKIP stays SKIP) |

**Placement rule for B and C items whose book section has no A item** (§1.1, §1.2, §1.3, §1.6): every book section still
gets its own notebook section. There the B/C items are written as short paragraphs tagged with their A parent, saying
which A block develops or uses them (e.g. §1.6's Laplace jump is stated in §1.6 and used by the C20 block's capillary
rise). Where the parent comes earlier, the B paragraph sits inside the parent's block, in the parent's section.

**Lapse-rate convention (user decision, binding).** Code computes with Kundu's **Γ ≡ dT/dz** (Γ_a = −g/C_p ≈ −9.8 K/km;
stable when Γ > Γ_a). The standard meteorology convention **Γ ≡ −dT/dz** (Γ_a ≈ +9.8 K/km; stable when Γ < Γ_a) is always
shown next to it, in a two-row convention table in the C54 block, in every lapse-rate readout of E4 and in the E4
convention panel. It is the standard in atmospheric science and the lapse-rate-feedback literature Shammunul reads, so
it is taught as a first-class form, and conversion (Γ_met = −Γ_Kundu, flip the inequality) is practised with numbers.

Other decisions kept from the first curation: four derivations the book never writes are kept (D34 kinetic pressure,
D35 continuum noise law, D36 N² from θ, D37 Archimedes from hydrostatics) and re-parented to A items. The analyst's
suspected typos (analysis §9) are taught as written there: Clausius–Duhem with the actual δq, (1.35) as "has the same
sign as".

**Totals:** 106 rows = **A 15 · B 70 · C 21**. Tier words: **CORE 15 · RECAP 3 · NOTE 86 · SKIP 2**. Derivations: **12
kept** (★ 6 · ★★ 4 · ★★★ 2) · **25 demoted** to stated results (§4c). Explainers 5 + 1 backup.

## 1. Teaching order (CORE IDs grouped by book section; one sentence each: "once you see X, Y follows")
A items are in **bold**; B and C items are listed where they are taught, with their A parent.

**§1.1 Fluid Mechanics**
- C01 [B → C06] What fluid mechanics is and its three routes (analysis, computation, experiment): once you see the map of the subject, you can see where every later idea in this chapter will be used.

**§1.2 Units of Measurement** (R01 [B → C64], R02 [C → C64], R03 [B → C40]: SI units, prefixes and kelvin with `pint`, so every later number carries its unit)

**§1.3 Solids, Liquids, and Gases**
- C02 [B → C12] Fluid vs solid: a fluid keeps deforming under *any* shear stress while a solid settles at a fixed strain (the rate law behind it is C12). N01 [C], N02 [C].
- C03 [B → C20] Normal stress: fluids resist being squeezed, but liquids break (cavitate) when pulled; the squeeze is the pressure of C20.
- C04 [B → C06] Liquid vs gas: molecules in a gas are about 10× farther apart than in a liquid, which is the spacing C06 averages over.

**§1.4 Continuum Hypothesis**
- C05 [B → C06] Pressure as the average of molecular impacts (D34 derived inside C06): once pressure is momentum delivered by many molecules, you can ask how many molecules an average needs.
- **C06 The continuum hypothesis** (density at a point as an average over a window of box sizes): once averages settle to a plateau as the box grows, the idea of a "value at a point" follows.
- C07 [B → C06] Knudsen number Kn = l/L (N03 [B] the value of l): once the continuum needs boxes much larger than molecular scales, l/L is the test.
- C24 [B → C06] Fluid particle (taught here with C06, although the book names it in §1.8): big enough to average, small enough to be a point, relaxes fast (C23).

**§1.5 Molecular Transport Phenomena**
- C08 [B → C12] Random molecular motion diffuses properties down gradients. C09 [C], C10 [B] Fick (1.1), N04 [C], C11 [B] Fourier (1.2).
- **C12 Newton's law of viscosity (1.3)**: the diffusion picture applied to momentum gives shear stress = μ du/dy, and a stress is a momentum flux.
- N05 [B], C13 [B] μ(T), C14 [B] ν = μ/ρ (1.4), N06 [C]: inside C12, once momentum diffuses, its diffusivity ν (m²/s), not μ, sets how fast motion spreads.

**§1.6 Surface Tension**
- C15 [B → C20] σ as force per length / energy per area; N07 [B] sphere Δp = 2σ/R; C16 [B] Laplace jump (1.5); N08 [C]; N09 [C]: a curved interface carries a pressure difference, which C20's block uses for capillary rise.

**§1.7 Fluid Statics**
- C17 [B → C20] absolute vs gauge; C18 [B] isotropy (1.6); N10 [C]; C19 [B] Pascal (1.7): once pressure is a scalar and gravity is vertical, only a vertical gradient can survive.
- **C20 Hydrostatic law dp/dz = −ρg (1.8)**: the vertical balance on a fluid cube gives the law everything in §1.10 builds on; D05 derives it, D37 turns it into Archimedes' buoyancy.
- C21 [B] uniform density (1.9) + buoyancy; N11 [B]; C22 [B] capillary rise (Ex. 1.1); N12 [C]; C03 [B] from §1.3.

**§1.8 Classical Thermodynamics**
- C23 [B → C25] system, equilibrium, relaxation time.
- **C25 First law δq + δw = Δe (1.10)**: once a fluid particle has an energy, heat in plus work done on it must change that energy, and drawing two paths between the same states shows q and w depend on the path while Δe does not.
- N13 [C], C26 [B] path vs state functions, C27 [B] de = dq − p dv (1.11), C28 [B] equations of state (1.12): inside C25.
- C29 [B → C35] h = e + pv (1.13), C30 [B] C_p (1.14), C31 [B] C_v (1.15), C32 [C], C33 [B] entropy (1.16), N15 [B] T ds = dq (1.17), C34 [B] Clausius–Duhem, N14 [C].
- **C35 Gibbs relations (1.18)**: once T ds = dq and de = dq − p dv hold, eliminating dq gives relations between state functions that hold for *any* process.
- **C36 Speed of sound c² = (∂p/∂ρ)_s (1.19)**: how stiff the equation of state is at constant entropy sets the speed of pressure signals; infinite stiffness is the incompressible limit.
- C37 [B → C36] α (1.20): the same kind of partial derivative at constant pressure (used again in C54).

**§1.9 Perfect Gas**
- C38 [B → C40] pV = nk_BT (1.21), C39 [B] constants.
- **C40 Continuum perfect-gas law p = ρRT (1.22)**: once n/V × m is the continuum density, (1.21) becomes the law used in the rest of the book.
- C41 [B] e = e(T), C48 [B] α = 1/T (1.28), R03 [B] kelvin, and (from §1.10) C62 [B] isothermal atmosphere, C63 [B] scale height H = RT/g: inside C40, because each is p = ρRT combined with one earlier idea.
- C42 [B → C45] R = C_p − C_v (1.23), C43 [B] γ (1.24), N16 [B], C44 [B] adiabatic vs isentropic.
- **C45 Isentropic perfect gas p/ρ^γ = const (1.25)**: once ds = 0 in the Gibbs relations, separating variables gives this power law.
- C46 [B] isentropic ratios (1.26), C47 [B] c = √(γRT) (1.27): inside C45, one line each from the power law (C47 closes the loop on C36).

**§1.10 Stability of Stratified Fluid Media**
- C49 [B → C54] static medium: hydrostatics + equation of state build p, ρ, T from one profile (the environment of E4).
- **C50 The displaced-parcel equation**: once a parcel keeps its own density while its surroundings change, Newton's law with buoyancy (D37) gives ζ'' + N²ζ = 0 (D18).
- **C51 Brunt–Väisälä frequency N² (1.29)**: once the parcel equation is linear, its coefficient is the squared frequency of vertical oscillation, and a negative N² turns cos into cosh.
- C52 [B] stable/neutral/unstable, C60 [B] seawater S, N21 [B], C61 [B] ocean criterion (1.35): inside C51.
- C53 [B → C54] Lapse rate: Kundu's Γ ≡ dT/dz (standard troposphere −6.5 K/km) and the meteorology convention Γ ≡ −dT/dz (+6.5 K/km), both shown.
- **C54 Adiabatic lapse rate Γ_a = −gαT/C_p (1.30)**: once a parcel rises isentropically through hydrostatic surroundings, its temperature changes at Γ_a = −g/C_p ≈ −9.8 K/km in Kundu's convention (= +9.8 K/km as the meteorology lapse rate): stable when dT/dz > Γ_a, i.e. when −dT/dz < +9.8 K/km. N17 [B] the value, in both conventions.
- **C55 Potential temperature θ (1.31)**: once adiabatic cooling is predictable, bringing a parcel to a reference pressure removes it and leaves θ; D36 shows N² = (g/θ)dθ/dz, so stability is θ increasing upward.
- C56 [B] (1.32), N19 [C], C57 [B] stability from dθ/dz, N18 [B], C58 [B] potential density (1.33), C59 [C] (1.34), N20 [C]: inside C55.

**§1.11 Dimensional Analysis**
- **C64 Dimensional homogeneity**: once laws cannot depend on the units we choose, every term in a correct equation must have the same dimensions (R01, R02 recapped here).
- C65 [B → C67] choose the variables (1.38), C66 [B] [q] and M, L, T, θ.
- **C67 The dimensional matrix (1.39)**: stacking dimension vectors as columns gives one matrix holding all the information; its rank r (C68 [B]) counts the independent dimensions.
- **C69 Buckingham's Π theorem (1.37)**: once dimensionless groups are the null space of the matrix, rank–nullity gives n − r of them (D28). (Taught after C67, not in book order, because its proof needs the matrix and its rank.)
- N22 [C], N23 [B] n − r, C70 [B] exponent algebra, N24 [C], C71 [B] combining groups, C72 [B] pipe law (1.40), N25 [B], C73 [B] Ex. 1.2, C74 [C] Ex. 1.3, C75 [B] Ex. 1.4, C76 [B] Ex. 1.5: inside C69.
- S01 [C], S02 [C]: pointer lines at the end.

## 2. Chapter map (depth)
One row per inventory row (106). `A parent` names the A block a B or C item lives in; `pointer` is for C items (where
the idea is used later, or the pointer line for SKIP).

| ID | Item | § | Depth | Tier | A parent | Pointer |
|---|---|---|---|---|---|---|
| C01 | What fluid mechanics is; the three routes (analysis, computation, experiment) [#1] | 1.1 | B | NOTE | C06 | — |
| R01 | SI base quantities and derived units N, Pa, J, W, Hz [#2] | 1.2 | B | RECAP | C64 | — |
| R02 | SI prefixes [#3] | 1.2 | C | RECAP | C64 | pint handles prefixes wherever a number is printed: hPa in §1.10, mN/m in §1.6, nm in §1.4 |
| R03 | Celsius–kelvin relation [#4] | 1.2 | B | RECAP | C40 | — |
| C02 | Fluid vs solid: a fluid deforms without limit under any shear stress [#5] | 1.3 | B | NOTE | C12 | — |
| N01 | Plastic (yield) and viscoelastic (partial memory) materials [#6] | 1.3 | C | NOTE | C12 | non-Newtonian behaviour: limits of the Newtonian constitutive law in Ch. 4 §4.5; blood in Ch. 16 §16.3 |
| C03 | Normal stress: compression and tension; liquids cavitate under tension [#7] | 1.3 | B | NOTE | C20 | — |
| C04 | Liquid vs gas (molecular spacing); free surface [#8] | 1.3 | B | NOTE | C06 | — |
| N02 | Fig. 1.1: solid element deflects to a fixed shape, fluid keeps deforming [#9] | 1.3 | C | NOTE | C12 | drawn (our drawing) as the γ(t) sketch beside the C02 statement in §1.3 |
| C05 | Pressure as the statistical average of molecular collision force per area [#10] | 1.4 | B | NOTE | C06 | — |
| C06 | Continuum hypothesis: properties at a point as averages over δV [#11] | 1.4 | A | CORE | — | — |
| C07 | Knudsen number Kn = l/L [#12] | 1.4 | B | NOTE | C06 | — |
| N03 | Mean free path of air at room conditions; water's is much smaller [#13] | 1.4 | B | NOTE | C06 | — |
| C08 | Random molecular motion diffuses species, heat and momentum down gradients [#14] | 1.5 | B | NOTE | C12 | — |
| C09 | Mass fraction Y and partial density ρY [#15] | 1.5 | C | NOTE | C12 | salinity S (C60) is a mass fraction; salt and heat diffusing at different rates drive double-diffusive instability, Ch. 11 §11.5 |
| C10 | Fick's law J_m = −ρκ_m∇Y (1.1) [#16] | 1.5 | B | NOTE | C12 | — |
| N04 | Fig. 1.2: Y increases with y, flux across AB points down [#17] | 1.5 | C | NOTE | C12 | the flux-arrow sketch is the "species" mode of explainer E2 |
| C11 | Fourier's law q = −k∇T (1.2) [#18] | 1.5 | B | NOTE | C12 | — |
| C12 | Newton's law of friction τ = μ du/dy (1.3); dynamic viscosity μ [#19] | 1.5 | A | CORE | — | — |
| N05 | Fig. 1.3: τ on AB, u(y) relaxes toward the dashed profile [#20] | 1.5 | B | NOTE | C12 | — |
| C13 | μ of a gas rises (≈ T^{1/2}), μ of a liquid falls with T [#21] | 1.5 | B | NOTE | C12 | — |
| C14 | Kinematic viscosity ν = μ/ρ (1.4) [#22] | 1.5 | B | NOTE | C12 | — |
| N06 | The transport laws are linear and contain first derivatives only [#23] | 1.5 | C | NOTE | C12 | generalised to a tensor law in Ch. 4 §4.5; the diffusion equations they lead to appear in Ch. 4 §4.8 (heat) and Ch. 8 §8.4 (momentum) |
| C15 | Surface tension σ: force per length = energy per area [#24] | 1.6 | B | NOTE | C20 | — |
| N07 | Spherical interface: Δp = 2σ/R [#25] | 1.6 | B | NOTE | C20 | — |
| C16 | Laplace pressure jump Δp = σ(1/R₁ + 1/R₂) (1.5) [#26] | 1.6 | B | NOTE | C20 | — |
| N08 | Capillary tube, capillarity [#27] | 1.6 | C | NOTE | C20 | surface tension as a boundary condition in Ch. 4 §4.10; capillary waves in Ch. 7 §7.3 |
| N09 | Fig. 1.4: hemispherical drop forces; patch with two principal radii [#28] | 1.6 | C | NOTE | C20 | principal radii return with capillary waves, Ch. 7 §7.3; one small sketch beside the C16 statement |
| C17 | Static pressure; absolute vs gauge; standard atmosphere; bar [#29] | 1.7 | B | NOTE | C20 | — |
| C18 | Pressure at a point is the same in all directions (1.6) [#30] | 1.7 | B | NOTE | C20 | — |
| N10 | Fig. 1.5: triangular element and pressure forces [#31] | 1.7 | C | NOTE | C20 | the isotropic stress becomes the −pδ_ij part of the stress tensor, Ch. 2 §2.6 and Ch. 4 §4.5 |
| C19 | Pascal's law ∂p/∂x = ∂p/∂y = 0 at rest (1.7) [#32] | 1.7 | B | NOTE | C20 | — |
| C20 | Hydrostatic law dp/dz = −ρg (1.8) [#33] | 1.7 | A | CORE | — | — |
| C21 | Uniform density p = p₀ − ρgz (1.9); pressure rise ρgh; buoyancy from pressure [#34] | 1.7 | B | NOTE | C20 | — |
| N11 | Fig. 1.6: fluid cube, top/bottom pressure difference balances weight [#35] | 1.7 | B | NOTE | C20 | — |
| C22 | Capillary rise h = 2σ sin α/(ρgR) (Ex. 1.1) [#36] | 1.7 | B | NOTE | C20 | — |
| N12 | Fig. 1.7: meniscus, pressure along E–F, force balance on ABCD [#37] | 1.7 | C | NOTE | C20 | the p(z) along E–F sketch beside the C22 statement; free-surface pressure jumps return in Ch. 7 §7.3 |
| C23 | Thermodynamic system, equilibrium, relaxation time [#38] | 1.8 | B | NOTE | C25 | — |
| C24 | Fluid particle: fixed molecules, big enough to average, relaxes fast [#39] | 1.8 | B | NOTE | C06 | — |
| C25 | First law δq + δw = Δe (1.10) [#40] | 1.8 | A | CORE | — | — |
| C26 | Heat and work are path functions, e is a state function; reversible process [#41] | 1.8 | B | NOTE | C25 | — |
| N13 | Specific volume v = 1/ρ [#42] | 1.8 | C | NOTE | C25 | used in (1.11)–(1.18) and in the energy equation, Ch. 4 §4.8 |
| C27 | Reversible first law de = dq − p dv (1.11) [#43] | 1.8 | B | NOTE | C25 | — |
| C28 | Equations of state p = p(v, T), e = e(p, T) (1.12); two properties fix the state [#44] | 1.8 | B | NOTE | C25 | — |
| C29 | Enthalpy h = e + pv (1.13) [#45] | 1.8 | B | NOTE | C35 | — |
| C30 | C_p = (∂h/∂T)_p (1.14) [#46] | 1.8 | B | NOTE | C35 | — |
| C31 | C_v = (∂e/∂T)_v (1.15) [#47] | 1.8 | B | NOTE | C35 | — |
| C32 | Heat per degree = C_v (constant v) or C_p (constant p) for reversible p dv work [#48] | 1.8 | C | NOTE | C35 | why C_p dT is the heat added at constant pressure; used when h is written C_pT in Ch. 4 §4.8 and in stagnation enthalpy, Ch. 15 §15.4 |
| C33 | Second law (i): entropy s₂ − s₁ = ∫dq_rev/T (1.16) [#49] | 1.8 | B | NOTE | C35 | — |
| C34 | Second law (ii): Clausius–Duhem inequality, s₂ − s₁ ≥ ∫δq/T [#50] | 1.8 | B | NOTE | C35 | — |
| N14 | Second law (iii): μ > 0 and k > 0 [#51] | 1.8 | C | NOTE | C35 | viscous dissipation is never negative, Ch. 4 §4.8 |
| N15 | T ds = dq for a reversible process (1.17) [#52] | 1.8 | B | NOTE | C35 | — |
| C35 | Gibbs relations T ds = de + p dv = dh − v dp (1.18) [#53] | 1.8 | A | CORE | — | — |
| C36 | Speed of sound c² = (∂p/∂ρ)_s (1.19); incompressible limit c → ∞ [#54] | 1.8 | A | CORE | — | — |
| C37 | Thermal expansion coefficient α = −(1/ρ)(∂ρ/∂T)_p (1.20) [#55] | 1.8 | B | NOTE | C36 | — |
| C38 | Molecular perfect-gas law pV = nk_BT (1.21) [#56] | 1.9 | B | NOTE | C40 | — |
| C39 | k_B, Avogadro per kmol, R_u = k_BA_o, M_w, R = R_u/M_w [#57] | 1.9 | B | NOTE | C40 | — |
| C40 | Continuum perfect-gas law p = ρRT (1.22) [#58] | 1.9 | A | CORE | — | — |
| C42 | R = C_p − C_v (1.23) [#59] | 1.9 | B | NOTE | C45 | — |
| C43 | γ = C_p/C_v (1.24) [#60] | 1.9 | B | NOTE | C45 | — |
| N16 | Air values of γ and C_p; C_p and C_v rise with T [#61] | 1.9 | B | NOTE | C45 | — |
| C41 | A perfect gas has e = e(T), h = h(T) (and conversely) [#62] | 1.9 | B | NOTE | C40 | — |
| C44 | Adiabatic (no heat) vs isentropic (adiabatic + frictionless) [#63] | 1.9 | B | NOTE | C45 | — |
| C45 | Isentropic perfect gas p/ρ^γ = const (1.25) [#64] | 1.9 | A | CORE | — | — |
| C46 | Isentropic ratios T/T₀ = (p/p₀)^{(γ−1)/γ}, ρ/ρ₀ = (p/p₀)^{1/γ} (1.26) [#65] | 1.9 | B | NOTE | C45 | — |
| C47 | Speed of sound in a perfect gas c = √(γRT) (1.27) [#66] | 1.9 | B | NOTE | C45 | — |
| C48 | α = 1/T for a perfect gas (1.28) [#67] | 1.9 | B | NOTE | C40 | — |
| C49 | Static medium: (1.8) with (1.12) tie p, ρ, T; one profile fixes the others [#68] | 1.10 | B | NOTE | C54 | — |
| C50 | Displaced-parcel equation ζ'' − (g/ρ)(dρ/dz − dρ_a/dz)ζ = 0 [#69] | 1.10 | A | CORE | — | — |
| C51 | Brunt–Väisälä frequency N² (1.29) [#70] | 1.10 | A | CORE | — | — |
| C52 | Stable (N² > 0), neutral (= 0), unstable (< 0) [#71] | 1.10 | B | NOTE | C51 | — |
| C53 | Lapse rate: Kundu Γ ≡ dT/dz and the meteorology convention Γ ≡ −dT/dz [#72] | 1.10 | B | NOTE | C54 | — |
| C54 | Adiabatic lapse rate Γ_a = −gαT/C_p (1.30) [#73] | 1.10 | A | CORE | — | — |
| N17 | Value of Γ_a for Earth's atmosphere: −9.8 K/km (Kundu) = +9.8 K/km (meteorology) [#74] | 1.10 | B | NOTE | C54 | — |
| N18 | Fig. 1.9a: T(z) with near-neutral layer, inversion, stable layer, neutral reference lines [#75] | 1.10 | B | NOTE | C55 | — |
| C55 | Potential temperature θ (1.31) [#76] | 1.10 | A | CORE | — | — |
| N19 | Log-derivative of (1.31) [#77] | 1.10 | C | NOTE | C55 | a step of D36 in the C55 block; the same move gives the atmosphere's θ profiles in Ch. 13 §13.2 |
| C56 | (T/θ)dθ/dz = dT/dz + g/C_p = Γ − Γ_a (1.32) [#78] | 1.10 | B | NOTE | C55 | — |
| C57 | Stability from the sign of dθ/dz; lab-scale T ≈ θ; N² = (g/θ)dθ/dz [#79] | 1.10 | B | NOTE | C55 | — |
| C58 | Potential density ρ_θ (1.33) [#80] | 1.10 | B | NOTE | C55 | — |
| N20 | θρ_θ = p_o/R = const [#81] | 1.10 | C | NOTE | C55 | the one line of algebra behind (1.34); not used again |
| C59 | −(1/ρ_θ)dρ_θ/dz = (1/θ)dθ/dz (1.34) [#82] | 1.10 | C | NOTE | C55 | density-based stability for the ocean in C61; vertical density variation in Ch. 13 §13.2 |
| C60 | Seawater: salinity S, ρ = ρ(T, p, S), potential density at constant S [#83] | 1.10 | B | NOTE | C51 | — |
| N21 | Parcel density gradient dρ_a/dz = −ρg/c² in the ocean [#84] | 1.10 | B | NOTE | C51 | — |
| C61 | Ocean stability criterion dρ_θ/dz ≅ dρ/dz + ρg/c² (1.35) [#85] | 1.10 | B | NOTE | C51 | — |
| C62 | Isothermal atmosphere p = p₀e^{−gz/RT} [#86] | 1.10 | B | NOTE | C40 | — |
| C63 | Scale height H = RT/g [#87] | 1.10 | B | NOTE | C40 | — |
| C64 | Dimensional homogeneity; laws can be written dimensionless [#88] | 1.11 | A | CORE | — | — |
| N22 | The relation f(q₁ … q_n) = 0 (1.36) [#89] | 1.11 | C | NOTE | C69 | written once as the input dictionary of `pi_groups`; Ch. 4 §4.11 starts from the same form |
| C69 | Buckingham's Π theorem φ(Π₁ … Π_{n−r}) = 0 (1.37) [#90] | 1.11 | A | CORE | — | — |
| C65 | Step 1: choose variables; pipe problem f(Δp, Δx, d, ε, U, ρ, μ) = 0 (1.38) [#91] | 1.11 | B | NOTE | C67 | — |
| C66 | Step 2: [q] notation, base dimensions M, L, T, θ [#92] | 1.11 | B | NOTE | C67 | — |
| C67 | Dimensional matrix of the pipe problem (1.39) [#93] | 1.11 | A | CORE | — | — |
| C68 | Step 3: rank as the largest nonzero minor [#94] | 1.11 | B | NOTE | C67 | — |
| N23 | Step 4: number of groups is n − r [#95] | 1.11 | B | NOTE | C69 | — |
| C70 | Step 5: exponent algebra with a repeating set [#96] | 1.11 | B | NOTE | C69 | — |
| N24 | Step 5 by inspection: remove M, L, T one at a time [#97] | 1.11 | C | NOTE | C69 | the ratio-by-ratio move is how Navier–Stokes is made dimensionless in Ch. 4 §4.11 |
| C71 | Combining groups: products and powers are groups; only n − r independent [#98] | 1.11 | B | NOTE | C69 | — |
| C72 | Step 6: pipe law Δp/ρU² = φ(Δx/d, ε/d, μ/ρUd) (1.40) [#99] | 1.11 | B | NOTE | C69 | — |
| N25 | Step 7: use physics (proportionality, linearity) to simplify [#100] | 1.11 | B | NOTE | C69 | — |
| C73 | Ex. 1.2: scale height by dimensional analysis [#101] | 1.11 | B | NOTE | C69 | — |
| C74 | Ex. 1.3: Pythagoras by dimensional analysis [#102] | 1.11 | C | NOTE | C69 | a curiosity: the same "a single group must be a constant" move as C73; not used later in the book |
| C75 | Ex. 1.4: G. I. Taylor's blast energy E = KρD⁵/t² [#103] | 1.11 | B | NOTE | C69 | — |
| C76 | Ex. 1.5: Rayleigh scattering S/I = V²φ₃(n_s)/(d²λ⁴) [#104] | 1.11 | B | NOTE | C69 | — |
| S01 | Exercises 1.1–1.30 [#105] | Ex | C | SKIP | — | "The chapter's 30 exercises are in the book and are not reproduced here. Of the derivations the book leaves to Exercises 1.10, 1.11, 1.13 and 1.14, D14 (isentropic law), D18 (parcel equation) and D19 (adiabatic lapse rate) are written out above; the result of Exercise 1.10 (e = e(T) for a perfect gas) is stated in the C40 block." |
| S02 | Literature cited and supplemental reading [#106] | Ref | C | SKIP | — | "The book's reading list is at the end of the printed chapter; the public data we used (CODATA, USSA-1976, IAPWS, Jennings, Taylor) are cited in reference/ch01/SOURCES.md." |

## 2a. The A list (15) — why each gets the full treatment
- **C06** continuum hypothesis — every field variable in the book (ρ, p, **u**, T at a point) rests on it; ★LB.
- **C12** Newton's law of viscosity — becomes the Newtonian constitutive law (Ch. 4 §4.5) behind Navier–Stokes, Ch. 8, Ch. 9; carries ν and the diffusion picture; ★LB.
- **C20** hydrostatic law — the base state of every stratified, wave and GFD problem (Ch. 4, 7, 11, 13) and the parent of buoyancy; ★LB.
- **C25** first law — becomes the energy equation (Ch. 4 §4.8); the path/state distinction it carries is the classic stumbling block of all thermodynamics that follows.
- **C35** Gibbs relations — the one tool from which isentropic laws, the adiabatic lapse rate and θ are derived; ★LB.
- **C36** speed of sound — sets when a flow may be treated as incompressible (Ch. 4 §4.9) and is the start of Ch. 15; ★LB.
- **C40** p = ρRT — the equation of state closing every atmospheric calculation in the book, and with (1.8) it gives the scale height; ★LB.
- **C45** isentropic law p/ρ^γ = const — the bridge from Gibbs to (1.26), θ and c = √(γRT); its derivation (D14) is the first long thermodynamic one.
- **C50** displaced-parcel equation — the first stability analysis of the book (a normal-mode argument in miniature, Ch. 11 §11.2); D18 is where N² comes from.
- **C51** N² (1.29) — internal waves (Ch. 7 §7.8, Ch. 13 §13.14), the Richardson number (Ch. 11 §11.7) and all of GFD; ★LB.
- **C54** adiabatic lapse rate (1.30) — the reference against which atmospheric stability is judged, and the place where the two sign conventions must be mastered; ★★★ derivation D19.
- **C55** potential temperature (1.31) — the conserved temperature of dry dynamics, used throughout Ch. 13 and in climate work; ★LB.
- **C64** dimensional homogeneity — every test and non-dimensional equation in the project depends on it; ★LB.
- **C67** dimensional matrix (1.39) — turns dimensional analysis into linear algebra (rank, minors, exponent solves), which is what makes the Π theorem computable.
- **C69** Buckingham's Π theorem (1.37) — the basis of dynamic similarity (Ch. 4 §4.11) and every dimensionless group; ★LB; ★★★ derivation D28.

## 2b. What each A block contains (for the lesson-designer)
Every A block: picture → question → derivation(s) → worked number → code (tested fluidpy function; from-scratch
version where §7 lists one) → figure → What you see / How to read it / What would change if. B items are paragraphs
with the equation and a number; C items are one sentence with the pointer from §2.

| A | Picture and question | Derivations (§4b) | Worked number | Code and figure | B statements inside (result stated; demoted D in §4c) | C sentences |
|---|---|---|---|---|---|---|
| **C06** | molecules in a box that grows: when does "density at a point" make sense? | D34 kinetic pressure, D35 noise ∝ N^{-1/2} | 10 µm air cube holds ≈ 2.5×10¹⁰ molecules, noise ≈ 6×10⁻⁶ | `ch01.sample_density` vs from-scratch count; log–log noise → plateau → drift figure; zoom animation; Kn slider; explainer E1 | C01 (chapter map from `book.yaml`, §1.1), C04 spacing (water 0.31 nm, air 3.4 nm), C05 pressure as impacts, C07 Kn = l/L, N03 l ≈ 67 nm, C24 fluid particle | — |
| **C12** | two plates, one suddenly moved: how does the fluid in between learn about it? | none (the book derives nothing in §1.5) | Couette gap h = 1 mm, U = 1 m/s, water μ = 1.0×10⁻³ Pa s → τ = 1 Pa; h²/ν = 1 s (water) vs ≈ 0.07 s (air) | `ch01.shear_stress_profile` vs central differences; u(y, t) diffusion animation with τ(y); explainer E2 | C02 fluid vs solid (γ = τt/μ), C08 diffusion picture, C10 Fick (1.1), C11 Fourier (1.2), N05 relaxation figure, C13 μ(T) (gas ↑, liquid ↓), C14 ν = μ/ρ (1.4) with air ≈ 15× water | N01, N02, C09, N04, N06 |
| **C20** | a fluid cube at rest: what must the pressure do for it not to move? | D05 hydrostatic law, D37 buoyancy = ρgV | 10 m of water ≈ 0.98×10⁵ Pa ≈ 1 atm | `core.statics.integrate_hydrostatic` vs Euler march; layered-tank p(z) figure; two-layer slider with block face forces | C03 normal stress/cavitation, C15 σ, N07 2σ/R, C16 Laplace (1.5), C17 gauge/absolute, C18 isotropy (1.6), C19 Pascal (1.7), C21 (1.9) + buoyancy, N11 cube, C22 capillary rise (tube radius 1 mm, fully wetting → h ≈ 15 mm) | N08, N09, N10, N12 |
| **C25** | the same gas taken between the same two states by two routes: what is the same, what differs? | none kept (D08 stated) | isothermal doubling of v for air at 300 K: w = −RT ln 2 ≈ −59.7 kJ/kg, Δe = 0, q = +59.7 kJ/kg; the isochoric→isobaric route gives a different q and w | `core.thermo.process_heat_work` vs trapezoid −∫p dv; p–v figure with shaded work; piston animation; explainer E3 | C23 equilibrium (collision time ~10⁻¹⁰ s), C26 path vs state functions, C27 (1.11), C28 equations of state (1.12) | N13 |
| **C35** | heat is path-dependent, but is there a heat-like quantity that is not? | D10 Gibbs relations | air 300 K → 600 K at constant p: Δs = C_p ln 2 ≈ 696 J/(kg K) by either Gibbs form | `core.thermo.perfect_gas_entropy_change`; sympy residual; T–s diagram with isotherms/isentropes; explainer E3 | C29 h (1.13), C30 C_p (1.14), C31 C_v (1.15), C33 entropy (1.16), N15 (1.17), C34 Clausius–Duhem (free expansion: q = 0, Δs = R ln 2) | C32, N14 |
| **C36** | why does sound cross air in a second per 340 m but water about 4× faster? | none (the book defers the proof to Ch. 15 §15.2) | p = Kρ^γ: c² = γp/ρ = 1.4 × 101325/1.225 → c ≈ 340 m/s; water K ≈ 2.2 GPa → c ≈ 1480 m/s | `core.thermo.sound_speed_from_eos` vs central difference; c vs bulk stiffness (log) with air and water marked | C37 α (1.20) (water: α changes sign at 4 °C) | — |
| **C40** | how do molecule counts become a law for a continuum? | D11 p = ρRT from pV = nk_BT | R_air = 8314/28.96 = 287 J/(kg K); ρ = 101325/(287 × 288.15) = 1.225 kg/m³ | `core.thermo.perfect_gas_density` vs the constants chain; isobars ρ(T); plotly 3-D p(v, T) surface; isothermal-atmosphere slider | R03 kelvin, C38 (1.21), C39 constants, C41 e = e(T) (van der Waals contrast; proof D13 stated), C48 α = 1/T (1.28), C62 p = p₀e^{−gz/RT}, C63 H = RT/g ≈ 7.3 km at 250 K | — |
| **C45** | a gas squeezed with no heat and no friction: how does p follow ρ? | D14 isentropic law | bicycle pump 1 → 2 bar: T = 288 × 2^{0.286} ≈ 351 K, ρ ratio 2^{0.714} ≈ 1.64 | `core.thermo.isentropic_pressure` vs step-by-step ds = 0 integration; p–ρ log–log (slope γ vs slope 1); explainer E3 | C42 R = C_p − C_v (1.23), C43 γ (1.24), N16 air γ = 1.4, C_p ≈ 1005, C44 adiabatic vs isentropic, C46 ratios (1.26), C47 c = √(γRT) = 340 m/s at 288 K | — |
| **C50** | push a parcel up a little and let go: what force brings it back? | D18 parcel equation | ocean thermocline dρ/dz = −0.01 kg m⁻⁴ (parcel incompressible), ρ = 1025 → N² ≈ 9.6×10⁻⁵ s⁻², period ≈ 11 min | `core.stratification.parcel_displacement` vs Euler–Cromer loop; three-parcel animation (stable/neutral/unstable); explainer E4 | (uses D37 from C20) | — |
| **C51** | what does the number N² tell you about a column? | none of its own (D18 in C50 gives it) | N = 0.01 s⁻¹ → period 2π/N ≈ 10.5 min; N² = −10⁻⁴ s⁻² → e-folding 1/√|N²| = 100 s | `core.stratification.brunt_vaisala_sq`, `stability_timescale`; ζ(t) for several N² with periods marked; live widget; explainer E4 | C52 regimes (cos / line / cosh), C60 seawater S (35 g/kg, linear EOS), N21 dρ_a/dz = −ρg/c², C61 ocean criterion (1.35) (ρg/c² ≈ 4.5×10⁻³ kg m⁻⁴) | — |
| **C54** | air cools with height; when is that stable? | D19 (★★★, sympy) | C_p = 1005 J/(kg K): Γ_a = −9.81/1005 = −9.76×10⁻³ K/m = **−9.8 K/km (Kundu)** = **+9.8 K/km (meteorology)**; standard troposphere dT/dz = −6.5 K/km > −9.8 K/km → stable; as meteorology lapse rates 6.5 < 9.8 → stable | `core.stratification.adiabatic_lapse_rate`; environment T(z) vs parcel adiabats figure with a convention table and both inequalities printed; lapse-rate slider; explainer E4 | C49 atmosphere from T(z) (USSA), C53 lapse rate in both conventions, N17 value in both conventions | — |
| **C55** | can we label air so its temperature label does not change when it rises? | D20 θ, D36 N² = (g/θ)dθ/dz | 500 hPa, 250 K → θ = 250 × 2^{0.286} ≈ 305 K | `core.stratification.potential_temperature` vs formula; inversion-strength slider (T, θ, N² coloured by regime); explainer E4 | C56 (1.32) (stated, D21), C57 stability from dθ/dz, N18 synthetic boundary-layer profile, C58 ρ_θ (1.33) | N19, N20, C59 |
| **C64** | can a law depend on whether we measure in metres or feet? | none | p = p₀ − ρgz: [ρgz] = kg m⁻³ · m s⁻² · m = Pa ✓; a wrong formula ρgz² fails | `core.units.dimensional_check`; exponent bars per term (match vs mismatch) | R01 SI units with `pint` | R02 |
| **C67** | how do we store "what dimensions does each variable have" so a computer can reason about it? | none kept (D26 stated) | pipe matrix (1.39): the 3×3 minor of columns (d, U, ρ) has determinant −1 ≠ 0 → r = 3 | `core.dimensional.dimensional_matrix`, `rank_by_minors` vs cofactor expansion; annotated heatmap of (1.39) with the minor highlighted; explainer E5 | C65 choosing variables (with and without μ: no collapse), C66 [q] and M, L, T, θ, C68 rank by minors | — |
| **C69** | why can 7 pipe variables be squeezed into 4 numbers? | D28 (★★★, sympy nullspace) | n − r = 7 − 3 = 4; Π₁ = Δp/ρU² for Δp = 100 Pa, ρ = 1000, U = 0.1 → 10, identical in cgs | `core.dimensional.pi_groups` vs `np.linalg.solve` exponent solve; units-invariance scatter; pipe-data collapse in Π coordinates; explainer E5 | N23 n − r, C70 exponent algebra (Π₁…Π₄ stated, D27), C71 combining groups (D33), C72 (1.40), N25 Δp ∝ Δx, C73 H ∝ R_uT/(gM_w) (D29), C75 E = KρD⁵/t², D ∝ t^{2/5} (D31), C76 S/I ∝ λ⁻⁴, blue/red ≈ 5.9 (D32) | N22, N24, C74 |

## 3. Section coverage
| § | Title | A (CORE) | B | C | CORE | RECAP | NOTE | SKIP |
|---|---|---|---|---|---|---|---|---|
| 1.1 | Fluid Mechanics | — | C01 | — | — | — | C01 | — |
| 1.2 | Units of Measurement | — | R01, R03 | R02 | — | R01, R02, R03 | — | — |
| 1.3 | Solids, Liquids, and Gases | — | C02, C03, C04 | N01, N02 | — | — | C02, C03, C04, N01, N02 | — |
| 1.4 | Continuum Hypothesis | C06 | C05, C07, N03 | — | C06 | — | C05, C07, N03 | — |
| 1.5 | Molecular Transport Phenomena | C12 | C08, C10, C11, N05, C13, C14 | C09, N04, N06 | C12 | — | C08, C09, C10, N04, C11, N05, C13, C14, N06 | — |
| 1.6 | Surface Tension | — | C15, N07, C16 | N08, N09 | — | — | C15, N07, C16, N08, N09 | — |
| 1.7 | Fluid Statics | C20 | C17, C18, C19, C21, N11, C22 | N10, N12 | C20 | — | C17, C18, N10, C19, C21, N11, C22, N12 | — |
| 1.8 | Classical Thermodynamics | C25, C35, C36 | C23, C24, C26, C27, C28, C29, C30, C31, C33, C34, N15, C37 | N13, C32, N14 | C25, C35, C36 | — | C23, C24, C26, N13, C27, C28, C29, C30, C31, C32, C33, C34, N14, N15, C37 | — |
| 1.9 | Perfect Gas | C40, C45 | C38, C39, C41, C42, C43, N16, C44, C46, C47, C48 | — | C40, C45 | — | C38, C39, C41, C42, C43, N16, C44, C46, C47, C48 | — |
| 1.10 | Stability of Stratified Fluid Media | C50, C51, C54, C55 | C49, C52, C53, N17, N18, C56, C57, C58, C60, N21, C61, C62, C63 | N19, N20, C59 | C50, C51, C54, C55 | — | C49, C52, C53, N17, N18, N19, C56, C57, C58, N20, C59, C60, N21, C61, C62, C63 | — |
| 1.11 | Dimensional Analysis | C64, C67, C69 | C65, C66, C68, N23, C70, C71, C72, N25, C73, C75, C76 | N22, N24, C74, S01, S02 | C64, C67, C69 | — | N22, C65, C66, C68, N23, C70, N24, C71, C72, N25, C73, C74, C75, C76 | S01, S02 (end of chapter) |

Totals: A 15 · B 70 · C 21 = 106. Every section has at least one row; sections §1.1, §1.2, §1.3 and §1.6 have no A item
and follow the placement rule in the header.

## 4. Prerequisites needing primers (concept or tool | needed by | why it is not CORE/RECAP)
Nothing is primed yet (`knowledge/primers.md` is empty). Rows marked **gloss** are needed only by B/C statements or by
demoted derivations: one plain sentence where the item is stated, no demo. All other rows are full 📎 primers (sentence
+ 2–4-line runnable demo) before the first A block that uses them.

| Concept or tool | Needed by | Why it is not CORE/RECAP |
|---|---|---|
| **Physics vocabulary** | | |
| stress = force per area; normal and tangential (shear) components; strain angle γ and strain rate | C12 (with C02, C03 stated) | pre-book mechanics, used without definition in §1.3 |
| elastic shear modulus G (solid: τ = Gγ) | C02 (B) | gloss: not in the book; our model for the solid in the C02 sketch |
| vapour pressure (why a liquid under tension boils/cavitates) | C03 (B) | gloss: mentioned in passing, taught in no chapter |
| Newton's second law F = ma and free-body diagrams | C20, C50, D05, D18, D34 | pre-book |
| weight = mg, gravitational acceleration g | C20 | pre-book |
| temperature as average molecular kinetic energy; Boltzmann constant k_B (informal, before C38) | C06 (D34) | used in §1.4 before §1.9 defines it |
| mole, kilomole, molecular weight, Avogadro's number (informal, before C39) | C06 (C04 spacing), C40, D11 | needed for molecular spacing before §1.9 |
| perfect-gas law as a working model (forward look at C40) | C25, C35 | the thermodynamics demos need a concrete e and p before §1.9 |
| specific heat (informal) and thermal diffusivity κ = k/ρC_p | C12 (C11 stated), E2 heat mode | gloss: C_p is defined only in §1.8 |
| boundary condition, no-slip at a wall, steady vs transient | C12, E2 | formal treatment in Ch. 4 §4.10 |
| internal energy and kinetic energy per unit mass | C25 | pre-book physics; the book uses e without defining it |
| light as waves: wavelength and colour; intensity ∝ amplitude² | C76 (B) | gloss: pre-book optics for Ex. 1.5 |
| similarity solution (why K comes from a full solution) | C75 (B) | gloss: full treatment is beyond Ch. 1 |
| **Mathematics** | | |
| ordinary derivative as a slope; d/dy of a profile | C12, C54 | pre-book calculus |
| partial derivative ∂/∂x | C20, D05 | pre-book calculus, first used in (1.7) |
| partial derivative with a variable held fixed, ( )_p notation | C35 (C30, C31 stated), C36 | thermodynamic notation, not school calculus |
| total vs partial derivative | D12 (demoted) | gloss |
| gradient vector ∇ and "down-gradient" | C12 (C10, C11 stated) | gloss: vector calculus comes formally in Ch. 2 §2.9 |
| dot product and unit normal vector | C03 (B), C10 (B) | gloss: Ch. 2 §2.8 formalises it |
| first-order Taylor expansion f(x + dx) ≈ f + f′dx | D05, D18 | pre-book |
| limits and orders of smallness (dz² ≪ dz) | C18 (B, D03 demoted), D05 | gloss for C18; used lightly in D05 |
| radius of curvature; principal radii | C16 (B) | gloss |
| definite integral; integrating a constant | D37 | pre-book |
| surface integral of pressure over a closed body (net force) | D37 | needed for Archimedes |
| natural log and exponential; e-folding | C40 (C62, C63 stated), C51, D14 | pre-book |
| separation of variables for a first-order ODE | D14 | pre-book ODEs |
| exponent rules (a^m a^n = a^{m+n}) | D20, C45 (C46 stated) | pre-book algebra |
| logarithmic differentiation d ln f = df/f | D36 | not school-standard |
| differentials; exact (d) vs inexact (δ) differentials | C25, C35 | the book's δ/d notation needs a primer |
| line integral along a path (∫p dv as area on a p–v diagram) | C25 | calculus-of-paths idea |
| product rule for differentials d(pv) = p dv + v dp | D10 | pre-book, easy to drop a term |
| chain rule (one and several variables) | D19 | pre-book |
| equality of mixed partials, exact differentials, Maxwell relations | D19 | the chapter never introduces them (analysis §9) |
| Gibbs free energy g = h − Ts (as a device) | D19 | not in Ch. 1, needed for the Maxwell relation |
| linear second-order ODE ζ″ + N²ζ = 0; cos, sinh/cosh solutions | C50, C51, D18 | pre-book ODEs |
| square root of a negative number → growth rate (imaginary N) | C51 | reading N² < 0 |
| inequalities under a sign change (multiplying by −1 flips >) | C54 (both lapse-rate conventions) | the convention conversion is exactly this move |
| Poisson counting statistics; relative noise ∝ N^{-1/2} | C06, D35 | statistics tool |
| mean square speed; Gaussian (Maxwell) velocity components | D34 | kinetic theory tool |
| power laws and log–log plots | C06, C45 | reading slopes as exponents |
| finite differences: central difference, FTCS scheme, stability limit DΔt/Δy² ≤ ½ | C12, C36 | numerical tool (Ch. 10 §10.2 formal) |
| trapezoid rule | C25 | numerical integration |
| explicit ODE time stepping (Euler, Euler–Cromer) | C20, C45, C50 | from-scratch versions |
| matrices, determinants, minors, cofactor expansion | C67 | linear algebra |
| linear independence and rank | C67, C69 | linear algebra |
| solving linear systems Ax = b | C69 (C70 stated) | linear algebra |
| null space and rank–nullity theorem | C69, D28 | linear algebra |
| similar triangles | C74 (C) | gloss |
| **Python** | | |
| numpy arrays, vectorised arithmetic, broadcasting | C06 onward | first Python cell |
| `np.linspace`, `np.logspace` | C06 | log-spaced box sizes |
| `np.random.default_rng(seed)`: `normal`, `poisson` | C06 | seeded sampling |
| `np.gradient` (and `edge_order=2`) | C12, C54 | finite differences on arrays |
| vectorised `np.where` / `np.select` | C51, C55 | regime labels |
| `matplotlib`: `subplots`, `plot`, `loglog`, `semilogy`, `fill_between`, `arrow`, `imshow` | C06 onward | first figure cells |
| `fluidpy.core.anim.animate` + `show_animation` | C06, C12, C25, C50 | first animation |
| `fluidpy.core.interact.slider_figure` (plotly) | C06 | first slider figure |
| `plotly.graph_objects.Surface` | C40 | 3-D p(v, T) surface |
| `pint`: `Q_`, `.to()`, `.dimensionality`, offset units | R01, R03, C64, C67 | units in code |
| `scipy.integrate.solve_ivp` | C50 | nonlinear parcel ODE |
| `scipy.integrate.cumulative_trapezoid` / `np.trapezoid` | C25 | path integrals |
| `sympy`: `symbols`, `diff`, `simplify`, `Matrix`, `nullspace`, `Rational` | D19, D28, C35, C67 | symbolic checks |
| `itertools.combinations` | C67 | enumerating minors |
| `np.linalg.det`, `matrix_rank`, `solve` | C67, C69 | numerical linear algebra |
| functions passed as arguments / `lambda` | C20, C36 | `integrate_hydrostatic(rho_fn)`, `sound_speed_from_eos(f)` |
| dictionaries as variable lists | C67, C69 | `pi_groups` input |
| f-strings with units; `assert np.allclose` | C06 onward | from-scratch checks |
| `show_viz(chapter, slug)` | first explainer (C06) | embedding explainers |

## 4b. Derivations (parsed by tools: ID first, CORE id in a column, ★★★ for hard, explainer slugs backticked in the LAST column)
Only derivations whose parent is an A item survive, plus the four added ones (D34–D37, re-parented). 12 rows.

| ID | Result (Eq.) | CORE | Difficulty | Steps | Tools used | Traps | Shown in |
|---|---|---|---|---|---|---|---|
| D05 | hydrostatic law (1.8) | C20 | ★ | 5 | first-order Taylor (primer), partial derivative (primer), Pascal (1.7) stated in the same block | sign of g with z upward; writing dp/dz before (1.7) shows p depends on z only | notebook |
| D10 | Gibbs relations (1.18) | C35 | ★ | 5 | product rule for differentials (primer), T ds = dq (N15, stated), de = dq − p dv (C27, stated in C25), h = e + pv (C29, stated) | concluding they hold only for reversible processes (they link state functions, so they hold for any process) | notebook · `heat_work_paths` |
| D11 | p = ρRT (1.22) from pV = nk_BT (1.21) | C40 | ★ | 6 | (1.21) and the constants (C38, C39, stated in the same block), unit bookkeeping with kmol (primer) | mixing mol and kmol (factor 1000); R_u vs R | notebook |
| D14 | isentropic law p/ρ^γ = const (1.25) (Exercise 1.11) | C45 | ★★ | 9 | C35, (1.22) from C40, (1.23)–(1.24) stated in the same block, separation of variables (primer), logarithms (primer) | dv/v = −dρ/ρ sign; assuming constant C_p, C_v silently; dividing the two Gibbs forms in the wrong order | notebook · `heat_work_paths` |
| D18 | parcel equation ζ″ + N²ζ = 0 and (1.29) (Exercise 1.13) | C50 | ★★ | 11 | Newton's second law (primer), buoyancy D37 (C20), first-order Taylor (primer), linear second-order ODE (primer) | sign of buoyancy; using the environment's gradient for the parcel; keeping O(ζ²) terms; ρ_p ≈ ρ(z_o) only in the inertia term | notebook · `parcel_stability` |
| D19 | adiabatic lapse rate (1.30) (Exercise 1.14, no perfect-gas relations) | C54 | ★★★ | 13 | Gibbs free energy (primer), Maxwell relations (primer), chain rule (primer), C35, C20, α (1.20, stated in C36), C_p (1.14, stated in C35) | Maxwell-relation sign; using dh = C_p dT (true only for a perfect gas); the sign of the result depends on the convention: with Kundu's Γ ≡ dT/dz, Γ_a = −gαT/C_p < 0 and stability is dT/dz > Γ_a; with the meteorology Γ ≡ −dT/dz the same result reads Γ_a = +gαT/C_p > 0 and stability is Γ < Γ_a, so negating Γ must also flip the inequality | notebook · `parcel_stability` |
| D20 | potential temperature (1.31) | C55 | ★ | 4 | exponent rules (primer), (1.26) stated in C45 | p_o (reference pressure) confused with p₀ at z = 0 | notebook · `parcel_stability` |
| D28 | Buckingham Π theorem (1.37) (book states it only) | C69 | ★★★ | 12 | null space and rank–nullity (primer), C64, C67 with rank (C68 stated there), scaling of base units (primer) | thinking the groups are unique; forgetting the relation must be dimensionally homogeneous and the variable list complete | notebook · `buckingham_pi_machine` |
| D34 | kinetic pressure p = ⅓ (N/V) m⟨u²⟩ (added: the book states the idea only) | C06 (re-parented; was under the pressure-as-impacts row) | ★★ | 8 | Newton's second law (primer), Gaussian velocity components and mean square speed (primer) | factor 2 from the bounce; ⟨u_x²⟩ = ⅓⟨u²⟩; counting only molecules moving toward the wall | notebook |
| D35 | relative density noise ≈ N^{-1/2} ∝ (δV)^{-1/2} (added) | C06 | ★ | 5 | Poisson counting statistics (primer), power laws (primer) | noise ∝ L^{-3/2} in box side, not L^{-1/2}; mixing absolute and relative noise | notebook · `continuum_averaging_volume` |
| D36 | N² = (g/θ)dθ/dz for a perfect-gas atmosphere (added: links (1.29) and (1.32)) | C55 (re-parented; was under the dθ/dz stability row) | ★★ | 8 | logarithmic differentiation (primer), C51, (1.26) and (1.23)–(1.24) stated in C45, (1.32) stated in the same block | using the environment's density gradient for the parcel; (γ − 1)/γ = R/C_p identity; sign of dp/dz | notebook · `parcel_stability` |
| D37 | buoyancy force = ρ_fluid gV from p = p₀ − ρgz (added: grounds D18) | C20 (re-parented; was under the uniform-density row) | ★ | 6 | surface integral of pressure over a closed body (primer), definite integral (primer), (1.9) stated in the same block | side-face forces do cancel; the net force is up; the body's own density is irrelevant to the buoyancy | notebook |

Counts: ★ 6 (D05, D10, D11, D20, D35, D37) · ★★ 4 (D14, D18, D34, D36) · ★★★ 2 (D19 in `parcel_stability`, D28 in
`buckingham_pi_machine`). Both surviving ★★★ derivations belong to explainers and are shown in them, with sympy checks.

## 4c. Derivations demoted to statements
These results are **given, not derived**, inside the named A block (a paragraph, the equation, a number). The first
column is the A parent (bold, so the parser does not read the row as an item). Former D ids are kept for traceability
only. 25 rows.

| A parent | Former D | Result stated (Eq.) | Was | Stated inside the A block as |
|---|---|---|---|---|
| **C20** | D01 | spherical jump p_i − p_o = 2σ/R | ★, 6 steps | N07: "tension around the rim 2πRσ balances the pressure on the projected disc πR²"; R = 1 mm droplet → Δp ≈ 146 Pa |
| **C20** | D02 | Laplace jump (1.5) | ★★, 9 steps | C16: equation with signed principal radii; sphere and cylinder special cases |
| **C20** | D03 | pressure is isotropic (1.6) | ★★, 10 steps | C18: the wedge's weight shrinks like dz², the face forces like dz, so in the limit the face pressures must match |
| **C20** | D04 | Pascal's law (1.7) | ★, 5 steps | C19: the two horizontal face forces on a cube can only cancel if ∂p/∂x = ∂p/∂y = 0 (the same move as D05, horizontally) |
| **C20** | D06 | uniform-density hydrostatics (1.9) and p − p₀ = ρgh | ★, 5 steps | C21: integrate (1.8) with ρ constant; 10 m of water ≈ 1 atm |
| **C20** | D07 | capillary rise h = 2σ sin α/(ρgR) (Ex. 1.1) | ★, 7 steps | C22: surface-tension pull 2πRσ sin α holds up the weight ρgπR²h; 1 mm radius water tube → h ≈ 15 mm; the book's α is the complement of the contact angle |
| **C25** | D08 | reversible first law (1.11) | ★, 6 steps | C27: slow boundary work on the particle is −p dv per unit mass, so de = dq − p dv |
| **C35** | D09 | heat per degree equals C_v (constant v) and C_p (constant p) | ★★, 8 steps | C32 (one sentence): only along reversible isochoric / isobaric paths; stirring breaks it |
| **C45** | D12 | R = C_p − C_v (1.23) | ★★, 7 steps | C42: h = e + RT differentiated with e, h functions of T only; air 1005 − 718 = 287 J/(kg K) |
| **C40** | D13 | perfect gas ⇔ e = e(T), h = h(T) (Exercise 1.10) | ★★★, 14 steps | C41: stated with the reason in words (Gibbs + p = ρRT force (∂e/∂v)_T = 0); a van der Waals gas as the counter-example figure |
| **C45** | D15 | isentropic ratios (1.26) | ★, 5 steps | C46: substitute p = ρRT into (1.25); bicycle-pump number |
| **C45** | D16 | c = √(γRT) (1.27) | ★, 5 steps | C47: (∂p/∂ρ)_s of p = Kρ^γ is γp/ρ = γRT; 340 m/s at 288 K; Newton's isothermal √(RT) = 287 m/s is the trap |
| **C40** | D17 | α = 1/T (1.28) | ★, 4 steps | C48: ρ = p/RT at fixed p gives −(1/ρ)∂ρ/∂T = 1/T; air at 300 K → 3.3×10⁻³ K⁻¹ |
| **C55** | D21 | dθ/dz and lapse rates (1.32) | ★★, 9 steps | C56: log-differentiate (1.31), use (1.8) and p = ρRT; stated in both conventions: (T/θ)dθ/dz = dT/dz − Γ_a (Kundu) = Γ_a,met − Γ_met (meteorology) |
| **C55** | D22 | potential density (1.33) | ★, 4 steps | C58: the same reference-pressure trick with the density ratio of (1.26) |
| **C55** | D23 | θρ_θ = p_o/R and (1.34) | ★★, 7 steps | N20/C59 (one sentence each): exponents add to 1, so the product is constant and the log-gradients are equal and opposite |
| **C51** | D24 | ocean criterion (1.35) with dρ_a/dz = −ρg/c² | ★★, 8 steps | N21, C61: the parcel only compresses, (∂ρ/∂p)_s = 1/c² and dp/dz = −ρg; ρg/c² ≈ 4.5×10⁻³ kg m⁻⁴; "same sign as", not equality |
| **C40** | D25 | isothermal p = p₀e^{−gz/RT} and H = RT/g | ★, 6 steps | C62, C63: (1.8) with ρ = p/RT separates; H ≈ 7.3 km at 250 K |
| **C67** | D26 | rank of the pipe matrix (1.39) is 3 | ★, 6 steps | C68: one nonzero 3×3 minor, columns (d, U, ρ), with determinant −1 is enough, shown as the block's worked number |
| **C69** | D27 | Π₁ = Δp/ρU², Π₂ = Δx/d, Π₃ = ε/d, Π₄ = μ/ρUd | ★★, 10 steps | C70: the four groups given with a dimension check; the exponent solve for Π₁ is run as code (`np.linalg.solve`), not derived |
| **C69** | D29 | Ex. 1.2: H = const·R_uT_o/(gM_w) | ★★, 9 steps | C73: one group → it must be a constant; the constant is 1 by comparison with C63 |
| **C69** | D30 | Ex. 1.3: A² + B² = C² | ★★, 8 steps | C74 (one sentence): area = C²φ(β) and the two similar sub-triangles add up |
| **C69** | D31 | Ex. 1.4: E = KρD⁵/t² and D ∝ t^{2/5} | ★, 6 steps | C75: one group from (E, ρ, D, t); K needs the full similarity solution |
| **C69** | D32 | Ex. 1.5: S/I = V²φ₃(n_s)/(d²λ⁴) | ★★, 8 steps | C76: r = 2, physics fixes V² and d⁻², so λ⁻⁴ follows; blue (450 nm) / red (700 nm) ≈ 5.9 |
| **C69** | D33 | combined groups Δp d²ρ/μ² = Π₁/Π₄², ε/Δx = Π₃/Π₂ | ★, 4 steps | C71: products and powers of groups are groups; a set stays at n − r independent members |

## 5. Interactive explainers (4–5 + backup)

### E1 · continuum_averaging_volume
- **CORE:** C06 · also shows C05 (molecules view), C07, N03 (Kn strip), C24 (fluid-particle window), C04 (medium chips)
- **Confusion removed:** "density at a point" sounds meaningless when matter is made of molecules. The reader sees that a well-defined value exists only inside a *window* of box sizes: larger than the molecular spacing, smaller than the scale over which the flow itself varies.
- **Why interactive:** the idea is a sweep over eight decades of box size with random noise. Dragging the box and resampling shows noise that no static curve conveys: the same box can read 30 % high one moment and 20 % low the next, then settles as it grows.
- **Stage:** (1) a zoomable cross-section of gas or liquid with molecules as dots and the sampling cube drawn on it, with a background density gradient visible at large zoom; (2) measured ρ vs box side L on log axes: a faint expected band (±N^{-1/2}), bold samples so far, the plateau and the macroscopic drift; (3) Kn = l/L strip with continuum/slip/free-molecular bands.
- **Controls:** box side L (log slider) · medium (chips: air at sea level, air at 80 km, water) · macroscopic gradient strength · "resample" button · body size for Kn (optional).
- **Equations shown:** ρ = δm/δV (continuum definition, §1.4), relative noise ≈ N^{-1/2} (D35), Kn = l/L (§1.4); kinetic pressure (D34) in the Equations tab only.
- **Mirrors:** `fluidpy.ch01_introduction.sample_density`, `density_noise_expected`, `density_expected`, `knudsen_number`, `mean_free_path_jennings`.
- **Derivations:** D35.
- **Depth features:** Explain tab (N = nL³ computed with your numbers → noise → verdict), synced Code tab, + linked views (3), presets (air / high altitude / water / microchannel), status verdict ("🎲 molecular noise ±x %" / "✅ continuum plateau" / "📈 box sees the flow's own variation"), inspector (click a sample: count, mass, δm/δV arithmetic), "Right now" notes with a highlighted table of real sizes (virus, aerosol, raindrop, pipe, cloud).
- **Follows reference:** `angular_frequency_explorer_1.html` (linked views, modes, highlighted real-world table, Right-now notes).
- **Aha:** a continuum value is not a point property of matter. It is the plateau of an average, and it exists because there is a wide range of box sizes between the molecules and the flow.

### E2 · viscosity_momentum_diffusion
- **CORE:** C12 · also shows C08 (diffusion picture), C14 (ν vs μ), C10 and C11 (species and heat modes), C02 (tracer lines keep shearing), C13 (fluid presets at their temperatures), N05
- **Confusion removed:** viscosity pictured as "stickiness" or friction between layers, and μ confused with ν. Momentum *diffuses* from the moving plate, and ν = μ/ρ, not μ, sets how fast. Air (small μ) spreads motion faster than water.
- **Why interactive:** the phenomenon is a transient: a profile creeping into the gap and settling to the linear Couette state with uniform τ. Playing, scrubbing and switching fluids shows the time scale h²/ν at work, and a static figure cannot show the process.
- **Stage:** (1) the gap between plates with dyed vertical tracer lines shearing and a few molecules hopping between layers carrying momentum colour; (2) u(y, t) with a faint steady linear profile (ghost) and the bold current profile; (3) wall shear stress τ_w(t) = μ ∂u/∂y at both plates, both tending to μU/h, with a marker at t = h²/ν.
- **Controls:** fluid preset or μ (log) · ρ · gap h · plate speed U · mode (momentum / heat / species).
- **Equations shown:** (1.3), (1.4), (1.1), (1.2); model diffusion equation ∂f/∂t = D∂²f/∂y² (labelled "ours; derived in Ch. 4").
- **Mirrors:** `fluidpy.core.diffusion.ftcs_diffusion_1d`, `couette_startup_profile`, `fluidpy.ch01_introduction.newton_shear_stress`, `kinematic_viscosity`, `thermal_diffusivity`, `wall_shear_history`.
- **Derivations:** none (the book derives nothing in §1.5; the Explain tab works τ, ν and h²/ν out with numbers).
- **Depth features:** Explain tab, synced Code tab, + transport (play/step/scrub, end-of-run summary card), linked views (3), modes (momentum / heat / species with D = ν, k/ρC_p, κ_m), presets (air, water, glycerine, honey), status verdict ("⏳ momentum still spreading: t/(h²/ν) = 0.12" / "✅ steady Couette: τ uniform").
- **Follows reference:** `forced_damped_vibrations.html` (system + graph on one clock, transient → steady state, "at the current time" values, regime-dependent interpretation).
- **Aha:** honey has about 10⁶ times the μ of air, but ν decides how fast motion spreads, and air's ν is about 15 times water's.

### E3 · heat_work_paths
- **CORE:** C25, C35, C45 · also shows C26 (path vs state), C27 (1.11), C33 (entropy), C34 (irreversible mode), C44 (adiabatic vs isentropic), C46
- **Confusion removed:** heat and work treated as things a gas "has"; δ vs d; why entropy is a property although heat is not; adiabatic vs isentropic.
- **Why interactive:** path dependence only clicks when the reader drags a path between the *same* two states and watches q and w change while Δe and Δs stay fixed. Term bars updating live beat any pair of static diagrams.
- **Stage:** (1) piston–cylinder with gas: the piston moves along the chosen path, a heater glows when δq > 0, an insulation jacket appears for adiabatic legs, a stirrer runs in irreversible mode; (2) p–v diagram with isotherms (faint), the chosen path bold "so far", and the work area shaded; (3) T–s diagram with the same path, with its heat area shaded.
- **Controls:** path (chips: isothermal · isochoric→isobaric · isobaric→isochoric · isentropic + isochoric · custom corner) · end state v₂/v₁ · T₂/T₁ · gas (γ = 5/3 or 7/5) · mode (reversible / irreversible: stirring or free expansion).
- **Equations shown:** (1.10), (1.11), (1.13), (1.16), (1.17), (1.18), (1.25), and the Clausius–Duhem inequality (with actual δq).
- **Mirrors:** `fluidpy.core.thermo.process_path`, `path_heat_work_totals`, `process_heat_work`, `entropy_change_reversible`, `perfect_gas_entropy_change`, `isentropic_pressure`, `irreversible_process`.
- **Derivations:** D10, D14. (The reversible first law (1.11) is stated in the Equations tab and worked with numbers in the Explain tab; its former derivation D08 is demoted, §4c.)
- **Depth features:** Explain tab (w = −∫p dv leg by leg with your numbers, q = Δe − w, Δs two ways), synced Code tab, + linked views (3), transport along the path, term bars (q, w, Δe, and Δs vs ∫δq/T, which sum or compare live), presets (the five classic paths), modes (reversible / irreversible), status ("🔁 same Δe, different q" / "⚠️ irreversible: Δs > ∫δq/T").
- **Follows reference:** `amplitude_phase_second_order_II_3.html` (several windows linked by one state, numbered live derivation in the explanation), with term bars after `fid_formula_lab.html`.
- **Aha:** along every path between the same two states, Δe and Δs are identical while q and w are not. That is what "state function" means.

### E4 · parcel_stability
- **CORE:** C50, C51, C54, C55 · also shows C49 (environment builder), C52 (status), C53 and N17 (both lapse-rate conventions), C56, C57 (θ view), C60, N21, C61 (ocean mode)
- **Confusion removed:** "colder air aloft means unstable" (wrong: the standard troposphere cools with height and is stable). Stability compares the environment's lapse rate with the parcel's own adiabatic cooling, which θ does in one step. N² < 0 means exponential runaway, not a faster oscillation. And a sign trap: the book and the meteorology literature write the same criterion with opposite signs and opposite inequality directions.
- **Why interactive:** the reader has to *push a parcel* and watch it oscillate, stay put or escape while the environment profile is reshaped. The link between a profile's slope, the parcel adiabat and ζ(t) is a three-way dependence that sliders make obvious; flipping the convention toggle while nothing physical changes shows that only the bookkeeping differs.
- **Stage:** (1) a vertical column coloured by θ (or ρ) with a parcel that moves on the clock, drawn with its own adiabat through the release height; (2) the profile panel: environment T(z) bold, parcel dry adiabat dashed (or ρ and ρ_a in ocean mode), θ(z) toggle; (3) ζ(t) with the linear solution (cos / linear / cosh) as ghost and the nonlinear `parcel_ode_atmosphere` path bold, period 2π/N or e-folding time marked.
- **Lapse-rate conventions (user decision, binding):** the physics is computed with Kundu's Γ ≡ dT/dz (slider in dT/dz, K/km, about −15 to +10; Γ_a ≈ −9.8 K/km). The criterion is displayed **in the form computed, with the inequality direction visible**: "stable ⇔ dT/dz > Γ_a: −6.5 > −9.8 K/km ✓". A **convention panel** (toggle Kundu ↔ meteorology, both rows always visible in the Explain tab and the Equations tab) shows the standard meteorology form alongside: Γ ≡ −dT/dz = +6.5 K/km, Γ_a = +9.8 K/km, "stable ⇔ Γ < Γ_a: 6.5 < 9.8 ✓", with the conversion Γ_met = −Γ_Kundu and the flipped inequality spelled out. Every lapse-rate readout shows both numbers; the toggle changes which form is primary, never the physics.
- **Controls:** environment lapse rate dT/dz in K/km (or dρ/dz in ocean/lab mode) · inversion toggle/strength · release displacement ζ₀ · mode (atmosphere / ocean / lab tank) · convention (Kundu dT/dz / meteorology −dT/dz).
- **Equations shown:** parcel equation (§1.10), (1.29), (1.30) in both conventions, (1.31), (1.32), (1.35) in ocean mode, N² = (g/θ)dθ/dz.
- **Mirrors:** `fluidpy.core.stratification.brunt_vaisala_sq`, `brunt_vaisala_sq_from_lapse`, `brunt_vaisala_sq_from_theta`, `parcel_displacement`, `parcel_ode_atmosphere`, `adiabatic_lapse_rate`, `potential_temperature`, `classify_stability`, `stability_timescale`, and the new convention helper (§8).
- **Derivations:** D18, D19 (★★★), D20, D36.
- **Depth features:** Explain tab (N² with your lapse rate → N → period or e-folding → θ gradient → the criterion in both conventions), synced Code tab, + linked views (3), transport (release/play/scrub, end-of-run card), presets (dry adiabatic = neutral, standard atmosphere −6.5 K/km, nocturnal inversion, superadiabatic surface layer, ocean thermocline), status verdict ("🌊 stable: period 11 min" / "⚖️ neutral" / "🚀 unstable: e-folds in 3 min"), modes (atmosphere / ocean / lab tank), "Right now" notes.
- **Follows reference:** `forced_damped_vibrations.html` (system + x(t) graph + regime-dependent interpretation, "reading the current setting").
- **Aha:** air that gets colder with height can still be stable. What matters is whether it cools more slowly than a rising parcel cools by expanding (dT/dz > −9.8 K/km, or equivalently a meteorology lapse rate below +9.8 K/km), and that is exactly whether θ increases upward.

### E5 · buckingham_pi_machine
- **CORE:** C64, C67, C69 · also shows C65, C66, C68, N23, C70, C71, C72 (pipe preset), C73, C75, C76 (presets)
- **Confusion removed:** the Π theorem as a magic recipe. The reader sees it as linear algebra: groups are null-space vectors of the dimensional matrix, their number is n − r, a different repeating set gives a different but equivalent basis, and forgetting a variable changes everything.
- **Why interactive:** toggling variables on and off and swapping the repeating set changes the matrix, its rank and the groups at once. The reader experiments with the method instead of watching it done once.
- **Stage:** (1) variable chips with their dimension vectors (on = column in the matrix); (2) the dimensional matrix with the chosen r × r repeating minor highlighted and its determinant shown (red when singular); (3) the resulting groups as monomials with exponents, plus a units-change check: switch SI → cgs → imperial and every Π value stays the same (C64).
- **Controls:** problem preset (pipe Δp · pendulum period · sphere drag · Ex. 1.2 scale height · Ex. 1.4 blast · Ex. 1.5 Rayleigh) · variable toggles · repeating-set selector · unit system.
- **Equations shown:** (1.36), (1.37), (1.38), (1.39), (1.40), n − r, the exponent system of Π₁.
- **Mirrors:** `fluidpy.core.dimensional.dimensional_matrix`, `rank_by_minors`, `minor_determinant`, `solve_exponents`, `pi_groups`, `groups_independent`, `group_value`, `rescale_units`.
- **Derivations:** D28 (★★★). (Minors, rank and the exponent solve are worked as live arithmetic in the Explain tab; their former derivations D26, D27 are demoted, §4c.)
- **Depth features:** Explain tab (matrix → minors → r → n − r → each exponent solve with numbers → final dimension check), synced Code tab, + presets (the classic problems), inspector (click a matrix cell or a group: its exponent arithmetic), live status ("✅ 4 independent groups" / "⚠️ singular repeating set" / "⚠️ a variable is dimensionless on its own"), linked views (3).
- **Follows reference:** `amplitude_phase_second_order_II_3.html` (numbered live derivation), with presets and "formula with numbers plugged in" after `stride_padding_playground.html`.
- **Aha:** dimensionless groups are the null space of a small matrix. Buckingham's n − r is just rank–nullity, and every group's value survives a change of units.

### B1 · pressure_in_still_fluid (backup)
- **CORE:** C20 · also shows C17, C19, C21, C22, C16 (capillary mode)
- **Confusion removed:** gauge vs absolute pressure; "pressure pushes down only"; why a submerged block feels a net upward force; why water climbs a thin tube.
- **Why interactive:** drag a fluid cube or block through a layered tank and watch the face-pressure arrows: horizontal ones cancel, vertical ones differ by the weight, and in a different fluid the difference is buoyancy. Shrinking the tube radius makes the meniscus lift water higher.
- **Stage:** tank with 1–2 layers and a draggable cube/block with pressure arrows · p(z) (absolute/gauge) · mode "capillary tube" with meniscus, p along E–F and h vs R.
- **Controls:** upper-layer density · layer thickness · probe depth · block density · mode (tank / capillary) with R and α.
- **Equations shown:** (1.5), (1.7), (1.8), (1.9), Ex. 1.1, buoyancy (D37).
- **Mirrors:** `fluidpy.core.statics.hydrostatic_pressure_uniform`, `integrate_hydrostatic`, `layered_pressure`, `gauge_pressure`, `buoyancy_force`, `net_pressure_force_on_box`; `fluidpy.ch01_introduction.capillary_rise`, `laplace_pressure_jump`.
- **Derivations:** D05, D37 (added to those rows' "Shown in" only if this backup is built).
- **Depth features:** Explain, Code, + linked views, presets (water, oil over water, mercury manometer, glass capillary), inspector (click a face: p × area), term bars (top force, bottom force, weight, buoyancy), modes.
- **Follows reference:** `angular_frequency_explorer_1.html` (modes, presets, Right-now notes).
- **Aha:** buoyancy is nothing but the pressure difference between a body's bottom and top faces. The fluid's hydrostatic law does all the work.

## 6. Python animations and interactive figures (CORE ID → what, why, player/figure kind)
Only A items get their own; B items reuse their parent's figure (an overlay or a panel) when that is cheap, otherwise
the equation and a number are enough. Every A item has at least one visual listed here or in §2b.

| Kind | CORE | What moves / what the slider controls | Why | Player / figure |
|---|---|---|---|---|
| Animation | C06 | zoom-in: the sampling cube grows from 1 nm to 1 m over a gas with a density gradient; ρ estimate trace | the noise → plateau → drift sequence is a process of scale | `player="frames"` (step through decades) |
| Animation | C12 (C02, C08, N05 as panels) | u(y, t) diffusing from a suddenly moved plate with a τ(y) panel (FTCS, D = ν); a small inset of a solid vs fluid element under the same stress | diffusion is a transient; the settling to a linear profile must be seen | `player="video"` |
| Animation | C25 (C26 inside) | piston following two paths between the same states with running q, w, Δe counters | path dependence seen as motion | `player="frames"` (leg by leg) |
| Animation | C50 | three parcels released at once in stable, neutral and unstable columns; ζ(t) traces | cos vs straight line vs cosh is a behaviour in time | `player="video"` |
| Plotly slider | C06 (C07, N03) | altitude (0–100 km): mean free path l and Kn vs body size with regime bands | the continuum's failure depends on two scales at once | `slider_figure`, 30 steps |
| Plotly slider | C20 (C21) | upper-layer density of a two-layer tank: p(z) kink and block face forces | how density sets the slope of p(z), and buoyancy as a face-force difference | `slider_figure`, 25 steps |
| Plotly slider | C40 (C62, C63) | temperature T: isothermal p(z) and H vs the USSA table | how good the isothermal model is | `slider_figure`, 20 steps |
| Plotly slider | C54 (C53, N17) | environment dT/dz from −15 to +10 K/km: T(z) vs the dry adiabat; the legend prints the lapse rate in both conventions and the verdict with its inequality ("dT/dz = −6.5 > −9.8 K/km ⇔ Γ_met = 6.5 < 9.8 K/km: stable") | the convention is learnt by watching both numbers and both inequalities move together | `slider_figure`, 25 steps |
| Plotly slider | C55 (C57, N18) | inversion strength: T(z), θ(z) and N²(z) coloured by regime | reading stability from θ rather than T | `slider_figure`, 20 steps |
| Plotly 3-D | C40 (C28) | p(v, T) surface with a state dot | "two properties fix the state" is a 3-D idea | `go.Surface`, height 520 |
| Static figures | C35, C36, C45, C51, C64, C67, C69 | T–s diagram (C35); c vs bulk stiffness (C36); p–ρ log–log slopes γ vs 1 (C45); ζ(t) for several N² with periods marked (C51); exponent bars per term (C64); annotated heatmap of (1.39) with the nonzero minor (C67); units-invariance scatter + pipe-data collapse (C69) | one relationship each | matplotlib |
| Live widget (kernel only, plus the C51 static figure for the page) | C51 | N² and ζ₀ with the nonlinear `parcel_ode_from_gradients` | free exploration beyond the precomputed figure | `live(...)` |

Optional overlays for B items (cheap, only if the builder has time): Bingham/Maxwell curves on the C12 inset (N01);
Δp vs R curve beside C16; blast-front D ∝ t^{2/5} log–log line beside C75; Rayleigh spectrum strip beside C76.

## 7. From-scratch moments
A items only. Every section with a computable A item has at least one; each is followed by
`assert np.allclose(mine, library)`.

| § | CORE | Hand-written version | Tested function it must agree with |
|---|---|---|---|
| 1.4 | C06 | count seeded random positions inside a cube, ×m/L³, repeat for the std | `ch01.sample_density` (same seed) |
| 1.5 | C12 | central differences (u[i+1] − u[i−1])/(2Δy) with one-sided second-order ends, ×μ | `ch01.shear_stress_profile` |
| 1.7 | C20 | explicit Euler march p[k+1] = p[k] − ρ(z_k)gΔz on a fine grid | `core.statics.integrate_hydrostatic` (to O(Δz)) |
| 1.8 | C25 | trapezoid sum of −p dv along each leg, Δe = c_vΔT, q = Δe − w | `core.thermo.process_heat_work` |
| 1.8 | C36 | (p(ρ + Δρ) − p(ρ − Δρ))/(2Δρ) on p = Kρ^γ, then √ | `core.thermo.sound_speed_from_eos` and `perfect_gas_sound_speed` |
| 1.9 | C40 | ρ = p/(k_B A_o/M_w · T) from the constants chain | `core.thermo.perfect_gas_density` |
| 1.9 | C45 | integrate dp/dρ = γp/ρ with small steps from (ρ₀, p₀) | `core.thermo.isentropic_pressure` |
| 1.10 | C50 | Euler–Cromer loop for ζ″ = −N²ζ | `core.stratification.parcel_displacement` |
| 1.10 | C54 | `np.gradient` of a dry parcel's T along its lifted path (T(p) from (1.26) + hydrostatic p(z)), compared with −g/C_p; then negate for the meteorology value | `core.stratification.adiabatic_lapse_rate` |
| 1.10 | C55 | θ = T(p_ref/p)^{(γ−1)/γ} | `core.stratification.potential_temperature` |
| 1.11 | C67 | cofactor expansion of every 3×3 minor of (1.39) with `itertools.combinations` | `core.dimensional.rank_by_minors` and `np.linalg.matrix_rank` |
| 1.11 | C69 | `np.linalg.solve` on the 3×3 exponent system for Π₁ | `core.dimensional.solve_exponents` / `pi_groups` |

(§1.1, §1.2, §1.3 and §1.6 have no A item, so no from-scratch moment; §1.2's R01 cell still converts units with `pint`.)

## 8. Notes for the implementer
`fluidpy/` is committed WIP (`ch01_introduction.py` + `core/{thermo,dimensional,stratification,statics,diffusion}.py`)
and already contains nearly every function the first curation planned. Only what A items and the explainers need is
listed as required; everything else is optional.

**Required (A items and explainers)**
1. **Parity namespace.** `tools/shot.py` exposes `ch01.<fn>`; the chapter module already re-exports from `core`.
   Confirm every function an explainer mirrors (§5 "Mirrors" lists) is reachable as `ch01.<name>`, including
   `parcel_ode_atmosphere`, `stability_timescale`, `path_heat_work_totals`, `irreversible_process`, `minor_determinant`,
   `group_value`, `rescale_units`, `couette_startup_profile`, `layered_pressure`.
2. **Lapse-rate conventions (user decision).** Keep computing with Kundu's Γ ≡ dT/dz: `adiabatic_lapse_rate()` returns
   Γ_a in K/m, negative for air (≈ −9.76e-3), and environment arguments are named `dT_dz`. Add a small, tested helper
   for the standard meteorology convention, which the notebook and E4 always show alongside:
   `lapse_rate_convention(dT_dz, convention="kundu" | "meteorology")` → Γ in that convention (Γ_met = −dT/dz), and
   `lapse_rate_stability(dT_dz, Gamma_a=None, convention="kundu")` → `(verdict, inequality_text)` where the text is the
   criterion in that convention with numbers ("−6.5 > −9.8 K/km" or "6.5 < 9.8 K/km"). E4 parity rows cover both
   conventions. Docstrings of `lapse_rate`, `adiabatic_lapse_rate`, `potential_temperature_gradient` state both
   conventions and the inequality flip in their first lines (the WIP docstrings already do most of this). Tests check
   Γ_a = −g/C_p (Kundu), +g/C_p (meteorology), and that both conventions give the same verdict on a profile sweep.
3. **E1:** `density_noise_expected` (deterministic parity for D35) and `sample_density` with a gradient — exist in WIP.
4. **E2:** `ftcs_diffusion_1d`, `couette_startup_profile` (analytic V1 check of the FTCS run), `wall_shear_history`,
   `kinematic_viscosity`, `thermal_diffusivity` — exist in WIP.
5. **E3:** `process_path`, `path_heat_work_totals`, `process_heat_work`, `entropy_change_reversible`,
   `perfect_gas_entropy_change`, `irreversible_process` / `stirred_isochoric_process` / `free_expansion`,
   `isentropic_pressure` — exist in WIP.
6. **E4:** `brunt_vaisala_sq*`, `parcel_displacement`, `parcel_ode_atmosphere`, `parcel_temperature`,
   `potential_temperature`, `classify_stability`, `stability_timescale`, `synthetic_boundary_layer_column`,
   `seawater_density_linear` — exist in WIP; plus the new convention helper (item 2).
7. **E5:** `dimensional_matrix`, `rank_by_minors`, `minor_determinant`, `solve_exponents`, `pi_groups`,
   `groups_independent`, `group_value`, `rescale_units`, presets `PIPE`, `PENDULUM`, `SPHERE_DRAG`, `SCALE_HEIGHT`,
   `BLAST`, `RAYLEIGH` — exist in WIP. E5's JS needs exact rational arithmetic (a small Fraction class) for the
   null-space solve; promotion candidate for `viz_lib.js` once Ch. 4 §4.11 needs it.
8. **A-item figures:** `ch01.poiseuille_pressure_drop` (data generator for the C69 collapse figure, labelled "result
   derived in Ch. 8"), `core.thermo.tait_pressure` (C36), `core.statics.standard_atmosphere` /
   `atmosphere_from_temperature` / `isothermal_pressure` / `scale_height` (C40, C54), `core.units.dimensional_check`
   (C64) — exist in WIP. `scripts/ch01_*.py` already hold the drawings; keep notebook cells short by calling them.
9. **Budget (< 5 min on Colab CPU):** cache FTCS runs and `sample_density` sweeps with `FAST` grids (N = 200 in FAST,
   400 full; ≤ 40 slider steps); cap `parcel_ode_*` with a `solve_ivp` event for the unstable case; sympy checks (D19,
   D28, D10 residual) must avoid `simplify` on large expressions (use `expand` + `cancel`).
10. **Teaching traps to carry into names and docstrings:** reference pressure `p_ref` ≠ surface pressure `p0`; kmol vs
    mol in every gas constant; (1.35) holds as "same sign"; Clausius–Duhem with the actual δq; `capillary_rise` takes
    the book's α (complement of the contact angle).

**Optional — already exists in WIP (B/C items; keep, test at ≥ 1 evidence level if the notebook calls them, no new work)**
`shear_deformation_history` incl. Bingham/Maxwell (C02, N01) · `traction_components` (C03) · `number_density`,
`mean_molecular_spacing` (C04) · `maxwellian_velocities`, `molecular_pressure`, `wall_impact_pressure` (C05, D34) ·
`knudsen_number`, `mean_free_path_jennings`, `mean_free_path_air` (C07, N03) · `mass_fractions` (C09) ·
`fick_mass_flux`, `fourier_heat_flux` (C10, C11) · `viscosity_power_law`, `sutherland_viscosity`, `water_viscosity`
(C13) · `diffusion_time` (C14) · `surface_tension_water`, `laplace_pressure_jump` (C15, C16, N07) ·
`gauge_pressure`, `absolute_pressure` (C17) · `wedge_pressure_difference`, `wedge_face_forces` (C18) ·
`hydrostatic_pressure_uniform`, `buoyancy_force`, `net_pressure_force_on_box` (C21, D37) · `capillary_rise`,
`capillary_rise_deg`, `alpha_from_contact_angle` (C22) · `mean_molecular_speed`, `collision_time` (C23) ·
`perfect_gas_state`, `specific_volume`, `enthalpy`, `specific_heat_cp`, `specific_heat_cv`, `partial_derivative`
(C28–C31, N13) · `thermal_expansion_coefficient`, `water_density` (C37) · `molecular_gas_pressure`, `gas_constant`
(C38, C39) · `van_der_waals_*` (C41) · `cv_from_cp`, `gamma_from_cp` (C42, C43) · `isentropic_ratios`,
`perfect_gas_sound_speed`, `perfect_gas_expansion_coefficient` (C46–C48) · `potential_temperature_gradient`,
`potential_density`, `isentropic_density_gradient`, `ocean_potential_density_gradient`, `seawater_density_eos80`
(C56, C58, N21, C61) · `isothermal_density`, `linear_lapse_pressure` (C62) · `groups_independent` (C71) ·
`pythagoras_phi` (C74) · `blast_energy`, `blast_radius` (C75) · `rayleigh_scattering_ratio`, `wavelength_to_rgb`
(C76). No new functions are needed for B/C items.
