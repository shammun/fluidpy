# Explainer patterns — what worked, what failed, what to promote

Read by the curator, designer, viz-builders and viz-reviewer before any explainer work. Appended by the knowledge-keeper.
Explainer short names: ch01 E1 `continuum_averaging_volume`, E2 `viscosity_momentum_diffusion`, E3 `heat_work_paths`,
E4 `parcel_stability`, E5 `buckingham_pi_machine`; ch02 E1 `rotation_of_axes`, E2 `cauchy_traction_principal_axes`,
E3 `strain_vs_rotation_split`, E4 `gauss_flux_box`, E5 `stokes_circulation_loop` (all PASS round 2, 127 parity rows);
ch03 E1 `flow_lines_unsteady`, E2 `material_derivative_probe`, E3 `galilean_frames_cylinder`, E4
`fluid_element_deformation`, E5 `spin_and_principal_axes`, E6 `vortex_paddle_wheels`, E7 `reynolds_transport_cv` (all
PASS round 2, 145 parity rows, every page of every pager audited); ch04 E1 `control_volume_budgets`, E2
`stream_function_spacing`, E3 `newtonian_stress_lab`, E4 `navier_stokes_term_balance`, E5 `rotating_frame_coriolis`, E6
`which_bernoulli`, E7 `viscous_dissipation_heating`, E8 `boussinesq_buoyancy`, E9 `dynamic_similarity_models` (all PASS
round 2, 255 parity rows, 2 744 views, 30 derivations matched to the notebook by script); ch05 E1
`vortex_tubes_cannot_end`, E2 `vortex_pressure_funnel`, E3 `kelvin_material_loop`, E4 `baroclinic_torque`, E5
`vorticity_stretching_tilting`, E6 `biot_savart_filament`, E7 `vorticity_equation_rotating`, E8 `point_vortex_lab`, E9
`vortex_sheet_rollup` (all PASS round 2, 268 parity rows, 2 776 views, 23 derivations matched); ch06 E1 `superposition_sandbox`, E2 `cylinder_circulation_lift`, E3 `vortex_wall_images`, E4 `complex_potential_corners`, E5 `blasius_kutta_contour`, E6 `conformal_joukowski`, E7 `laplace_relaxation`, E8 `axial_singularity_bodies`, E9 `added_mass_sphere` (all PASS round 2, 291 parity rows, 2 864 views, 31 derivations matched). ch07 E1 `dispersion_relation`, E2 `particle_orbits`, E3 `capillary_gravity_waves`, E4 `seiche_standing_waves`, E5 `group_velocity_packets`, E6 `wave_rays_refraction`, E7 `hydraulic_jump`, E8 `two_layer_modes`, E9 `internal_wave_beams` (all PASS round 2, 222 parity rows of which 172 JS ↔ fluidpy, 3 312 views, 29 derivations matched). ch08 E1 `couette_poiseuille_backflow`, E2 `lubrication_scaling`, E3 `slider_bearing`, E4 `viscous_gravity_current`, E5 `stokes_first_problem`, E6 `similarity_exponents`, E7 `oscillating_plate`, E8 `stokes_sphere_flow`, E9 `stokes_drag_settling` (all PASS by round 3, 278 parity rows of which 207 JS ↔ fluidpy, 3 312 views, 26 derivations matched step by step against the notebook's `--dump`). ch09 E1 `bl_scaling_thicknesses`, E2 `blasius_similarity_collapse`, E3 `falkner_skan_family`, E4 `thwaites_marching`, E5 `cylinder_drag_crisis`, E6 `karman_street_stability`, E7 `free_jet_similarity`, E8 `wall_jet_invariant`, E9 `teacup_secondary_flow` (all PASS round 2 after a library fix, 248 selftest rows of which 178 JS ↔ fluidpy; 20 of 22 derivations match the final notebook — **D06 and D14 drifted after the lesson review**, see Failures). ch10 E1 `fd_stencil_order`, E2 `von_neumann_amplification`, E3 `upwind_cfl_advection`, E4 `cell_peclet_wiggles`, E5 `fem_hat_assembly`, E6 `mac_projection_staggered`, E7 `lid_driven_cavity`, E8 `mixed_fe_lbb` (all PASS round 2, 310 selftest rows of which 270 JS ↔ fluidpy, 20 derivations matched step for step against the notebook's `--dump`; no library change; heavy local numerics — dense/banded LU, lstsq, tridiagonal, generalised eigen — see ch10 candidates). ch11 E1 `normal_mode_growth`, E2 `kelvin_helmholtz_boundary`, E3 `benard_neutral_curve`, E4 `salt_fingers`, E5 `taylor_couette_onset`, E6 `richardson_shear_instability`, E7 `inviscid_shear_criteria`, E8 `orr_sommerfeld_neutral_curve`, E9 `lorenz_attractor` (all PASS round 2 — round 1 failed four on text/CSS items — 417 selftest rows of which 350 JS ↔ fluidpy, 2 824 views, 19 derivations / 207 steps matched against the notebook's `--dump`; no library JS change, two CSS rules added to `assets/viz_base.css` by the orchestrator; tables of fluidpy values embedded in E3, E5–E8 — see ch11 candidates and generators). **Rule since ch03: 5–10 explainers per chapter, as
many as the CORE ideas need** (`book.yaml → project.min/max_explainers_per_chapter`), and **every book equation cited
in tour, Explain, Derivation, quiz, notes or status text is written out in TeX next to its number** (`tools/eq_refs.py`
lists offenders; `ref:` labels next to shown TeX, metadata and selftest names are exempt).
ch12 E1 `reynolds_averaging_window`, E2 `correlation_and_spectrum`, E3 `reynolds_stress_parcels`, E4
`energy_cascade_spectrum`, E5 `turbulent_energy_budget`, E6 `turbulent_jet_similarity`, E7 `law_of_the_wall`, E8
`mixing_length_closure`, E9 `stratified_surface_layer`, E10 `taylor_dispersion` (all PASS round 2 — five failed round 1;
597 selftest rows of which 444 JS ↔ fluidpy, worst 8.0e-07; 24 derivation instances / 220 steps equal to the design; the
first chapter with ten explainers).

## Patterns that worked
| Pattern | Where proven | Why it works |
|---|---|---|
| Walkthrough step that sets the parameters itself, then shows one inline slider ("Drag U and watch…") | `templates/viz_example.html` step 3 | the reader sees the change before being asked to act |
| Live equation card in a walkthrough step (`eq:` + `live:`) | `templates/viz_example.html` step 4 | symbols become numbers the reader controls |
| Tracer particles + streamlines + colour field on one stage | `templates/viz_example.html` | motion shows what a static streamline plot hides |
| Draggable probe showing the local velocity vector | `templates/viz_example.html` | direct manipulation of "the field at a point" |
| Explain tab in numbered sections: what the colours are → each quantity from the controls (boxed results) → at the current time → "Reading the current setting" | `templates/viz_example.html` (after `forced_damped_vibrations.html`); ch01 E1–E5 (7–8 sections each) | every number on screen is accounted for, then interpreted for the regime the reader chose |
| Derivation step that moves the picture (`set`) and shows the line with the reader's numbers (`live`) | `templates/viz_example.html` derivations `energy` step 4, `omegad` step 2 | the symbols in the derivation become the numbers and curves on screen |
| Walkthrough step quoting one derivation step with an "All steps →" button | `templates/viz_example.html` tour step 4; ch01 E3 → D14 step 10, E4 → D18 step 11, E5 → D28 step 12 | the story stays short while the full derivation is one click away |
| **Two-convention status badge**: verdict word first (with emoji), then the library's exact criterion text in the chosen convention, then the other convention's bare relation ("🌊 stable · stable ⇔ dT/dz > Γa: −6.5 > −9.8 K/km · met: 6.5 < 9.8 K/km"); a selftest row rebuilds the badge from `cfg.status` and compares it with the fluidpy text | ch01 E4 (48/48 browser cases match `lapse_rate_stability`, 13 exact-text parity rows) | both sign conventions are on screen at every size and every step; the verdict never reads as the opposite of the inequality text |
| **Exact-text parity rows**: when the Python function returns a sentence (a criterion with numbers), the JS mirror's string is compared byte for byte, not only its numbers | ch01 E4 | catches rounding, minus-sign glyph and inequality-direction drift that a numeric row misses |
| **Separate y-title strip on narrow views**: below 420 px shift the plot rect 18 px right and draw the rotated y title in that strip | ch01 E4 (`vplot`) | the library clamps the left pad to 40 px there, and a rotated title otherwise covers 4-character ticks ("−150" read as "150") |
| **Put a hidden view's key number in a visible title**: on portrait phones the Kn strip is hidden, so the curve title reads "Kn(body) 0.064 → slip flow" and the step adds a Kn readout; on landscape phones the strip gets its own row | ch01 E1 step 6 | "turn your phone" hints are not an explanation; the number the step is about is always visible |
| **A trace step sets the global parameter to the sample's own value**: the step that follows one sample's arithmetic sets L to that box (21 nm), so badge, Eq. card and inspector all describe the same box | ch01 E1 step 4 | two different boxes on one screen (the round-1 bug) make the arithmetic untraceable |
| Steady-state ghost + fading trail of earlier profiles + a thin dashed exact series on top of the numerical profile | ch01 E2 | shows the transient, its end state and the correctness check in one picture |
| Modes that change only the diffusivity (momentum ν / heat κ / species κ_m) on the same stage | ch01 E2 | "same equation, different D" is the idea; switching modes proves it |
| Linked p–v and T–s diagrams with shaded work and heat areas (negative areas hatched), plus term bars q, w, Δe and Δs vs ∫δq/T | ch01 E3 | path functions change area while the state functions' bars stay put; irreversible mode separates Δs from ∫δq/T |
| Draw only what is defined: no moving state dot on a non-equilibrium (stirring, free expansion) route | ch01 E3 | a dot would claim equilibrium states that do not exist |
| On narrow views move the legend into the view title | ch01 E3 p–v view | frees the plot area; nothing overlaps the route |
| Labels that cross curves sit on solid card-coloured boxes (now `Viz.text(…, {bg: true})`) | ch01 E4 ΔT/Δρ labels | readable over any curve without shrinking the font |
| A derivation step that moves the picture to the special case it treats (double root: Γ = Γa, flat ζ(t)) | ch01 E4 D18 step 14 | the degenerate case is seen, not just stated |
| Unit-system switch (SI → cgs → imperial) on a table where every raw value changes and every Π value stays | ch01 E5 Equations tab | invariance is watched, not asserted |
| Exact rational arithmetic in JS (small `Frac` class mirroring `fractions.Fraction`) for the exponent solve | ch01 E5 | parity with `pi_groups` is exact, fractions print as ½, not 0.49999 |
| Presets at the classic cases: dry adiabatic (neutral), standard atmosphere, nocturnal inversion, superadiabatic surface layer, ocean thermocline, salt tank; fluids at their temperatures; six Π problems | ch01 E4, E2, E5 | each preset is a regime the text talks about |
| Pointer interaction through `stage.onPointer` (drag a gas state, click a diagram to probe) | ch01 E3 round 2 | no local capture-phase workaround; the engine routes `ev.viewId` |
| **Trace steps hold still**: `autoplay: false` in the app config, and `play: false` on every walkthrough step that quotes numbers; only a step that says "press ▶" plays; a derivation's result page gets `play: false` too | ch02 E1 (round-1 Must), E2 (D17 result page), E4 | a running transport overwrites the value the step just set, so the text and the picture disagree (1.866 in the text vs 2.215 on screen) |
| Long live equation lines split with `\begin{aligned}` (one relation per row) | all ch02 explainers (Equations tab, Explain) | fits the 1000×700 notebook frame without `equation-too-wide` |
| **One view row on landscape phones**: scoped `<style data-chapter>` hides `.viz-views > .viz-row:nth-child(2)` and the preset strip under `[data-layout="landscape"]` (and under `@media (max-height: 700–720px)` for short portrait phones), with `:not([data-tab="derive"])` so the Derivation tab keeps its own view; the hidden row's key number moves into the kept view's title | ch02 E3, E4, E5 | the engine has no per-row hide flag yet; the scoped CSS keeps every view ≥ 60 px without touching the library |
| **Shrink transport**: the transport runs a parameter σ = −log₁₀h, so "play" shrinks a box or loop geometrically toward a point; the limit view uses a log-h axis | ch02 E4 (`hOf = 10^−shrink`), E5 | a linear sweep spends almost all its time at large h; the limit is only visible on a log scale |
| **Shortened chip labels** for fields and presets ("(x², 0)", "b×x", "vortex") with the full name in the title/tooltip; `.viz-seg` chips wrap (`flex-wrap: wrap`) when there are six | ch02 E4, E5 | six field chips fit on a phone without overflow |
| **Diverging colormap with a white zero** (blue–white–orange for divergence, purple–white–teal for curl), limits symmetric ±max | ch02 E4 (`Viz.CMAPS.bwo`), E5 (`curlImage`) | the sign of ∇·u or ∇×u reads at a glance; "zero" is visibly empty |
| **Counter-example table**: four candidate triples, each with its (2.8) residual under the same C, plus a "fixed array (not a tensor)" preset whose (2.12) residual is in the status badge | ch02 E1 | the definition is taught by what fails it, by order one, not by round-off |
| **Passive/active table with the current row lit** and a toggle labelled "Wikipedia R" that computes Cx | ch02 E1 | the classic confusion is settled with a picture; `activeRotate(θ) = transformVector(−θ)` pinned by a parity row |
| Clickable matrix cells with an inspector that draws the two unit vectors and their angle (cos of the angle = the cell) | ch02 E1 | "what does one number C_ij measure?" is answered by pointing |
| A mirror toggle wired to a derivation step (det C = −1 while CᵀC − I stays 0; status turns amber) | ch02 E1, D02 step 8 | the step's special case (reflection) is seen |
| **Element + decomposition + term bars**: the removed half greyed, n, f and its blue normal / rose shear parts; σn = τ₁₁n₁² + (τ₁₂+τ₂₁)n₁n₂ + τ₂₂n₂² as bars that add exactly | ch02 E2 | the stress tensor is visible as "what pushes on a cut" |
| Three linked views on one angle φ: element · σn/τs curves with λ bounds and special-cut ticks · Mohr circle with the double angle | ch02 E2 | principal and max-shear cuts are the same event in three pictures |
| A 2-D sketch that shrinks with the derivation step (h = 1 → 0.5 → 0.25) while face arrows stay and the dashed volume term shrinks | ch02 E2 D05 | the orders-of-smallness argument is seen, not stated |
| Live numbers on the derivation steps where numbers help (6 of 15 in D17: discriminant, b·b = 0, τ′ = diag, c_k, σn = Σλc², shear bound) | ch02 E2 D17, E4 D22 (5 of 9) | a ★★★ derivation stays concrete |
| **Contrast toggle for an index convention**: non-symmetric τ draws τ·n in amber next to n·τ; Mohr says "not a circle" | ch02 E2 | the first-index rule becomes a visible difference (√3) |
| Adaptive status length (phone: verdict + one number; landscape: medium; desktop: full) | ch02 E2 | no three-line badges on phones |
| Mode-dependent axis names (stress: σn, τs [Pa]; strain rate: n·S·n, s·S·n [1/s]) | ch02 E2 (round-1 Should) | the same stage serves two tensors without mislabelled units |
| **One clock, three squares**: G leans, S (teal) stretches along ±45°, A (orange) spins; walkthrough G → split → S → A → all → blind → your turn | ch02 E3 | the cleanest story of ch02: the decomposition is watched, not asserted |
| Closed-form matrix exponential `expm2` (cosh / cos / shear-limit branches) matched to scipy `expm` at three regimes (1e-10) + tracer inspector showing e^{Mt} and x(t) under G, S, A | ch02 E3 | exact trajectories of linear flows in JS; reusable for Ch. 3 pathlines |
| Invariant pinned as a bar: S:A always 0 while S:S + 2S:A + A:A = G:G | ch02 E3 | Eq. (2.29) as a live invariant; name the factor 2 in the bars title |
| Conventions pinned by parity rows ("vorticity ω₃ = vector(R)", "R = 2A") and by one readout text ("Spin ½ω₃ = A₂₁") everywhere | ch02 E3 | the factor-2 trap cannot reappear silently |
| **Waterfall bars** left → bottom → right → top → ∮ next to ∬, with an interior-faces bar when tiled; face arrows coloured by sign (out orange, in blue) | ch02 E4 | Gauss' theorem as bookkeeping the reader can watch |
| Tiling with shared faces drawn as opposite pairs + a shared-face inspector ("+0.250 … from the left tile, −… from the right, sum = 0 — D25 step 8") | ch02 E4, E5 (edges) | the cancellation step of the proof is clickable |
| Limit view with three regimes in one picture: exact at every h (polynomial), slope-2 log–log inset (smooth), 1/h² blow-up (singular source) | ch02 E4 | "limit", "order" and "hypothesis" in one plot |
| **A theorem failing on purpose**: vortex core inside the loop → Γ = 2πK vs curl flux 0, bars "core inside: sides differ by 6.283", status ⚠️, `hypothesis_ok = False` parity row | ch02 E5, E4 (point source) | hypotheses are taught by breaking them |
| "Unrolled loop" view: u·t against arc length, shaded area = Γ, split into Ex. 2.6's four sides as (hatched when negative) term bars; loop-point inspector for one of 400 terms | ch02 E5 | a line integral becomes an area |
| Paddle wheel turning at ½(∇×u)₃ on a curl heat image, plus RK4 tracers; shear preset shows "curl without curved streamlines" | ch02 E5 | the curl's meaning in one glance |
| Staircase → circle refinement in a derivation (8×8 tiles, gap printed live) | ch02 E5 D26 step 9 | the tiling limit is measured, not claimed |
| Flip-n toggle that changes both sides of the theorem together, with the status saying so | ch02 E5 | orientation is a convention, the equality is not |
| **Measured ◇ against the formula**: an independent measurement (a float's sampled rate, tracked segment lengths and corner angles, a finite difference of ∫F) drawn as a diamond on the formula's bar or curve | ch03 E2 (float rate on the DT/Dt bar), E4 (stretching/closing/area dots on n·S·n, 2n₁·S·n₂, e^{t tr G}), E7 (FD of ∫F on the budget total) | the formula is checked in front of the reader, every frame; the chapter's signature move |
| **Answer-sheet view beside the experiment**: the closed-form result (Ex. 3.1's circles, line, centres) in its own view next to the simulated dye/particle/streamline, same clock; a ringed dot ties one simulated point to the inspector's arithmetic | ch03 E1 | the reader compares simulation and derivation without switching tabs |
| **Presets at the regime threshold** (loops / no loops at U₀ = ωξ_o; stalled front; halfway and overtaking observers) | ch03 E1, E2, E3 | each preset is a qualitatively different picture the text names |
| **Observer slider**: one parameter moves the observer continuously from one frame to another (lake → towed body); local and advective bars trade places while the total and a dashed body-frame reference stay put | ch03 E3 | frame dependence of the split vs frame independence of the sum is watched, not asserted |
| u = U + u′ vector triangle next to the flow | ch03 E3 | the Galilean transformation as a picture |
| **Single threads vs the flat pair average**: each material line's turning rate as a curve over θ, the pair average as a flat line at ½ω₃ | ch03 E5 | the "why ½ω, and why any pair" aha in one plot |
| **Co-rotating observer** (Ω = ½ω₃ stops the paddle wheel; ω′ = ω − 2Ω readout) | ch03 E5 | frame dependence of vorticity; seed for Ch. 13 relative vs absolute vorticity |
| **Mode-neutral slider name**: "Rate k" plus a meaning line next to G (shear → γ, solid body → ω₀, pure strain → S₁₁ = k/2) | ch03 E5 (round-1 Should) | one slider can drive several flows without mislabelling a quantity |
| Split view with a small table (strain part, rotation part, exact ellipse axes) for the probe point | ch03 E5 | (3.19) as numbers per row |
| **Orbit without turning**: paddle wheels carried round the line vortex keep their orientation while the draggable loop still reads Γ = 2πB | ch03 E6 | "going round is not spinning" is unmistakable |
| Log–log mean vorticity in a disc vs radius (slope −2 for the line vortex, flat for solid body) | ch03 E6 | the δ-function core (3.27) as a measurable slope |
| Derivation sector drawn on the phenomenon (D17's polar sector on the vortex, its four legs coloured) | ch03 E6 | the derivation's geometry is the picture |
| **Real-world table with the current row lit** (bathtub, tornado, tropical cyclone), sliders switch to log scales in metres and the status leads with Γ and σ | ch03 E6 (after the round-1 fix) | non-dimensional picture, dimensional numbers |
| **Swept band coloured by the sign of b·n** (advancing orange, retreating rose) + budget waterfall (volume term → swept in → swept out → total) landing on the measured ◇ | ch03 E7 | the signed wall term is visible; the budget closes on screen |
| **Dropped-term log–log panel**: the (3.32) term T4 against Δt with slope 2 next to the O(Δt) terms | ch03 E7 | "orders of smallness" becomes a slope the reader measures |
| Three geometries (1-D interval, deforming ellipse, growing cone) on one engine via modes; per-mode control hiding in Explore (`.rtt-off`) | ch03 E7 (Explore 15 → 11 pages at 360×640) | one theorem, three pictures, without a long control list |
| `\textstyle` integrals inside `aligned` derivation rows (display style only on desktop via a width class) | ch03 E7 | avoids the tall ∫ glyph being clipped at a pager slice break (library bug 3) |
| **One budget object, many scenes**: five scenes (wake, bore, jet on a plate, rocket, expanding balloon) all feed the same `cv_scenario`-shaped dict (faces, mass and momentum fluxes, body, surface, residual); the rule changes the box, never the law; residual bar at round-off | ch04 E1 | the reader learns that (4.5)/(4.17) is one tool, not five tricks |
| **"Drag the box until storage vanishes"**: the bore's CV speed b is a slider; at b = U the contents are steady and the unsteady wave becomes a steady budget | ch04 E1 | a genuine discovery interaction for "choose the frame that makes it steady" |
| Measured ◇ (finite-difference storage d/dt∫ρ dV) on the formula bar in the bore, rocket and balloon scenes (ch03 E7 pattern, again) | ch04 E1 | the budget closes in front of the reader |
| **Explain derives a hidden result from the same law**: the jet split Q(1 ± cos θ)/2 is derived in Explain §4 from along-plate momentum (after the review caught a 50/50 split) | ch04 E1 | the explainer teaches the method on a case the book does not do |
| **Spacing = speed bars**: amber bars of width Δn and height Δψ/Δn next to \|u\|, with a note that width × height of every bar is Δψ; a draggable, bendable gate whose flux stays ψ₂ − ψ₁ | ch04 E2 | "crowded contours = fast flow" is measured, not asserted |
| Continuous θ unwrapping across the source's branch cut when integrating ψ along a gate | ch04 E2 | a gate crossing the cut would otherwise jump by the source strength |
| Axisymmetric mode with 2πΔψ ring flux and the reading "contours crowd away from the axis" | ch04 E2 | the Stokes ψ's R-weighting is seen (it silently corrected the storyboard's "toward the axis") |
| **Source-coloured stress grid**: each τ_ij cell stacked by source (−pδ orange, 2μ dev S rose, μ_vS_mmδ purple), switching to λ/μ/γ colours during D09; "only μ + γ acts" watch | ch04 E3 | the constitutive law is read cell by cell; the γ = μ naming step is visible |
| Spinning cube mode: log–log α ∝ 1/h² with an end card; symmetry (4.25) as a divergence the reader watches | ch04 E3 | D08's orders-of-smallness argument becomes a slope |
| **Two-column "ρDu/Dt = forces" bar next to log-scaled term rows** and a regime word in the status (pressure ↔ viscous, local ↔ viscous, inertial, inviscid-idle) over 7 exact solutions | ch04 E4 | "which terms are awake" is read at a clicked point, not from a formula |
| The three viscous forms of (4.40) drawn as coincident arrows; a compressible toy field used only in D12 so the ∇(∇·u) term has something to show before it vanishes | ch04 E4 | an identity is seen as three arrows landing on one point |
| **Two observers on one clock** (inertial view with a "chalk mark" curve on the turning table; rotating view with the bent path) | ch04 E5 | the Coriolis force as the same motion seen twice |
| **Side toggle "forces (−) / accelerations (+)"** that flips every bar, arrow and side-panel label; labels carry the signed form ("−2Ω×u′ Coriolis"), pinned by 4 exact-text selftest rows | ch04 E5 (round-2 fix) | the (4.43) vs (4.45) sign trap is taught by a switch, and cannot silently come back |
| **The "why the 2" split arrow**: D15 draws the Coriolis arrow in two halves at the two steps that produce them (turning basis, d(Ω × x′)/dt) | ch04 E5 | the factor 2 clicks |
| Pole mode reproduces 9.34 km (exact) vs 9.45 km (small-angle Ωut²); f-plane high/low parcels; effective-gravity cross-section with an honest "×-exaggerated" note | ch04 E5 | the climate hooks, with numbers |
| **TRUTH vs your hypotheses decision table**: toggles for steady / inviscid / barotropic / irrotational; teal "holds here", amber "you claimed something this flow does not have"; the status names the valid Bernoulli form with its TeX | ch04 E6 | four Bernoulli equations sorted by their hypotheses, per flow |
| B along vs across streamlines with stacked term bars; whirlpool dip read from B(0); U-tube's amber ∂φ/∂t term makes B + ∫∂u/∂t ds flat | ch04 E6 | "constant along what?" answered by two probes |
| **A book slip shown with live numbers**: D26 step 6 evaluates both gauges (φ − ∫B dt′ → 0, the printed + → 2B) | ch04 E6 | the correction is checked, not just claimed |
| **Energy budget in = stored + out** with the share-of-work curve over κt/h²; "μ < 0 ⛔" preset turns ε and the entropy production negative and trips the status | ch04 E7 | the second law as a preset that breaks it |
| 2-D dissipation keeps the missing z deviator (−S_mm/3) so compressible parity rows are exact | ch04 E7 | a plane-flow shortcut that would silently be wrong |
| **Grey "dropped" bar with its ×αδT ratio** next to the kept buoyancy bar; buoyancy-off ghost circle; four-number validity chart with a five-setting table, current row lit | ch04 E8 | an approximation's size is on screen next to what it keeps |
| **Prototype and model drawn in their own l on one t* clock** with paired log bars and "✔ matched / ×5 / ÷125" badges | ch04 E9 | similarity (and the impossibility of matching Re and Fr) at a glance |
| 40 synthetic spheres of every size and fluid collapsing onto Morrison's C_D(Re) curve, with a toggle back to raw F vs U | ch04 E9 | data collapse is the definition of similarity |
| A forward-pointer mode clearly flagged as a preview (Ro = U/(2Ωl) map, "Ch. 13 preview" in tour, Equations, map and header) | ch04 E9 | later-chapter physics can appear without being taught here |
| Per-mode control **and readout** hiding via a scoped `.xxx-off` class toggled in `applyMode` (ch03 E7 pattern) | ch04 E1, E2, E3, E5, E6, E9 | keeps Explore short; now used by 7 explainers (engine request below) |
| Exact-text parity rows for rendered strings (status verdicts, signed term labels, `which_bernoulli_text`) | ch04 E5, E6 (after ch01 E4) | catches sign and wording drift a numeric row misses |
| Scenario tables mirrored byte for byte from fluidpy (`SCEN` = `BOUSSINESQ_SCENARIOS`, 20 rows parity-checked) | ch04 E8 | the explainer cannot drift from the notebook's numbers |
| **A contrast chip where the theorem fails**: four vorticity fields plus "broken ⚠" (∇·ω ≠ 0) whose Gauss sum −0.40 m²/s equals ∫∇·ω dV (selftest row) | ch05 E1 | "cannot end" is taught by the one field where it does end — and fails by exactly the violated hypothesis |
| **Key numbers in the scene caption in the bar colours** (lower / wall / upper / sum) so hiding the budget view on phones loses nothing | ch05 E1 | a cheaper alternative to "repeat the number in a visible title" (ch01 lesson) |
| **Needed = supplied bars with a measured slope**: "forces on the probe" (centripetal ρu²/r needed vs −∂p/∂r supplied) with the finite-difference slope of p(r) drawn beside the formula | ch05 E2 | "the pressure field really does the pushing" is measured, not asserted |
| **A deliberate miss**: the measured ◇ dΓ/dt lands on the (5.9) formula for a material loop and visibly misses it for a fixed loop ("the bars are for material loops") | ch05 E3 | the hypothesis (material loop) is taught by the mode that violates it |
| ✓/✗ hypothesis table (inviscid, barotropic, conservative, inertial) with the surviving (5.10) term as a bar, six flows on one clock | ch05 E3 | ch04 E6's decision table, now with the *consequence* (which term survives) drawn |
| **Finite element → point law with an R² gap panel**: the torque route ◇ (2M_G/I_G) rides the formula's sine; a log–log panel shows (route − formula) falling as R² | ch05 E4 | why a law derived on a finite disc is exact at a point, as a slope |
| "Order matters" step: θ = 270° flips the spin; quiz and a selftest invariant pin ∇ρ × ∇p (not ∇p × ∇ρ) | ch05 E4 | a sign convention becomes a visible reversal |
| **Colour-coded decomposition on a 3-D line**: purple stretching (along e_s), blue tilting (across e_n) arrows, bars split the same way; closed-form e^{Gt} presets + scaling-and-squaring `expm` for a custom G (parity 1e-11) | ch05 E5 | (5.32) read off the picture; the same colours in notebook, bars and Derivation |
| A transient that *settles* into a balance (Burgers core 4 → 2 mm) instead of starting at the steady state | ch05 E5 | "stretching balances diffusion" is watched happening |
| **"Build the sum" transport**: integrand pieces laid tip to tail onto the closed form, plus an unrolled integrand whose shaded area equals the result | ch05 E6 | an integral (Biot–Savart) becomes a running sum the reader can scrub (best screenshot `reports/viz/ch05/biot_savart_filament/desktop__tour-step4.png`) |
| **Book-slip toggle that flips the whole curve**: the printed −1/(4π) reverses u_θ(r) (rose) with a "⚠ printed" status; the second slip in rose beside the correct line in D11 step 5 | ch05 E6 | a sign slip is impossible to miss and its cancellation is visible |
| **Conserved quantity drawn flat in amber** against the changing ζ or Γ (column over a ridge, ring moved poleward), f-by-latitude table lit | ch05 E7 | (ζ + f)/h and Γ_a are seen staying put while everything else moves — the PV idea |
| Budget scenes that each isolate one new term (planetary in the stretched column, baroclinic in the lock, inertial in Burgers) | ch05 E7 | a 7-term equation is taught one term at a time |
| **Click-to-place, drag and delete point vortices with singular guards** (≥ 5 cm apart, clamped inside walls/bucket, halt at 4 mm, sub-stepping near close approaches) | ch05 E8 | free play without NaNs; reusable for any N-body or image system |
| Invariants plotted as Q/Q(0) flat lines + the measured orbit rate ◇ next to (Γ₁ + Γ₂)/2πh² in the title ("rate 0.637 · formula 0.637") | ch05 E8 | the integrator is checked live; the formula is checked by a measurement |
| A book-caption correction with a number and an arrow (fluid at G moves at −1.70 m/s) | ch05 E8 | "G is not a stagnation point" is seen |
| **Draggable circuit whose side bars trade while the total stays γ ds** | ch05 E9 | "circulation = vorticity inside" as bookkeeping (ch02 E4's waterfall on a sheet) |
| Convention toggle that flips only the sign (caption u₁ − u₂ vs text u₂ − u₁), with a status line | ch05 E9 | two conventions, one physics |
| Result page of a long derivation = the boxed result + the 5 key lines + "all steps" (D15: 25 → 5 pages at 360×640) | ch05 E7 (round 2) | the whole chain is one click away, the phone page stays short |
| Short plots (< 60 px): only the 0 tick and the current value at the dot ("Γ 0.02335 m²/s"); region labels in opposite corners on `bg: true` boxes ("ρ₂ heavy" top-left, "ρ₁ light" top-right) | ch05 E3 (round 2) | the ch01 narrow-view lesson, applied to labels as well as ticks |
| Collapse zero term rows below ~25 px per row ("zero: advective, baroclinic, diffusion") and hide a schematic row that carries no number on portrait | ch05 E7 (round 2) | seven term rows legible in ~120 px |
| Lazy trajectory cache keyed by the parameters (compute once; scrubbing reads the cache) | ch05 E3, E6 | smooth scrubbing of RK4/DOP853 runs without recomputing per frame |
| **Per-term bars "faint = size on the circle, solid = what survives the trip round"** for each power of z in the far-field Laurent series; only the 1/z bar survives | ch06 E5 (reviewer: best in chapter) | the residue theorem becomes visible bookkeeping; the (6.61) slip sits on a bar that vanishes anyway |
| **Result vs contour size, flat where the contour encloses the body, rose band where it cuts it** (force vs R) + teal per-arc shares that change while the amber total stays pinned | ch06 E5 | Cauchy's "move the contour freely" shown as a flat line and its failure as a band |
| **Term bars for the expanded formula with ◇ on the total** (the four pieces of (6.39), cross term amber, a dotted "without it" C_p curve) | ch06 E2 | the one term that makes lift is seen, not asserted |
| Keep the physically meaningful quantity when a slider moves another (Γ kept when U or a is dragged) — "L stays put when a changes" is discovered | ch06 E2 | a formula's missing variable taught by an invariance |
| **Both sign conventions in the status at every size** ("Γ_cw 2 ⇔ Γ_ccw −2 m²/s") + exact-text selftest rows for regime words (merged / free stagnation) | ch06 E2, E3, E6 | two conventions inside one chapter never confuse the picture (E1, E5 still lack it) |
| **Tracers that carry a physical label** (teal source fluid fills the half-body and never leaves it) | ch06 E1 | D07's mass-balance step made visible |
| Element kits with **ψ-bars tip-to-tail onto the velocity arrow**; click a bar to isolate that element | ch06 E1 | superposition read off one probe |
| A **"squeeze" button that animates a limit with a product held fixed** (2mε fixed) onto grey ghosts of the limit | ch06 E1 | P150's limit as motion; reusable for any singular limit |
| **Wrong-branch field drawn in rose** next to the right one (streamlines through numpy's principal root) + a ghost second root inside the circle | ch06 E6 | a branch cut becomes a visible failure, not a footnote |
| **Same particles in two planes** (ζ and z) and a cross with live β = α that doubles at the critical point | ch06 E6 | conformality and its exception on one clock |
| **Stencil stepping node by node** with the four-equation matrix (ψ now · b · exact · b − Aψ) and ghosts of the other two methods on the residual plot + predicted-vs-measured ρ table | ch06 E7 | an iterative method in slow motion (the fast.ai "pixels as parameters" move) |
| **"Moments converge, bars don't"**: Σk_nΔξ and the dipole against each target's exact moment while alternating strengths are shown, twin axes with a "cond > 10¹²: noise" band | ch06 E8 | honest numerics as a picture; fits every ill-conditioned inverse problem |
| Inspector that traces **one row of the matrix** with its two largest terms | ch06 E8 | collocation made concrete |
| **Speed part + acceleration part = total** in colour (teal + orange = black) with force bars showing the speed part at 0; added mass three ways (pressure, energy, formula) | ch06 E9 | d'Alembert and added mass in one picture; the "sideways" preset shows force following acceleration |
| **Sensor trace with faint whole + bold so far + markers at special times** (0, ±4πh²/Γ, ±4√3πh²/Γ) and a wall washed by p − p∞; "same-sign image" toggle leaks the wall in rose | ch06 E3 | the ch05 "deliberate miss" pattern applied to an image rule |
| Branch cut drawn on the plane and a status "⚠️ the step crosses the branch cut" | ch06 E4 | the reader learns where the formula stops being analytic |
| Honest refinement panel: order 4/3 printed with its r^{2/3} corner reason, not "≈ 2" | ch06 E7 | "state what converges and what doesn't" in an explainer |
| Per-state transport range (`syncRange`: the sweep slider's max follows method, problem and refinement) | ch06 E7 | works around the fixed transport range (library quirk) |
| **Dispersion diagram with ghost limits and shaded regime bands**: c(λ) with the deep √(g/k) (dashed) and shallow √(gH) (dotted) asymptotes and the 2 %/3 % bands shaded, a moving dot, beside the tank with orbits and the depth profile — three views on one clock | ch07 E1 | the regimes are areas on the curve, not words; the template for every Ch. 13 dispersion diagram |
| A **banner that carries a hidden view's key number** onto portrait phones (sea-bed pressure fraction when the profile view is hidden) | ch07 E1 | the ch01 lesson "a hidden view must not carry the step's number", solved without shrinking anything |
| **Linear vs exact path-line modes** with the measured drift per period next to the formula (7.85)/(7.86), and a dyed line leaning forward | ch07 E2 | the linearisation's failure (orbits that do not close) is seen and measured |
| **Restoring pushes drawn at the crest** (gravity blue, tension rose, with their pressures) + term bars that sum to c²/tanh kH with a c²_min floor; liquids as modes | ch07 E3 | two restoring forces add visibly; the minimum is a floor the bars never cross |
| **Shared "our computed numbers" constants** (`WM`, `WG`, `HM` from `cmin()`/`cgmin()` with σ = 0.07274, ρ = 998.2) feed every prose string, check answer and legend card at 4 s.f. | ch07 E3, E5 (round 2) | never type a number a function can print — and no exercise answer leaks at 3 s.f. (rule 9) |
| Right-going, left-going and sum drawn as ghosts; period-vs-mode plot with exact, shallow and deep curves | ch07 E4 | superposition and the mode ladder on one screen; parity rows prove right + left = standing |
| **Chord vs tangent on ω(k)** (beats → c_g), a dashed "group moving at c" ghost, an energy strip with conserved ∫E dx, and a printed-slip toggle (½Δω x) with its own parity row | ch07 E5 | the difference quotient becomes the derivative in front of the reader; the slip visibly freezes the envelope |
| **Pond step with a live "Right now" callout** (calm out to r = 33.1 cm at t = 1.86 s; edge ripples λ ≈ 4.354 cm) tied to the purple minimum on c_g(k) | ch07 E5 | a minimum of c_g becomes a place on the water |
| **Hamilton's ray equations integrated in the page** with the ω drift in the status, RK4 vs `ray_trace` parity rows (beach and island), Snell closed form vs RK4, Snell arithmetic in the inspector | ch07 E6 | rays are computed, not drawn; ω conservation is a live number |
| **One Code block per geometry, traced off-screen** (`codeRay()` runs the chosen geometry independently of the displayed mode) so every live comment holds a real number | ch07 E6 (round 2) | no "…" placeholder, no live value that depends on what happens to be displayed |
| Momentum and energy budgets as term bars; a **"forbidden" preset** (Fr₁ < 1 would create energy); bore and solitary modes with KdV invariants | ch07 E7 | the second law as a visible sign of a bar; reviewer "very good" |
| **Separate magnifications for surface and interface**, a vortex-sheet halo, ω(k) of both roots with their limits, a "which g′" note | ch07 E8 | a mode invisible at the surface (η/ζ ≈ −0.002) is still seen; reviewer "very good" |
| **Physical space + wavenumber space as two square views side by side on portrait phones** (tank cropped to \|x\|, \|z\| ≤ 5 m; K plane a square ±1.1\|K\| window; the third view hidden) | ch07 E9 (round 2) | replaces two squashed strips; K, c, c_g and the right angle stay legible at 360×640 |
| **Sign-safe preset "K up-left (k < 0)"** with parity rows per quadrant, ω/N-vs-θ view, live w-equation residual, E_k = E_p and F = c_g E invariants | ch07 E9 | the book's k > 0 assumption becomes a thing to try |
| **Per-mode control hiding by a data attribute + CSS** (`markMode(s, app)` sets `html[data-…-mode]`, CSS hides `[data-viz-key="param:H"]` etc., then `app.fit()`) | ch07 E5, E8, E9 | works around library quirk Q4 without touching the engine; Explore pages shrink |
| **Line + parabola = sum in colour** (rose Couette + orange Poiseuille = blue u(y)) with a black floor tangent that stands vertical exactly at backflow onset, an amber reversed layer with tracers and a bending dye line, and a Q waterfall (Q_c + Q_p = Q) | ch08 E1 | superposition and the onset criterion du/dy(0) = 0 are one picture; "backflow with Q > 0" is obvious from the bars |
| **Extremum marker driven by the derivation** (u_max dot at y*, live "5 mm + 2.5 mm = 7.5 mm" on D05 step 9) | ch08 E1 (round 2) | a derivation step has a picture to point at; the dense-grid-verified `state()` keys are finally shown |
| **Log term bars with a grey "dropped" band** and pressure-scale chips (P_a vs μUL/h²) showing that Λ only rescales p*; the cross-gap verdict "∂p*/∂y* = O(ε²)" moved into the title on phones; a "thick gap ✗" preset that breaks the theory | ch08 E2 | an ordering argument becomes bar lengths; the failing preset shows where the approximation ends |
| **Pad-frame picture** (u − U profiles, tracers, orange pressure arrows) with a black dp/dx tangent on the pressure hump, the **printed formula as a ghost with a measurable miss** (invariant row: slope error 1.398), golden-section optimum = fluidpy's Brent, Explain splitting C₁ into pressure part + drag part | ch08 E3 | the frame that makes the flow steady is the one drawn; a book slip is seen *and* measured |
| **Brute-force regime check in the explainer's own JS** (`recirc()` scanning both ends of the gap) independent of the mirrored fluidpy flag | ch08 E3 | it exposed that fluidpy's flag tested x = L only — an explainer can be the second, independent route |
| **Conservative finite-volume march with the volume in the status** (10⁻¹⁵), the front on log–log bending onto slope ⅕ beside diffusion's ½ ghost, two humps forgetting their start; parity with the implicit solver and Huppert's η_N from the Γ function | ch08 E4 | self-similarity is reached, not assumed; conservation is a live number |
| **Log clock (0.1–1000 s) + raw/rescaled mode** with Crank–Nicolson dots riding the erfc curve; "stop at T" and "second wall" presets that break the collapse (spread 1); δ₉₉ = 3.643 computed from erfcinv; ghosts of air and honey beside the reader's fluid | ch08 E5 | "profiles at all times are one curve" is shown *and* its failure when a scale is imposed |
| **Exponent plane** (n, m) with the 2m = 1 and n = m lines and a star; bracket-power chips that turn green when the powers of t match; a conserved-integral inset (∫ω dy ∝ t^0.1) that shows why exponents can match yet be wrong; four systems as modes on one control set; displayed spread = `similarity_collapse_error` | ch08 E6 | exponent matching as a search the reader performs; the second condition (a conserved quantity) is visible |
| **The discarded root drawn**: D24 step 7 plots the growing solution e^{+y/δ} as a rose dotted curve that explodes, so "B = 0 by boundedness" is seen; a "start from rest" CN march that lands on the periodic solution within a period; the book depth 4√(ν/ω) with its 5.9 % on both views; real presets (tide, molecular vs eddy ν, 1 kHz air, Ekman look-alike) | ch08 E7 | reviewer "excellent"; the reason a solution is rejected becomes a picture; "where is the initial condition?" is answered |
| **Three models on one control set** (Stokes, Oseen, ideal) with a dashed 1 % disturbance contour, an amber "inertia = friction" circle found by brentq on the exact ratio (41.5a at Re = 0.1) and the ½Re_a r/a asymptote, side-line speeds 0.9248U vs 1.0005U | ch08 E8 | an order-of-magnitude claim ("r ~ a/Re") becomes a measured circle with its prefactor |
| **Traction sweep**: orange pressure + rose friction arrows round the sphere summing to the same black x-push everywhere; running-drag curves reaching 2πμaU and 4πμaU; end card "2πμaU (⅓) + 4πμaU (⅔) = 6πμaU"; a real-particle table with the current row; ring-sum vs Gauss–Legendre parity row | ch08 E9 | reviewer "excellent"; a surface integral is built in front of the reader and its split is read off |
| **Derivation step counts checked against the notebook by script** (Playwright dump of `app.cfg.derivations` vs `build_chNN.py --dump`) | ch08 viz review | caught four explainers one step short (D05, D14, D19, D29) where the notebook had inserted moves |
| **Real vs ideal streamlines lifted by exactly δ*** on the plate view, with **shape chips at a fixed δ₉₉** (Blasius, linear, sine, cubic, exponential) so "same edge, different δ* and θ" is seen; D01 carries the scaled y-momentum with a live 1/Re; Explain §5 closes D04 numerically (ρU²θ = 2τ₀x) | ch09 E1 | three thicknesses become three integrals of one picture; a derivation's result is checked with the reader's numbers |
| **Raw ↔ rescaled collapse + a shooting overlay** ("guess too small / too large / converged" for f″(0)) and a truncation preset; τ₀(x) panel with "area = F_D" | ch09 E2 | similarity is the mode switch; the boundary-value problem becomes something the reader can miss by hand |
| **One dial through a fold**: Falkner–Skan n with the attached branch, the dashed reversed branch, a "reversed" preset and the fold marked; rose/orange term bars at the wall that stay equal (μu_yy = dp/dx) for every n; an independent RK4 shooting check printed in Explain §2; wall-curvature parabola on click | ch09 E3 | a bifurcation (saddle-node) is a place on a curve, and the wall relation is visibly an identity, not a coincidence |
| **Outer-flow picker → integral → criterion crossing**: θ² bar split into "∫U_e⁵ term" (teal) and "inlet memory" (rose); an inspector that redoes (9.50) → λ → l → τ₀ at the cursor; **two criteria labelled "book λ = −0.090" and "exact FS −0.0681", the book's first** in status, Explain and chips; diffuser closed form vs march in a selftest row | ch09 E4 (round 2) | one integral replaces a PDE in front of the reader; a book value and a computed exact value coexist without confusion |
| **Log-Re dial driving four linked things** (cartoon regime, separation angle, C_p(φ), C_D(Re)) with a highlighted regime table and honest "illustrative / rounded" badges; a **"set C_b by hand" check** that shows later separation alone does not lower the drag | ch09 E5 | a causal chain (transition → separation → wake pressure → drag) instead of four disconnected figures; a model's hidden assumption becomes a control |
| **Configuration + growth curve + complex-plane spectrum** on one control (b/a, stagger): eigenvalue dots reach the imaginary axis exactly at 0.2805; kick-and-watch growth; facing-rows preset (π/4 at every spacing) | ch09 E6 | the template for every normal-mode stability problem (Ch. 11, Ch. 13 baroclinic growth) |
| **Wrong-exponent toggle** that makes the conserved integral visibly change (ρ∫u²dy falls when δ ∝ x^{1/2}); momentum check in the status at three stations; a **level slider exposing the printed half-width as the 4 % point** | ch09 E7 | "the exponents are forced by the invariant" is shown by breaking them; a book slip becomes something the reader verifies by dragging |
| **Two invariant bars: one falls, one stays flat** — rose ρ∫u²dy drops 41 % downstream while the teal Ψ bar is flat to 4×10⁻⁷; printed-ODE toggle with hollow dots leaving the curves; gauge preset (C quartering with an unchanged dimensional profile) | ch09 E8 | reviewer's best "which quantity is conserved" picture; a redundant constant (gauge symmetry) is demonstrated, not asserted |
| **Force bars that add with height**: orange pressure push constant in z, rose centrifugal shrinking toward the floor, teal net = orange − rose; river-bend mode with the same balance; live Ekman-preview scales | ch09 E9 | the teacup's cause is one subtraction the reader sees at every height; the seed of the Ch. 13 Ekman-layer explainer |
| **Log–log error vs h on a transport that shrinks h** (tenfold every 0.8 s) with a rose round-off band, a ◆ best-h marker and a faint "so far" curve; **Taylor-term bars with exact fractions** (c₃ = 1/6, "cancels" labels) beside the measured error; a ×1/h zoom inset that keeps the stencil nodes visible below h = 10⁻² | ch10 E1 | "order p" becomes a slope the reader watches form, and the floor where round-off wins is part of the same curve |
| **Mode that turns stencil bars into scheme truncation bars**: the three terms of (10.17) add to the measured one-step residual (0.0317 vs 0.0317 1/s); a β = 1/6 preset makes the diffusive Δx² term cancel | ch10 E1 | a modified equation is decomposed like a budget; a special value is a preset the reader can test |
| **Complex-plane G(θ) + parameter-plane region drawn from the closed form *and* confirmed by a brute-force dot scan + a live 10⁻¹⁰ kick marched to a ×10⁶ cap with the predicted slope max\|G\|ⁿ dashed**; draggable (α, β) dot; classic-settings table; exact fluidpy reason strings pinned by parity rows | ch10 E2 | the causal chain point → curve leaves the unit circle → growth of one mode is on screen at once; the template for any growth factor (Ch. 11 dispersion relations, Ch. 13 model stability) |
| **x–t stencil view with the characteristic's foot ◆ and its interpolation weights** (×1.05 / ×−0.0526 when CFL fails) beside three schemes on one ring and an end-of-run card; \|G\| and phase-speed panels with the pulse's main wave marked; a "u = −1, printed" mode for the book's sign slip | ch10 E3 | CFL is geometry (interpolation vs extrapolation), and "first order pays in amplitude, second in phase" is two panels on one clock |
| **Discrete-root plot with a rose r < 0 band and a ± sign strip of rʲ** + a sweep transport that stops at the first wiggle (n = 19) + an inspector tracing T_j = (rʲ − 1)/(rⁿ − 1) + a dashed modified-equation ghost with its measured gap | ch10 E4 | "wiggles iff R_cell > 2" is read off a sign, not asserted; upwind "solves a different problem exactly" is a ghost on the dots |
| **Assembly as a transport**: each element's 2 × 2 block lands in the global matrix, outlined in amber; scatter-add table with the current element lit; **matrix-cell inspector listing its contributions** ("K₂₂ = k⁽²⁾₂₂ + k⁽³⁾₁₁ = 1.5 + 0.5 = 2"); FE = FD (10⁻¹⁷) as status and invariant row; natural-end mode with the learned slope (0.667 on 4 elements) | ch10 E5 | a sum sign becomes a sequence of visible additions; "natural = approached" is a number the reader sees |
| **An algorithm in slow motion on a transport** (u* → Poisson → push → uⁿ⁺¹) with cell-border divergence colours turning white; cell flux inspector with the wall face marked "wall" (never corrected); 3-cell pipe tiny example; SOR history vs direct and "divergence left = Δt × residual" invariant; **checkerboard mode with a row view and null-space bars** (collocated 4, staggered 1) | ch10 E6 | the projection deletes exactly the divergent part; the C-grid's reason to exist is a mode switch (template for Ch. 13 shallow water) |
| **A real solver running live in the page** (MAC 16²–32², banded LU factorised once, advanced by a background pump between frames) + cached finer runs from a parity-checked table + benchmark points + a Richardson panel + a benchmark table with the current row lit; interpretation "128² is further from Ghia than 64² — the comparison measures Ghia" | ch10 E7 | verification is watched, not told; the reader sees convergence saturate at the benchmark's own accuracy |
| **Counting that becomes an argument** (n_u vs n_p − 1; "at least 6 invisible patterns" at n = 2 = the 6 spurious modes found), an "invisible pattern" picture (max\|Bq\| = 4×10⁻¹⁶), β_h vs n with P1–P1's zeros on the floor; live dense FE ≤ 6 × 6 with a Cholesky + Jacobi generalised eigen-solver; one-element mode (six P2 shapes, 7-point rule) | ch10 E8 | an abstract inf–sup condition becomes a count and a curve that falls or stays flat |
| **A browser pass that changes state with the Explain tab open and compares with a fresh page at the same state** (25 cases), plus a scan of every config string and live output for control characters and `.katex-error` nodes | ch10 viz review | catches stale Explain text (the ch09 R2-M1 class) and lost backslashes that `shot.py` cannot see |
| **Three linked views on one state: the motion, the curve it is a point of, and the eigenvalue plane** — interface/shear/Bénard wave + σ(k) with the growing band shaded, ▲ cut-off, ◆ fastest mode and a ghost + the σ-plane with both roots; per-system controls hidden by a data attribute | ch11 E1 (also E2: wave + (k, ΔU) boundary + c-plane where the two roots collide ✕ and leave as a pair; E6; E7) | "which wavelength goes first" is visible only when one point on the curve and the motion it means move together |
| **Term bars under the square root that add to the discriminant, with the sign deciding** (gravity, surface tension, shear) | ch11 E2 | a stability criterion becomes "which bar wins"; presets at the minimum ◆ (6.70 m/s at 1.73 cm) and in the thermocline |
| **Table of fluidpy values with exact zero on the live neutral curve**: a small embedded table of `ch11.benard_growth_rate` (4 s.f.), interpolated, and **snapped to σ = 0 on the neutral curve the page computes live** (determinant) so table and curve can never disagree about the sign | ch11 E3 (replaced a one-mode estimate; 60 off-grid points within 7.7e-5; σ = 0 exactly at 0.999, 1, 1.001 × Ra(K)) | the sign (every verdict and band edge) is exact, the magnitude is the tested function's; parity rows at 2e-3 instead of rtol 0.25 |
| **Interpolation that ends on the exact neutral point**: the exact neutral point (k_n, c_i = 0, c_r known: k_n = 2, 1, √3/2 …) is appended to the tabulated spectrum as its last node and the table is interpolated cubically up to it, so growth reaches 0 exactly there and nowhere before | ch11 E7 (10 off-grid wavenumbers within 1.6e-5 of `rayleigh_eigs_contour`), E6 (growth map ending on J = k(1 − k)) | removes the false "stable" sliver that a bilinear table leaves next to a neutral curve |
| **Text-checksum parity rows for mirrored fluidpy strings**: the JS mirror of a verdict text is compared with the fluidpy string through a position-weighted sum of character codes computed on both sides, `rtol: 0` (plus regime-word checks on a grid of states) | ch11 E4 (seven checksum rows, one per regime text; regime word = `ch11.salt_finger_regime` on 1 120 states), E6 (nine rows for ±∞ / nan Richardson cases) | a mirrored sentence cannot go stale silently (the ch10 E2 failure); works inside shot.py's restricted `py:` builtins |
| **Status computed from the run, not from the parameters**: the badge reads what the integrated trajectory actually did ("still switches lobes at t = 50", "window too short (2.7 decades): no fitted line") | ch11 E9 | at r = 24 a stable fixed point coexists with a chaotic set — a parameter-only verdict would be wrong for half the starts |
| **Canvas orthographic 3-D projector with drag-to-orbit** (no three.js, no CDN): rotate, project, depth-sort, draw on the 2-D view canvas | ch11 E9 (works at 1366×768 and 390×844) | a 3-D orbit that survives offline and in a notebook frame; cheaper than `Viz.three` for line data |
| **Dimensionless transport parameter when the run length depends on the state**: the transport scrubs a fraction of the run (`run` or `tau` from 0 to 1) and the stage converts it to physical time with the current growth rate or period; where a clock in the flow's own unit is kept (E3 d²/κ, E6 L/U₀, E8 λ/U₀) the title prints the amplitude reached ("×43.7 (drawn capped)") | ch11 E1, E2, E4, E5 (0–1 parameter); E3, E6, E8 (scaled clock) | one transport serves presets whose natural durations differ by orders of magnitude (71.6 s thermocline vs 2.09 min cloud) |
| **A slider whose displayed value depends on another parameter**: the Ta slider prints "−175 (< 0)" past Rayleigh's line, the plane marks a hollow dot "slider value: not reachable here", and the status states the actual Ta | ch11 E5 (after review) | keeps one control honest where part of its range is physically unreachable; also the two-convention gradient slider (dT/dz and Γ = −dT/dz) of E4 |
| **Criterion band on a profile + map + three-way status**: Ri(z) on a log axis with the ¼ line and a rose band where N² − ¼U′² < 0; growth map with the exact neutral curve and a slice at the reader's J; status "🛡️ stable for every k" / "🌀 allowed and this wave grows" / "⚖️ allowed, not happening" | ch11 E6 (best of the chapter) | teaches a one-way theorem: "allowed" and "happening" are separate lights; a counter-example preset (R = 3, J = 0.4) shows Ri_centre > ¼ and growth |
| **Necessary conditions as bars evaluated on a computed mode** (Fjørtoft's two integrals; the two lobes of (11.84) cancelling) with a "Rayleigh ✓ Fjørtoft ✗" preset and a "both met, nothing grows" preset | ch11 E7 | the proof's integral becomes a picture; the reader sees why the theorem is only necessary |
| **Click a point on a neutral-curve map → the mode over the profile + production/dissipation bars adding to dE/dt = 2kc_iE**, with a Squire ghost point and honest labels ("(table)", "(table, interpolated)", "shape interpolated", "lower edge below k = 0.02, not resolved", "within the table's resolution of a neutral curve") | ch11 E8 | an energy budget explains the sign of the growth; the labels say exactly what is computed and what is not |
| **"How far to trust it" paragraph in Explain**: what a computed 0 means near a neutral curve, which numbers are interpolated, which finite-time fits are not long-time limits | ch11 E6, E8, E9 | readers stop reading resolution limits as physics |
| **A review-phase bug taught as a walkthrough step** ("top-heavy is not enough": Rs − Ra = 438 < 657.5, flips at −0.0076 K/m) | ch11 E4 | a boundary the code got wrong is exactly the boundary readers get wrong |
| **A DOM sweep for `.katex-error` nodes and control characters over every tab, tour step, derivation page and preset** (46–70 states per explainer), plus the regime word against fluidpy on a grid of reachable states | ch11 viz review (found E2's two KaTeX errors; 0 in round 2) | catches what `shot.py` passes (L3) |
| **Linked views on one cursor or one clock**: profile + paired term bars + whole-domain summary on one cursor height (E5); wind profile + temperature profile + budget bars on one cursor (E9); signal with a draggable lag + correlation + spectrum (E2); particles + log–log ⟨X²⟩ + D_T on one clock (E10) | ch12 E2, E5, E9, E10 | the reader moves one thing and sees every representation follow; phones show two views at a time and swap them by walkthrough step (E2) |
| **A draggable fit window on real (DNS) data** — drag the lower edge of the window and the fitted κ and B change live (model built with 0.41 returns 0.382 from y⁺ = 30; DNS goes 0.400 → 0.384 from 30 → 350); status "a quoted κ always belongs to a window" | ch12 E7 `law_of_the_wall` (the chapter's "aha"; replaced a wrong design expectation) | turns "constants are data" from a sentence into an experiment; the same stage fits a spectral slope or an e-folding scale in Ch. 13 |
| **Two-convention status badge in two lines**: line 1 verdict word · L_M · Rf; line 2 the criterion in the chosen convention with numbers, then the other convention's relation ("dT/dz > Γa: −6.5 > −9.8 K/km · Γ_met < Γ_d: 6.5 < 9.8 K/km"), identical at every size, pinned by exact-text rows; the verdict word comes from `["verdict"]`, never from the regime name | ch12 E9 `stratified_surface_layer` (151 rows, 134 py) | the project's lapse-rate rule is visible in every state, including phones where the temperature view is hidden |
| **Paired budgets that share one mirrored term, bars that sum to zero**: the production bar is drawn twice with opposite signs and joined by a dashed link; a "sink (remainder)" bar is labelled as the remainder; a whole-channel flow diagram (work → direct 44.7 % + via turbulence 55.3 %) with a table against Re_τ | ch12 E5 `turbulent_energy_budget` | "one term seen from two sides" in one glance; closing bars make a missing or wrong-signed term visible (→ Ch. 13 geostrophic and Ekman balances, energy cycle) |
| **A physical parameter as the transport**: the surface heat flux H sweeps through zero and L_M jumps through infinity while the wind crosses the neutral logarithm | ch12 E9 | a sign change is something to watch happen; the reader finds the singular point without being told |
| **Log-paced transport**: a 0…1 clock mapped by hand to time on a logarithmic scale, so that the ballistic and the diffusive regime are both seen in one run; first frame opens at t = 0.3 Λ_t with rays already drawn and "↺ replays from the release" | ch12 E10 `taylor_dispersion` | phenomena that span decades of time (dispersion, spin-up, adjustment) cannot be watched on a linear clock |
| **"Printed vs corrected" toggle for a book slip, with parity rows on both** (`printed=True` ↔ corrected): the printed +5/3 as a ghost that fails a units inspector and a code comment "the book's printed +5/3 (slip #1)"; the condition of the eddy diffusivity and the plume-width caption shown the right way round in the house wording; the printed exponential family worked with numbers beside a family that passes both tests | ch12 E4 (slip #1), E10 (slips #12, #13), E6 (slip #5) | a slip is taught by computing both versions; the house wording "⚠️ slip #k — the book prints …; the correct form is …" |
| **Constants computed from invariants instead of typed**: C₅ from the momentum-flux invariant (2.577 for ξ½ = 0.10), C₆ = 1/(C₅∫HF dξ) = 2.19 shown as the arithmetic, pinned by selftest rows; a scan for decimal literals in reader strings finds none | ch12 E6 `turbulent_jet_similarity` (after a stale typed 2.240) | a typed number goes stale the moment a default changes; a computed one cannot |
| **Live Code tabs that were actually run**: every Code block executed as displayed in three or four states by the reviewer (11 blocks in the five re-reviewed files); a block that builds its own runs says so ("N → ∞ here …" vs "the page's own 8 runs …") | ch12 viz review (found E1's `NameError`) | displayed code that raises is worse than no code |
| **Raw ↔ rescaled toggle plus a wrong-exponent slider that keeps every profile plausible and turns the invariant line rose** | ch12 E6 | the argument is checked by breaking it (same idea as ch09's collapse, with the invariant as the judge) |
| **One clicked sample → its arithmetic → thousands of samples with a ± 5 s.e. band around the running estimate**; every sampled number labelled "estimate" and held to 5 standard errors in the selftest | ch12 E3 `reynolds_stress_parcels` | statistics become a count the reader can follow; the band says when to stop trusting a digit |
| **Error view that decomposes an estimate's error into named parts** (noise · drift · leak) with ◆ formula beside ◇ measured | ch12 E1 `reynolds_averaging_window` | "how long must I average?" has a trade-off the reader can see; hook: a 30-year climate normal |
| **Real-case table with the current row highlighted + presets that are real flows** (kitchen mixer, wind tunnel, atmospheric boundary layer, ocean thermocline) and an inspector that checks units at the clicked point | ch12 E4 `energy_cascade_spectrum` | scale separation is a number in four familiar settings; the units inspector shows the printed slip failing |
| **A formula that ends with ✕ at the edge of its range, with the reason on the canvas** (log-linear wind stops at z = \|L_M\|/5) and the continuation labelled "a commonly used form", no source claimed | ch12 E9 | validity limits are drawn, not footnoted; honest about an unread citation |
| **Reviewer's `.katex-error` DOM probe + throwing-mode string probe** over default · every preset · every quiz state · both ends of every slider · every chip · every tab · every derivation page · every walkthrough step (218 state loads, 3 624 tab/page visits; 335 states, ≈ 4 240 fragments) | ch12 viz review (found E1 and E4; 0 in round 2) | catches what `shot.py` passes; the selector found the two round-1 failures, so the zeros can be trusted |

## Failures and fixes
| Problem | Where | Fix |
|---|---|---|
| App shell not 100 % high → panels clipped on phones but the audit saw no overflow | engine, before chapter 1 | `#app` height 100 %; audit fails `app-does-not-fill-window` |
| A walkthrough step with text + controls + readouts + equation overflowed on 360×640 | example step 5 | ≤ 2 extras per step |
| Explore intro + callouts pushed the controls off a phone screen | example | callouts optional by default; intro hidden at density 3 on phones |
| 7 tabs squeezed the title into a one-letter column on a 768 px tablet (no overflow, but views shrank below 60 px) | example, Explain + Derivation added | header falls back to short labels, then to a tab row of its own |
| Derivation tab on phones left 19 px views for a 3-view stage | example | phones keep one view (`derivation.view`) and hide the preset strip and transport on that tab |
| A derivation line "r = … = …" and a long live line overflowed the side panel in the 1000×700 notebook frame | example `omegad` step 4, `energy` step 4 | one relation per line; definitions move to their own step or into *why* |
| Term bars in a walkthrough step showed "–" | engine | term rows index their own item (bars from several blocks) and stale rows are dropped |
| **Every `onPointer` handler threw** (`ev.view = v` assigns a read-only `UIEvent` getter in strict mode); `tools/shot.py` never clicks, so all explainers passed | engine, found during ch01 E1–E3 | lib sets `ev.vizView` and shadows `view` with `defineProperty` (b45cbef); the viz-reviewer's 2 141-action click/drag pass is now the model; machinery TODO: shot.py must click/drag every view and fail on page errors |
| A view the walkthrough talks about is hidden on portrait phones ("turn the phone sideways") | ch01 E1 Kn strip (round-1 Must) | own row on landscape phones; key number in a visible title on portrait phones |
| One step showed two different sampling boxes (badge vs inspector) | ch01 E1 step 4 (round 1) | the step sets the global L to the traced box |
| Badge and inspector use different bases for the same quantity (mean count vs crest count: ±6.7 % vs 0.0607) | ch01 E1 (open) | compute both from the same count, or label the badge "at the mean density" |
| The other sign convention was not always on screen | ch01 E4 (round-1 Must) | badge carries both conventions at every size and step; wraps to two lines on phones |
| Unstable badge led with "stable ⇔ …" (the criterion text) | ch01 E4 (round-1 Must) | verdict word first, criterion second |
| Rotated y title covered minus signs of ticks on portrait phones | ch01 E4 (round-1 Must) | separate y-title strip (pattern above) |
| Tick labels collide with a marker in a 45 px ζ(t) view at 360×640 when the badge wraps | ch01 E4 (open) | below 60 px draw only the 0 tick or move ±ζ_max into the title |
| Library 3-significant-figure tick labels repeat on narrow ranges (301.9 … 302.3 K) | ch01 E4 | local `xTicks` with decimals from the tick step (promotion candidate below) |
| Derivation tab step counts drifted from the notebook after the notebook split steps (D14 10 vs 11, D18 13 vs 14, D28 12 vs 13) | ch01 E3, E4, E5 (round 1) | reviewer compares `app.cfg.derivations` step counts and move titles with `notebooks/build_chNN.py`; renumber walkthrough deep links in the same edit |
| *Why* text too long (46 and 53 words; limit 35) | ch01 E3 D10 step 7 (fixed), E5 D28 step 9 (open) | keep the rule in *why*; move the consequence into *in words* or *watch* |
| Raw TeX ("e^{±λt}") in a plain-text *why* | ch01 E4 D18 step 13 (open) | plain text uses "e^(±λt)"; TeX only inside `$…$` |
| WATCH line pointed at a view hidden on phones; "switch the medium" hint where the chips are hidden | ch01 E4 (round 1) | name the visible substitute ("the coloured bracket in T(z)", "use the thermocline preset") |
| Code tab used names it never defined (`eps`, `L_flow`, `T`, `p`) | ch01 E1 (round 1) | every name in a snippet is assigned in the snippet, with its live value in the comment |
| Legend overlapped the route on a phone p–v view | ch01 E3 (round 1) | legend in the view title |
| One-line end-of-run card ran past the right edge of a portrait-phone view and covered a state label | ch01 E3 (open) | `Viz.card` (wraps within the canvas) or a shorter card below 420 px |
| KaTeX warnings from Unicode superscripts ("m³") inside TeX | ch01 E5 (round 1) | map units to TeX powers (`m^{3}`) before rendering |
| SI / cgs / imperial columns dropped in the 1000×700 notebook frame | ch01 E5 (round 1) | narrower column labels + footer "columns: SI · cgs · imperial" |
| Clicking a starred (solution) header added it to the custom repeating set | ch01 E5 (open, optional) | ignore header clicks on the solution variable |
| Eq. card value lags the canvas while playing | ch01 E2 step 3 (open) | refresh the live equation from the same frame state as the canvas |
| `Viz.fmt`/`Viz.tnum` printed molecular masses (~5e-26 kg) as 0 | ch01 E1 | local `fmtK`/`T`; library now has `keepTiny` in both |
| **Transport autoplayed under trace steps**: the step text quoted 1.866 at θ = 30° while the plane showed θ = 55° | ch02 E1 (round-1 Must) | `autoplay: false` + `play: false` on every trace step (pattern above); the inspector re-renders or pauses on cell click |
| Labels collided at the default angles (arc label on shadow feet, θ label on x₁, primed labels on one face) | ch02 E1 (round-1 Must) | inspector arc label across the origin, θ label beyond the arrowhead, shadow labels offset perpendicular to their axis, primed normal/shear labels on different faces |
| Unprimed "τ₁₂ = 1.000" label crossed by its own rose arrow in tensor mode | ch02 E1 (open, cosmetic) | offset along the face or put it on a `bg: true` box |
| Third matrix-assembly line clipped at the canvas bottom on 360×640 and 844×390 | ch02 E1 | below ~180 px keep one fitted line (the rest is in the badge) |
| The bars gap label gave the wrong reason in the singular regime ("midpoint rule ∝ 1/n²" for a delta-function discrepancy) | ch02 E4 (round-1 Must) | regime-dependent label: "the delta at the origin carries the flux m — Gauss needs a smooth Q" |
| A derivation step's key number was cut off in a truncated phone title (and the limit view is hidden there) | ch02 E4 D21 step 4 (round-1 Must) | lead the narrow title with the key number; add a `live` line to the step |
| A walkthrough step said "the box shrinks" while the sweep first grew it (start at h = 2.5) | ch02 E4 | start the sweep at the stated value |
| Portrait field views tiny (square ~40 px, box ~70 px) or with a wide x-range (±5/±4) because equal-aspect plots in a short row | ch02 E3, E4, E5 (partly open) | `rows: [2, 1]` on portrait or clamp the range; let the heat image fill the width |
| End-of-run card covered axis labels of a panel | ch02 E3 | place it in the empty title band or a free corner |
| Single-panel mode used a third of a wide view | ch02 E3 | `pw = min(v.w, v.h·1.6)` |
| A step referred to a view hidden on portrait phones ("the dot on the flat curve") | ch02 E5, E1 step 7 | add the readouts to the step ("watch the Γ/A readout …"); repeated ch01 lesson |
| D26 step 9 showed a 13 % staircase-vs-circle gap without saying why | ch02 E5 | use the finest tiling and print the gap live |
| *Why* over 35 words (D01 36, D02 40, D17 39–42) | ch02 E1, E2 | move the consequence into *in words* (repeated ch01 lesson) |
| **Two builders sharing one scratchpad clobbered each other's `splice.py`** | ch02 viz phase (parallel builders) | every builder uses a private subfolder (`<scratchpad>/<slug>/`) for helper scripts |
| `shot.py` `py:` parity expressions have **no builtins** (`abs`, `len`, `float` fail) | ch02 all | index dicts/tuples down to one float in the expression; put `abs` on the JS side or return the magnitude from fluidpy |
| Library 3-significant-figure ticks and manual log ticks written twice | ch02 E4 `drawLimit`, E5 `drawBars` | promotion candidate `P.axes({xlog: true})` |
| `test_viz_library_inlined_and_template_lints` failed while builders were mid-build ("STALE viz/ch02/…") | ch02 verify phase | expected while the viz phase runs; re-run `tools/viz_inline.py --all` before the merge gate |
| **Pager: oversized items.** One paged item taller than its page (a walkthrough derivation quote, an Explain §0 paragraph, a "Right now" table, a code chunk with wrapped lines) was clipped by 24–209 px at 360×640; `shot.py` only saw page 1 of each pager, so all 10 ch01–ch02 explainers passed with it | review L1 (eq retrofit) | `Pager.layout` now slices such an item at natural breaks (paragraphs, text lines, list items, table rows, code rows, rows of an `aligned` formula; never through a line, a formula or a control), one slice per page with the cue "continued on the next page ›" and a "continued" rule at the top of the next page; a formula taller than a whole page sets `P.tall` and `fit()` moves the portrait split towards the text (`data-room="text"`). Equations on hidden pages are scaled too. The portrait stage's floor now includes `--viz-stage-need` (every visible view ≥ 60 px) and stays at that floor while a page would not fit (`data-paged`), so a page is never shorter than it was packed for (before, a short page let the stage grow and a long one was squeezed — the real cause of most L1 clips); live notes/inspector re-pack when their height changes; a text page < 110 px raises the density. `shot.py` pages through every page of every pager at every size and fails on clipping (`CLIP__*.png`). Builders: keep display lines short enough for 0.78× scaling at 360 px, and do not rely on paging to hold a 10-line formula |
| Views squeezed below 60 px on a phone by a full text page (parcel_stability Explain needed a local floor) | ch01 E4 | library: when a view is < 60 px in portrait, `fit()` retries with `data-room="stage"` (multi-view stage floor `max(300px, 62%)`) and keeps it if fewer problems remain |
| **Equation retrofit regressions** (ch01/ch02, "show every equation"): (a) an equation put into a heading/label that is also shown just below duplicated it and overflowed the page (gauss "Right now" +163 px); (b) TeX attached to the wrong number ((2.28)/(2.29) labelled with their consequences) | ch02 E4 R1, E3 R2 | (a) keep headings short, show the equation once in the body; (b) the TeX next to a number must be that numbered equation — the reviewer compares with the rendered page |
| Long live lines in the Equations tab (428 px product in a 332 px card) | ch01 E1 P2 | split into `\begin{aligned}` (product on line 1, result on line 2) or show only the formula and the result |
| **Presets overrode Γ and σ while the sliders still showed lab values** | ch03 E6 (round-1 Must) | a preset must set the slider values it uses (log sliders in metres for real vortices) and the status leads with them; a selftest row pins the status text |
| One slider ("γ") meant a shear rate, a rotation rate or a strain rate depending on the flow | ch03 E5 (round 1) | "Rate k" + `kMeaning()` line (pattern above) |
| Step 1's question fell onto page 2 at 360×640 | ch03 E4 (round 1) | step text ≤ 24 words when the stage is tall; check `phone__tour-step1.png` page 1 |
| Explore paged to 15 pages on phones (every mode's controls listed) | ch03 E7 (round 1) | per-mode control hiding via a class toggled by the chapter script |
| The "now" box of a derivation step on page 2 | ch03 E7 (round 1) | put the live line first on the step |
| Negative waterfall segment drawn pointing the same way as a positive one | ch03 E2 (open) | arrowhead or hatch for negative segments (E7's hatched "swept out" bar) |
| Legend strip covers the particle the step is about (portrait) | ch03 E3 (open) | legend into the view title on portrait (ch01 lesson, again) |
| Explain 21 pages at 360×640; D22 result 19 pages (whole 12-line chain) | ch03 E3, E7 (open) | merge sections on phones; collapse long result chains to first + last line on phones |
| Precomposed "ḃ" loses its dot at phone size; "36000 m" instead of "36 km"; labels cut inside ~45 px bars; touching y ticks in a ~70 px view; live numbers for a different particle than the view shows | ch03 E7, E6, E1 (open) | TeX `\dot b`; `Viz.fmt` with a unit switch; labels above narrow bars; draw only −1, 0, 1 below ~90 px; label which particle the live line is for |
| **Physics bug the tests missed**: the oblique jet split 50/50 for every θ (should be Q(1 ± cos θ)/2) — the scenario's *totals* closed, only the along-plate momentum was wrong | ch04 E1 via `cv_scenario("jet")` (derivation review Must-fix 1) | every scenario reports its full momentum balance per direction (`residual_momentum_along`), and a test derives the split by hand |
| **Side-panel term labels contradicted their signed values** ("Coriolis 2Ω×u′ = 3.35" on the force side, where the bar said −2Ω×u′) — the chapter's key sign trap | ch04 E5 (viz round-1 Must) | labels follow the side toggle; exact-text rows "term label sign = sign convention of the value" |
| **Garbled KaTeX**: `'…\rho … \tfrac … \partial…'` in single-quoted JS strings became CR, TAB and "p" (the file was written through a shell heredoc / single backslashes) | ch04 E7 D22 goal page and step 1 (viz round-1 Must) | double every backslash (or reuse an `EQ55` constant); scan every file for single-backslash TeX; write explainers with the Write/Edit tools, never a heredoc |
| Library `Viz.num.erf/erfc` accurate only to ~1.2e-7: too coarse for 1e-8 parity rows and for second differences of the Stokes erfc profile | ch04 E1 (`erfHP`), E4 (`erfcHi`) | local double-precision series/continued fraction (~1e-15); promotion candidate |
| A walkthrough card that needs 5–6 pages at 360×640 (inline display-wide equation + 50-word text + derive quote) | ch04 E1 s4, E2 s4, E4 s5–6, E5 s2/4/5, E6 s2, E7 s2/5, E8 s3, E9 s2 (Should, most still open) | move the equation into the step's `eq:` card, keep the text ≈ 30 words, let the `derive` quote carry the formula (E5, E7 went from 5–6 to ≤ 4 pages) |
| Explain pagers of 24–28 pages at 360×640 (fine at 390×844: 5–9) with a lonely §0 heading on page 1 | ch04 E3 (28), E8 (24), E4 (14) (Should, open) | two-sentence §0; open Explain at the section of the current mode; drop long lists on phones |
| A view title throttled during play shows a different time than the transport (126 ms vs 210 ms) | ch04 E4 step 4 (Should, open) | keep the time in the transport only, or refresh titles every frame |
| A stacked-bar label showed B + ∂φ/∂t while the view title said B; the total marker crossed its label | ch04 E6 U-tube (Should, open) | label the stack with the quantity it sums; offset in-bar labels from markers |
| Model panel empty while its label gives "wave length 64 m" (16 l, off the picture) | ch04 E9 ship (Should, open) | say "= 16 l — off the picture" |
| A selftest row tested a JS literal instead of an explainer function | ch04 E5 centrifugal potential (Should, open) | every parity row calls the explainer's own function (`centPot`) |
| A half arrow barely visible at t = 0.5 s on phones | ch04 E5 D15 step 6 (Should, open) | minimum drawn length 18 px, or set the step's time larger |
| Status "work in = heat out = 0.998 W/m²" while 0.2 % is still stored; wall fluxes without units | ch04 E7 (round 1, fixed) | "nearly steady: 99.8 % of the work in leaves as heat"; units on every number |
| The only `eq_refs` hit was a `<meta name="viz:concept">` string citing numbers bare | ch04 E6 (Should, open) | reword the meta without numbers, or teach `tools/eq_refs.py` to skip `<meta>` (it also flagged selftest names — see library/tool quirks) |
| **Rotated y title covered the minus sign of tick labels** ("−100" read as "100") at 390×844 — the ch01 E4 lesson, again | ch05 E2 (round-1 Must) | y title in its own strip on narrow views; **the engine should do this by default** (library candidate since ch01) |
| Two region labels drawn on top of each other ("heavy light 1000 kg/m³" suggested the heavy side was 1000) and overlapping ticks in a ~45 px Γ(t) plot | ch05 E3 (round-1 Must) | labels in opposite corners on `bg: true` boxes; below 60 px only the 0 tick + the value at the dot |
| **Seven term rows in ~120 px** at 390×844 — the step's key bar unreadable | ch05 E7 (round-1 Must) | hide the schematic row in budget mode on portrait; collapse zero rows; ≥ ~25 px per row |
| A derivation result page ("the whole chain", 15 steps) paged to 25 pages at 360×640 and 8 in the notebook frame | ch05 E7 D15 (fixed round 2) | result page = boxed result + 5 key lines + "all steps" |
| Walkthrough cards of 5–6 pages at 360×640 from a derivation quote + a 3-line code excerpt + readouts | ch05 E1 s5 (open), E2 s3/5/6 (fixed: 5 → 4) | one extra per step; quote a one-line derivation step |
| Round-off printed as a result ("Γ = −1.036×10⁻¹⁶") | ch05 E3 Helmholtz title (open) | print "Γ ≈ 0 (round-off 10⁻¹⁶)" (`termRO` in E3's bars already does this) |
| Combining diacritic in a canvas label ("ω̄" at 12 px renders "ω¯" with the bar offset) | ch05 E1 (open) | write "mean ω" or draw the bar |
| Isobar labels "+200 / +400 / +600 Pa" stack at 360×640 | ch05 E2 (open) | label only the middle isobar below 200 px, or stagger in r |
| Two readouts of one running sum computed from different counters (0.001843 vs 0.003868 m/s while playing) | ch05 E6 (open) | derive every "sum so far" from one `builtCount` |
| A parity row at rtol 1e-3 where JS and Python agree to 5e-11 | ch05 E9 (open) | set the tolerance near the achieved agreement (1e-8; 1e-4 for the N = 1000 Simpson row) — loose rows hide regressions |
| **Design claims copied would have shipped wrong physics**: the baroclinic angle convention reversed, D12's kernel variation 3a/\|x − x′\| (really 2a), D23 "short sides cancel" (really O(dn·ds), vanishing as dn → 0) | ch05 E4, E6, E9 builders (caught by computing) | builders recompute every design number and sign before coding it and report the correction so the design itself is fixed |
| Labels clipped at a view edge (Burgers "outflow u_z = αz" at the bottom; "arc length s along the filament [m]" at 844×390) | ch05 E5, E6 (open) | keep labels inside the plot rectangle; short titles below ~420 px |
| **Single-backslash TeX in JS strings**: `'$u_r=\frac1r\partial_\theta\psi$'` and `'$L=\rho U\Gamma$'` — `\f`, `\t`, `\r` became form-feed, tab and carriage return and KaTeX rendered "♀rac1r…" and "hoUΓ" | ch06 E2 (round-1 Must) | double every backslash; **candidate lint rule** `(?<!\\)\\(f\|t\|r\|b\|v)[a-z]` inside JS string literals (reviewer scan found only these two); the audit screenshots of derivation pages caught it |
| Bare equation numbers in `<meta name="viz:concept">` (published on the gallery card) — the chapter's only `eq_refs` hit | ch06 E6 (round-1 Must; same as ch04 E6) | meta text without numbers ("the Zhukhovsky map z = ζ + b²/ζ") |
| Unicode inside `\text{}` (m²/s, m³/s, ✚) → KaTeX "No character metrics" warnings at every size | ch06 E1 (open) | `\text{m}^{2}\text{/s}`; words instead of symbols in `\text{}` (ch01 lesson, again) |
| The walkthrough's named idea invisible at the chosen preset: with w = z² the x and iy quotients are exact at every h, so "two directions converge to one number" showed two coincident flat lines | ch06 E4 (open) | pick a preset where the effect is non-trivial (z³); check each step's picture shows the step's claim |
| **Check-question premises must be checked numerically**: a plate question in E6 and E5's "cut-cylinder" claim needed correcting by computation (orchestrator report) | ch06 E5, E6 | compute the answer of every quiz item with the fluidpy function before shipping |
| Computed ψ = 0 contour hooks onto the axis at both ends (marching squares through the R = 0 row) | ch06 E8 (open) | clip the contour to R ≥ 0 or stop at the axis row |
| Tracers drawn outside the plot rectangle (into tick labels and gutters) | ch06 E6 (fixed round 2) | draw particles inside `P.clip(...)` |
| A tolerance on one derivation page different from every other number on screen (10⁻⁸ vs 10⁻¹⁰ from r₀ = 3: 27/14 vs 34/19 sweeps) | ch06 E7 D23 step 7 (open) | one tolerance and one r₀ per explainer, stated once |
| Text that no longer matches the library after a fluidpy change (odd-N axial solver now falls back to least squares; Explain still said "the fit fails") | ch06 E8 (fixed round 2) | re-read explainer text after every fluidpy change the explainer mirrors |
| Explain §0 says a view is "hidden on phones" when it is not; colour names that differ from the drawn colours | ch06 E2, E9 (open) | generate hints from the view config (`hidePortrait`) and colour names from the palette |
| `play: true` on a question step or on a step that quotes numbers (the trace runs before the reader reads) | ch06 E3 s1, E9 s6 (open) | start paused; let a later step press play |
| 15–16 Explore pages of one control each at 360×640 | ch06 E1, E5 (open) | mark mode-irrelevant controls `optional` (library quirk Q4 again) |
| **Rule-9 near-miss**: a value we compute (minimum group speed) printed at 3 s.f. is an exercise's printed answer — in a tour step, a check question and a legend (two explainers) | ch07 E3, E5 (round-1 Must) | print 4 s.f. (17.76 cm/s at 4.354 cm) from the function; shared constants; grep explainer strings against `tests/book_values_chNN.json` before review |
| A walkthrough step about the island highlighted the beach code and showed a "…" placeholder in a live comment | ch07 E6 (round-1 Must) | its own Code block, traced independently of the displayed mode, with three island parity rows |
| Two stacked strips (tank, K plane) unreadable on portrait phones — K, c, c_g arrows and the right angle lost | ch07 E9 (round-1 Must) | two square views side by side at full stage height; hide the third view; crop the tank |
| **Intermittent 3 px clip** of a term list in a walkthrough card when KaTeX finished late (one `--quick` run in four) | ch07 E7 (open) | keep one line of slack in cards whose content sits near the page height; the pager should re-pack after KaTeX renders |
| A pager label "12/18" wrapping to two lines changed the header height on phones | ch07 E2, E6–E9 | chapter-local `.viz-pager-label { white-space: nowrap; min-width: 52px }` (promotion candidate for `viz_base.css`) |
| A truncated value in a derivation *live* line (0.969 for 0.96961) | ch07 E1 D09 step 4 (open) | round, never truncate (→ 0.970) |
| A status criterion looser than its words: "lands on the lee side" tested `x_end < 0`, so a ray grazing the flank counted as sheltered | ch07 E6 (open) | a strict criterion (x < −R/2 or polar angle > 120° from the incoming direction), or say "flank / lee" |
| A computed value near a book-rounded one prints the book's digits at 3 s.f. (the live c label at λ = λ_m reproduces the rounded (7.59)) | ch07 E3 (open; computed, not a breach) | 4 s.f. when \|λ − λ_m\|/λ_m < 1 % |
| Phone page counts: Explain 16–19 pages, Explore 10–12 pages, a tour step of 7 pages | ch07 E2, E7, E8 (open) | merge Explain sections on portrait; trim a step's extras on phones; mark mode-irrelevant controls hidden (markMode) |
| Explore intro names a view hidden on portrait ("click the period plot") | ch07 E4 (open) | phone-specific intro (the ch01 lesson again: never send phone readers to hidden controls) |
| A stray ray label "0.25" without "ω/N =", a phone K plane without a k-axis title, a pond x-axis without "cm" | ch07 E9, E5 (open) | on views < 420 px drop the label or write the symbol; units in the view title |
| A hard-coded live value in the Code tab (`lam0: '0.0435'`) next to computed ones | ch07 E5 (open) | compute it from the same constant (`CGM.lam`) |
| Derivation step titles that relied on a bare equation number (eq_refs hits) | ch07 E9 (round 2) | rename the step by its move ("Apply ∇_H² to the second link"); `ref:` badges are accepted |
| **Explainer Derivation tabs built from design Part F drifted from the notebook**: the notebook inserted moves the design skipped (D05 extremum steps 8–9, D14 step 12 "substitute C₁ and C₂", D19 step 10 "insert A and B", D29 step 4 split) — four explainers one step short; one kept a retracted *why* (D14 step 15) | ch08 E1, E3, E5, E8 (round-1 Must) | **the notebook (`build_chNN.py --dump`) is the reference for step counts and titles**; build Derivation tabs after the notebook's derivations are final, re-dump after every lesson-review round, renumber every walkthrough `derive: {id, step}` link in the same edit |
| **An explainer mirrored a fluidpy flag that was itself wrong** (slider recirculation at x = L only) and no parity row covered α < −½ although the slider reached −0.9 | ch08 E3 (round-1 Must) | parity rows at every regime of the slider's range, not only the default; when fluidpy changes, re-read the explainer's mirror, Explain text, quiz and code comments ("inlet backflow" → "recirculation at the wide end") |
| **Design presets that were already past the regime they meant to show**: an "8 Pa/m" favourable/adverse preset already beyond zero net flow (6μU/h²); ice at μ = 10⁴ instead of 10¹³ Pa s; an inconsistent oil density (ν = 10⁻⁴ vs ρ = 870, μ = 0.05); a "Re = 2" bead that is not creeping flow | ch08 E1, E4, E2, E9 (found by the builders) | compute every preset's regime numbers (threshold multiples, Re, Λ) with the fluidpy function before building; a preset's label states its regime and a parity row pins it |
| **Intermittent phone-land clip** of a "Vorticity content" equation card: a display-size ∫ in a *live* row rendered after the pager had measured it (failed one full run in three) | ch08 E5 (round-2 Must) | **text-size integrals (`\textstyle`) in live rows**; move one-line notes into *why*; run shot twice to confirm (the ch07 E7 KaTeX-late lesson again) |
| Full-HD Equations clip of a tall Γ-function card whose live line the pager's height estimate missed; raw TeX "\rho g/3\mu" in a symbols list | ch08 E4 (round-1 Must) | split the card or move the formula into `note`; symbol descriptions are plain text or wrapped in `$…$` |
| Caret text in titles and axis labels ("t^(1/5)", "h τ^m") | ch08 E4 (fixed round 2) | superscript glyphs for integer powers (t², 10⁻³). **Never fake a superscript slash with 'ᐟ' (U+141F, a Canadian-syllabics letter)**: it rendered misaligned on the site's chapter list (user report, 2026-09-29). Fractional powers: `$t^{1/5}$` in any HTML text (mathify/KaTeX) and in `viz:title` / `viz:summary` (the site, gallery, notebook toolbar and page render `$…$` via `fluidpy.core.embed.title_html`; `title_text` gives the plain form); a raised small-font exponent on canvas, or U+2044 '⁄' as a fallback; plain `t^(1/5)` only in `<title>` |
| **Library pager places a tall "now" box as a unit** and leaves ~450 px of blank space on phone-tall derivation pages ("we had" + move title on page 1, "now" on page 2) | ch08 E3 D14 s14, E8 D29 s5 (open; library) | library: allow the "now" box to slice at `aligned` rows, or a derivation-result "chain" option; locally: shorten "we had" (right-hand side only on phones) |
| Ghost labels at fixed offsets sat on the wrong lines ("water" above the air line) | ch08 E5 (fixed round 2) | anchor each label at its own line's end |
| Curves identified by colour only on a portrait view; literal "y/δ_e" underscore; a narrow-view "Q 3.33" label over the bars; phone-land term bars that lose their names | ch08 E9, E7, E1, E2 (open) | legend words in the view title on portrait; Unicode subscripts (δₑ); drop value labels below ~200 px; short term names ("inert.", "press.", "fric.") |
| A minimum labelled "optimum" when the sign flips it into the strongest suction (U < 0) | ch08 E3 (open) | regime-dependent labels for special points |
| Shot runs under memory pressure fail or time out when nine explainers are audited in parallel | ch08 viz review | run `tools/shot.py` serially (one file at a time) when memory is short; re-run a failing size alone before calling it a physics failure |
| **Lost backslashes became control characters again**: a plain single-quoted JS string `'…-\frac1\rho\frac{dp}{dx}+\nu u_{yy}'` turned `\f` → form feed, `\r` → CR, `\n` → newline (D12 step 1 read "♠rac1 ho♠rac…"); `'$\eta_i'` → "eta_i"; and **inside a `String.raw` hint the file itself contained a literal form feed and tab** (`$St=\x0crac{…}`, `2\pi\x09imes0.2`) — introduced when the file was written | ch09 E3, E5 (round-1 Must; ch04 E7 and ch06 E2 before) | double every backslash in plain strings or use `String.raw` consistently; scan the written file for bytes < 0x20 other than `\n` (a script found none after the fix). **Candidate: `tools/viz_lint.py` fails on any control character other than newline** (the fourth chapter with this bug) |
| **Explain went stale when the transport parameter was dragged by hand** (library, R2-M1): `app.set` treated a change of only the transport parameter as "time only" and refreshed just `explain.live`, so Explain §1–§4 showed x = 0.0757 m while the view said 0.087 m (thwaites_marching) and n = 1.89 while the slider read 3.56 (falkner_skan_family); it looked right while playing (the loop rebuilds every 110 ms) | ch09 E1, E2, E3, E4, E7, E8 (transport = x, n, a cursor); not E5, E6, E9 (a real clock) | **fixed in `assets/viz_lib.js`** (aa34f17): a user-driven `set()` always calls `updateDynamic(true, true)`; `viz_inline.py --all` (75 files), manual keyboard-slider check in Edge, `test_machinery.py` 9 passed, template and ch09 `--quick` PASS. `shot.py` cannot see it (it sets state, it does not drag): **candidate machinery test that moves a transport slider and asserts the Explain text changes** |
| **Example 9.2 preset led with the non-book criterion** (x/L = 0.126 at the exact-FS λ = −0.0681) while the book's example uses λ = −0.090 (0.158); the −0.09 line was labelled only "Thwaites fit" | ch09 E4 (round-1 Must) | preset sets `closure: 'white'`; status "book criterion λ = −0.090: x/L = 0.158 · exact Falkner–Skan λ = −0.0681: x/L = 0.126"; chips "book / exact FS"; 4 parity rows against `example_9_2()["criteria"]`; cylinder preset 103.1° first. Rule: **when the book's value and the exact value differ, show both everywhere, the book's first** |
| Label "fold n = −0.0904" right-aligned onto the λ = −0.09 line — exactly the n ≈ λ ≈ −0.09 coincidence a reader confuses | ch09 E4 (round 1, fixed) | "fold n = −0.0904 (λ = −0.0681)" beside the red line, "book λ = −0.090" at the amber one |
| Two boxed `aligned` results 1–2 px too wide at 768×1024 (tablet) in Explain §6–§7 | ch09 E8 (round-1 Must) | shorter rows (tail formula and "(correct: …)" on their own rows) or `W.line` instead of a box |
| Bare "(9.83)" in a derivation result, check and live line (eq_refs hit) | ch09 E8 D21 (round-1 Must) | write the implicit solution (or its short form) next to the number |
| `\&` inside `aligned` rendered a literal "&" ("f → λf(λη)&f∞ → λf∞"); dimension slip K = (ṁ₁/ρ)²/ν (m²/s) where K and C are m^{3/2}/s; "centre speed" 183 m/s for a velocity **scale** whose peak is 14.4 m/s | ch09 E8 (round 1, fixed) | `\\&`; K = (ṁ₁/ρ)²/(ν x₁^{1/2}) with x₁ = 1 m; call u₀ the "velocity scale" in wall mode |
| **Derivation tabs drifted from the notebook after the viz review**: the lesson review split D06 step 6 into five moves (notebook 14 steps; explainer 10 with an extra first step) and the notebook's D14 has 16 steps (explainer 15: "Conjugate the lattice sums" missing) — the viz round-2 audit compared against the pre-split state | ch09 E2, E6 (found by the knowledge pass, open) | the ch08 rule "notebook = reference" needs a **re-check after the notebook phase, not only in the viz phase**: re-dump `build_chNN.py --dump` after every lesson-review round and diff against `app.cfg.derivations` (script it in `embed_check`) |
| A number copied from a Python table checked at rtol 1e-4 while it agrees to 1e-12 (`m_sep`) | ch09 E3 (Should, fixed) | tighten to near the achieved agreement (1e-9) so the row can catch a stale table |
| A mode-blind Explain §0 (sphere mode still describing the cylinder's C_p view); a default (U = 10 m/s) that made the street preset 13 kHz for a 0.15 mm wire | ch09 E5 (§0 open; default fixed to U = 1 m/s: 266 Hz "audible") | branch §0 on the mode; pick defaults whose numbers mean something to a first reader |
| A convention differing between a derivation and its picture (+Γ upper row in D14 step 1, −Γ in the view and code) | ch09 E6 (open) | one sign convention across derivation, picture, Explain and code; if the derivation must differ, say so on the step |
| A curve drawn on another quantity's axis without its own scale (blue u(z) on the force axis) and unlabelled grey dashed curves | ch09 E9 (open) | every curve gets a legend word and a scale (top axis u/u_e) or is removed |
| The rescaled collapse drawn as one function three times (the collapse assumed, not shown) | ch09 E2 (open) | plot each station's own u(y) samples divided by its own δ(x_k) |
| **An explainer mirroring a fluidpy string went stale when fluidpy changed**: `stability_verdict` gained "\|C\| + 2b > 1" for upwind with diffusion; E2's JS `verdict()` still said "2β > 1, the zigzag grows" at α = β = 0.3 (2β = 0.6) and its Code-tab comment no longer matched the Python output; E3's copy was dead code with the same stale text | ch10 E2 (round-1 Must), E3 (Should) | mirror the full function (every scheme, β > 0), extend the pretty-printer, and **pin every reason string with an exact-text parity row** (`rtol: 0`) so the next change fails the audit |
| **Explain template assumed signs**: `1-${t3(n.C)}(1-${t3(cs)}+${t3(sn)}\mathrm i)` printed "= 1 − 0.6(1 − −1 + 0i)" at θ = π, and round-off sin θ gave "−0 i" | ch10 E2 (Should, fixed in rounds 2–3) | sign-aware complex formatter (`ctex`), print the real number alone when \|Im\| < 1e-12 |
| **A derivation line with `\quad(10.88)` was 17 px too wide in the tablet walkthrough mini-card** | ch10 E4 (round-1 Must) | two `aligned` rows with the number on the second |
| **A bare equation number inside a template string** ("the parent triangle of (10.184)") — `tools/eq_refs.py` missed it because other `$…$` followed in the same string | ch10 E8 (round-1 Must) | write the map next to its number; eq_refs should scan template strings per sentence |
| **State-dependent phone clip**: the audit reached the Equations tab with the β = 100 BTCS preset still set, so the (10.26)–(10.27) live lines were long ("1.59 × 10⁵ + 0 = …") and page 3 clipped the (10.28) card by 26 px — reproducible, but only after that preset | ch10 E2 (round-2 Must, fixed by the orchestrator) | shorten the card (one relation per line, live text on one line) or mark long live lines `optional`; **audit Equations/Explain at the extreme presets, not only the default**; shot should reset presets per tab or run a preset sweep |
| **Intermittent phone clips under machine load**: a 1 px phone-land Explore "Right now" note (FAIL in one of two runs) and an 11 px phone clip of the (10.88) large-R card that appeared only during the loaded merge-gate run | ch10 E3, E4 (fixed by the orchestrator) | trim a clause / shorten the card and its live line and trim neighbours (leave margin); **run shot 3× in a row and once under load** before calling a clip real or fixed |
| **Illegible picture on phones at the largest n** (16² staggered grid: 1–2 px face arrows, divergence borders a speckle) | ch10 E6 (Should, open) | cap n per layout (12 on portrait) or draw every other arrow — another case for `views[i].hideOn`/per-layout parameter limits |
| **Three clocks disagreed on a phone while a live solver spun up** (transport 0.6 L/U, view title t 0.26, readout t 0.51) | ch10 E7 (Should, open) | status "computing … t = 0.26" whenever the picture lags the transport by more than a snapshot, or hold the transport until the run catches up |
| **A quiz question not answerable in the UI** ("why is the central stencil exact for x² but not x³?" with only sin, eˣ, e^(−x²) offered) | ch10 E1 (Should, open) | add the function (a polynomial chip) or reword to an offered one |
| **A legend missing a series** (MacCormack ◆ in the convergence panel and Explain §0) | ch10 E7 (Should, open) | every marker style in the view title's legend and in Explain §0 |
| **`py:` parity rows cannot call ndarray methods** (`.min()`, `.sum()`) under shot.py's restricted builtins | ch10 builders | use `np.min(...)`, `np.add.reduce(...)`; tooling: expose them |
| **Code-tab `{{placeholders}}` share one namespace across code blocks** (a name reused in two blocks shows the last value in both) | ch10 builders (library quirk) | unique placeholder names per block; library: per-block namespaces |
| **Shell heredocs still mangle backslashes** (`\approx` lost in an explainer string written by a heredoc; `\to` → TAB + "o" in the notebook) | ch10 builders + notebook | Write/Edit only; scan every written file for bytes < 0x20 (the builder's `self_check_ctrl`); **viz_lint rule still missing** (five chapters) |
| **Parallel builders collided in one shared scratchpad folder** | ch10 viz phase | each agent writes its own scratchpad subfolder |
| **A one-mode estimate of σ presented as the rigid-wall growth rate was 5–35 % off** (+19.6 % at K = 3.1163, Ra = 2500; +35.1 % at K = 2, Ra = 2500), while the text said "about 10–20 %" and the parity row carried rtol 0.25 | ch11 E3 (round-1 Must) | embed a table of the tested fluidpy function and interpolate; never ship a parity row with a tolerance wider than the claim in the text |
| **"≈" bridges hide solver defects**: a status that says "≈ neutral" or a row with a loose tolerance near a neutral point looked like honest rounding, but the underlying fluidpy tables had false zeros (Rayleigh jet k = 1.8, 1.9; TG tongue at k ≥ 0.65) | ch11 E7, E6 builders (found the defects) | when the page disagrees with a table near a boundary, report it to the implementer instead of widening the tolerance or adding "≈"; the builders were the accidental testers |
| **Bilinear interpolation across a neutral curve** puts the zero contour in the wrong place and gives the wrong sign in a sliver (and 65 % error in a c_i of size 1e-4 inside the Bickley gap) | ch11 E6, E7 (fixed), E8 (accepted, labelled) | interpolate in a variable that vanishes on the known neutral curve, or snap to the live root-found curve; where neither is possible, say "within the table's resolution of a neutral curve" and do not claim a sign |
| **A reference grid that switches mode family cannot be interpolated**: `reference/ch11/os_grid_<flow>.csv` holds the least-damped discrete mode, which jumps between wall and centre modes outside the neutral curve | ch11 E8 | generate the explainer table by **continuation of one mode** from the most unstable node (`tools/viz_tables/ch11/orr_sommerfeld_neutral_curve_table.py`); for the jet's two bands (two families) tabulate the least-damped discrete mode on a finer log-k grid and say so |
| **A normal quoted JS string with single-backslash TeX** (`'$c=\bar U\pm[\,\cdot\,]^{1/2}$'`: `\b` became a backspace) rendered two red KaTeX errors on a Derivation result page; `shot.py` passed | ch11 E2 (round-1 Must) | raw template strings (or doubled backslashes) for every TeX string; shot.py must fail on `.katex-error` (L3, open) |
| **Coloured legend words in view titles printed in plain text colour** (`.viz-view-title b` outranked `.c-teal` …): curves had no other label on the stage | ch11 E1, E4, E8, E9 (E4, E9 Must) | fixed in `assets/viz_base.css` (L1: `.viz-view-title b:not([class*="c-"])`), re-inlined everywhere; E6's local workaround is now redundant |
| **Explain text assumed the sign of the other roots** ("both decay" while the conjugate root grows equally, or a second real root is positive) | ch11 E4 (round-1 Must) | derive every such phrase from the computed signs: "its conjugate grows equally; the third decays" / "one more grows, one decays" / "both decay" |
| **A fitted slope that depended on an unrelated control** (Lyapunov fit 0.89 at δ₀ = 1e-8 but 0.59 at 1e-4 because the fitting window shrank) | ch11 E9 | print "no fitted line: window too short (2.7 decades)" below two decades; explain the flat start of the separation curve (≈ 12 time units) |
| **Slider value vs readout disagreed** (Ta slider 1708, readout −175 past Rayleigh's line) | ch11 E5 | show the reachable value on the slider, mark the requested one as "not reachable here" on the plane |
| **Code-tab lines wrap at 1366×768 above ≈ 44 characters** (E5 18/20, E9 16/23, E2 14/17 lines; comments detach from their lines) | ch11 E1, E2, E4, E5, E6, E9 | reflow to ≤ 44 characters (the skill's "≤ 70" is too generous for the side panel); lint candidate L5 |
| **`#<param>=…` deep links are overridden by tour step 1** (boot applies the hash, then `goStep(0)` applies the step's `set`) | all nine ch11 explainers (library behaviour) | L4: apply hash parameters after `goStep(0)` |
| **`shot.py` crashed in the parity stage while a solver module was being edited** (`AttributeError: … no attribute 'rayleigh_eigs_contour'`) | ch11 viz review run 1 | never run the implementer and a shot audit at the same time; a crashed run is not a result |
| **Builders sharing one scratchpad overwrote each other's generator scripts** (second chapter in a row) | ch11 viz phase | one private scratch subfolder per agent; generators that produce embedded tables are saved in `tools/viz_tables/chNN/` |
| **A mirrored fluidpy text containing "(11.46)" tripped `eq_refs.py`**, so the builder split the string (`const EQTAG = ' of (' + '11.46)'`) to hide it | ch11 E4 (open) | L6: an allow-marker in `tools/eq_refs.py` for strings mirrored byte for byte from fluidpy; do not split strings |
| **Thin strip on phones**: seven Bénard cells in a ~45 px view; a rotated y title losing its last character at 360×640 | ch11 E1, E8 (open) | draw fewer wavelengths when the view is short; shorten axis titles ("k L (log)") |
| **KaTeX double superscript**: `Viz.tnum(x) + '^2'` (or `^{2/3}`) is invalid TeX when x prints as `a\times10^{b}` — a boxed Explain line shown as red source text at the default state (E1) and on a preset and at slider ends (E4); `shot.py` passed both | ch12 E1, E4 (round-1 Must) | bracket every formatted number before a power (local `tsq`/`tpow`: `\big(3.27\times10^{-5}\big)^{2}`); run a `.katex-error` probe over presets × slider ends × tabs × pager pages; library `Viz.tpow` and a shot.py check are the top ch12 candidates |
| **A typed number went stale** ("C₆ = 2.240" where the page computes 2.190 after a default changed) | ch12 E6 (round-1 Must) | compute every displayed constant from the page's own functions and show the arithmetic; scan reader strings for decimal literals before review |
| **A placeholder left in a derivation step** (`\overline{u_i\times(\text{step 3})}` instead of the eight-term equation of D10 step 4) | ch12 E5 (round-1 Must) | a derivation step always shows its own equation; long lines go in `aligned` rows (three rows; one row per page on a phone) |
| **A Code block that did not run** (`prod` used undefined `ens_u`, `ens_v`, `k`) | ch12 E1 (round-1 Must) | run every Code block as displayed in several states before review |
| **Raw library strings as badges or readouts**: fluidpy's verdict strings start with the criterion ("stable ⇔ dT/dz > Γa: −12.0 < −9.8 K/km" for an unstable layer) and write "Γa" for +9.8 K/km in the meteorological line; a regime name ("forced convection") read as a stability verdict | ch12 E9 (design; avoided in the build, still visible in the Code tab) | show the bare verdict word beside the string, translate the other convention's symbol (Γ_d), and never derive stability from a regime label |
| **"The page prints" used for the book** and slips not in the house wording | ch12 E10 (round-1 Must) | "⚠️ slip #k — the book prints …; the correct form is …" with both equations written out |
| **Three stacked views on a 360×640 phone**: Explore 10 pages (E7, E8), a walkthrough step over 4 card pages (E9 step 1), D10 step 4 over 9 pages (E5) | ch12 E5, E7, E8, E9 (open, Should) | show two views at a time and swap by step (E2 did); cut Explore readouts to five; drop `notes: true` from long steps |
| **`hidePortrait` leaves a hole**: a hidden view keeps its flex share of the row, so its neighbour does not grow | ch12 E8, E9 (local workaround: hide the whole row by scoped CSS and set `flex-grow` by hand) | library: per-layout `hideOn` that removes the view from the flex row |
| **`<sub>` below 12 px** in view titles and status text | ch12 E2, E7, E8 (local `sub { font-size: max(12px, .8em) }`) | put the rule in `assets/viz_base.css` |
| **A derivation line cut over two pages was re-centred in its slice** (`.viz-der-chain` is `align-items: center`) | ch12 E8, E9 (local `.viz-der-chain.viz-sliced { align-items: flex-start; }`) | same rule into `assets/viz_base.css` |
| **Values below 1e-12 printed as 0** (a spectrum far beyond the cut-off read "= 0") | ch12 E4 | pass `keepTiny` wherever a physically tiny number is displayed |
| **View titles truncated with "…"** instead of wrapping (E5 at 1000×700; E9 at 1000×700 and 844×390) | ch12 E5, E9 (open) | shorter titles per layout; library: wrap or a `titleShort` |
| **The "All steps →" button of a quoted derivation step cut by a page slice on phones** | ch12 E2 (open) | shorten the quote card or drop the button from the quote on phones |
| **Two estimators or two models, two numbers for "the same" quantity**: DNS intercept 4.29 (mean at fixed κ) vs 4.28 (least squares); equal stresses at y⁺ = 9.9 (Spalding) vs 10.46 (mixing length); model fit 0.382 / 3.94 (page) vs 0.383 / 3.99 (notebook grid); Ri 0.12 vs 0.0752 from two formulas on one screen | ch12 E5, E7, E8, E9 (open, Should) | name the estimator or the model beside each number, or use one across the chapter |
| **Colour code reused with another meaning** (rose and blue for parts of a product; blue for uv > 0) while the chapter fixes blue = buoyancy, rose = dissipation | ch12 E2, E3 (open) | keep the chapter colour code or say the exception in Explain section 0; muted grey for "the other sign" |
| **Builders' probe scripts collided by file name in the shared scratch folder** (third chapter) | ch12 viz phase | per-slug scratch subfolders and slug-prefixed file names |

## Promotion candidates (helpers duplicated across explainers)
| Helper | Found in | Proposed library name | Status |
|---|---|---|---|
| TeX number with mantissa renormalised (10×10⁻⁵ → 1×10⁻⁴) | ch01 E1 `T`, E2 `tn`, E5 `tn` | `Viz.tnum` fix | **promoted** (ch01 knowledge pass); + `keepTiny` option for E1 |
| canvas font string in the page font | ch01 E1 `font`, E5 `font` (inline in E2, E3, E4) | `Viz.font(px, weight, mono)` | **promoted** |
| rounded-rectangle path | ch01 E5 `rrect`, E2 inline `roundRect` fallback | `Viz.roundRect(ctx, x, y, w, h, r)` | **promoted** |
| pixel text with halo / background box | ch01 E3 `canvasText`, E4 `drawText` | `Viz.text(ctx, str, x, y, opt)` | **promoted** |
| end-of-run summary card | ch01 E1 (inline rect card), E2 `card`; E3 needs it (clipped card) | `Viz.card(ctx, cx, cy, lines, opt)` | **promoted** (wraps, stays in canvas) |
| duration in everyday units | ch01 E2 `humanTime`, E4 `fmtT` (s/min) | `Viz.fmtTime(seconds)` | **promoted** |
| exact rationals | ch01 E5 `Frac`/`Q` | `Viz.Frac` | candidate — promote when Ch. 4 §4.11 (similarity) needs it in a second explainer |
| Poisson / normal / lgamma samplers matching numpy | ch01 E1 | `Viz.num.poisson`, `Viz.num.normal`, `Viz.num.lgamma` | candidate — one explainer so far (Ch. 12 turbulence statistics may reuse) |
| plot with a y-title gutter on narrow views | ch01 E4 `vplot` | `Plot` option `{ylabelGutter: true}` (or default below 420 px) | candidate — really a library defect; fix in the engine when a second explainer hits it |
| x ticks with decimals from the tick step | ch01 E4 `xTicks` | `P.axes({xdecimals: 'auto'})` | candidate — engine defect on narrow ranges |
| length in everyday units (pm … km) | ch01 E1 `lenLabel` | `Viz.fmtLength(m)` | candidate |
| number keeping trailing zeros (1.00) | ch01 E2 `fm`, `tn` | `keepZeros` option for `Viz.fmt`/`Viz.tnum` | candidate |
| diagonal hatch pattern for negative areas | ch01 E3 `hatchPattern` | `Viz.hatch(ctx, color)` | candidate |
| regime colour/word maps (stable teal / neutral amber / unstable rose) | ch01 E4 | keep local; record the colour meaning in the chapter | not a library helper |

The ch01 explainers still carry their local copies; the next rebuild of any of them should switch to the library
helpers (and E3's end card to `Viz.card`, which also closes its open Should-fix).

### ch02 candidates (listed, **not promoted**: the ch02 knowledge pass ran while the site-publisher was reading `viz/`)
Promote in a pass where nothing else reads `viz/`: edit `assets/viz_lib.js` (documented, backwards compatible) →
`tools/viz_inline.py --all` → `tools/shot.py --chapter ch01 --quick`, `--chapter ch02 --quick` and
`tools/shot.py templates/viz_example.html --quick` must PASS, otherwise revert.

| Helper | Found in (file · local name) | Proposed library name | Priority |
|---|---|---|---|
| fixed-decimal formatters, HTML minus vs TeX minus, tiny → 0 | ch02 E2 `fx`/`tx`, E3 `f2z`/`tz`, E4 `fz`/`tz`, E5 `fz`/`tz` (4 copies) | `Viz.fixed(v, d, {tex})` | **high** (4 explainers) |
| 2×2 symmetric principal values/angle | ch02 E2 and E3 `principal2d` | `Viz.num.principal2d(S)` | high (2 copies) |
| matrix cells with row/column highlight and fitted assembly lines | ch02 E1 `drawMatrix`, E3 `matrixLayout`/`drawMatrix`; ch01 E5 matrix view | `Viz.matrix(ctx, rect, M, {hot, colors, fmt})` | high (3 explainers) |
| log-x axis with manual ticks | ch02 E4 `drawLimit`, E5 `drawBars` | `P.axes({xlog: true})` | medium |
| hatched (negative) bars | ch02 E5; ch01 E3 `hatchPattern` | `Viz.hatch(ctx, color)` | medium (2 chapters) |
| tiling helpers (rectangle or disc into k×k tiles, per-tile sums) | ch02 E4 `tiled`, E5 `tilesOf`/`tileSums` | `Viz.num.tiles(shape, k)` | medium |
| diverging colormap with white zero | ch02 E4 `Viz.CMAPS.bwo` (local extension), E5 `curlImage` | `Viz.CMAPS.bwo` + `Viz.field.heatmap({diverging: true})` (vmin = −vmax) | medium |
| text clamped inside the plot rectangle | ch02 E4 `labIn` | `Viz.textIn(P, str, x, y, opt)` | low (1 copy) |
| angle arc with a label | ch02 E1 `drawArc` | `P.arc(a0, a1, r, {color, label, labelRadius})` | low |
| label pushed past an arrow tip | ch02 E1 `tipLabel`, E2 inline `tipLabel` | `P.tipLabel(from, to, str, {color, offset})` (needs collision-avoiding offset) | low |
| try progressively shorter strings until one fits | ch02 E1 (inside `drawMatrix`) | `Viz.fitText(ctx, alternatives, maxW)` | low |
| filled polygon in plot coordinates | ch02 E2 `fillPoly` (+ `halfSquare`) | `P.fillPoly(pts, color, alpha)` | low |
| closed-form 2×2 matrix exponential; ε array | ch02 E3 `expm2`, `eps3` | `Viz.num.expm2(M, t)`, `Viz.num.eps3` | low now; **Ch. 3 pathlines will want it** |
| `Viz.card` with `align: 'left'` | ch02 E3 | option on `Viz.card` | low |
| adaptive status (short / mid / full by width) | ch02 E2 | engine: `status: s => ({short, mid, text, tone})` | engine request |
| `terms.unit` as a function of state | ch02 E2 (stress Pa vs strain 1/s) | engine: `terms.unit: s => …` | engine request |
| readout label wrap CSS | ch02 E2 | engine CSS | engine request |
| `.viz-seg` wrap when there are many chips | ch02 E5 | engine CSS default | engine request |
| hide a view row on landscape phones (and on short portrait phones), Derivation tab exempt | ch02 E1 (request), E3/E4/E5 scoped CSS | engine: `views[i].hideLandscape` / `rows` flag | engine request |
| `Viz.num.poisson/normal/lgamma`, `Viz.Frac`, y-title gutter, `xdecimals: 'auto'`, `fmtLength`, `keepZeros` | ch01 (see table above) | unchanged | carried over |

Python promotion candidate (same rule, `fluidpy/core/`, then `pytest -q`): `fluidpy.core.anim` could save animations
inside `matplotlib.rc_context({"figure.constrained_layout.use": False})`. matplotlib 3.11's `print_figure`
re-installs constrained layout after the first saved frame, so builders currently toggle the rcParam by hand around
`subplots_adjust` animations (`notebooks/build_ch02.py` l. 3482/3507). The ch02 kinematics helpers
(`velocity_gradient_preset`, `linear_flow_map`, `deform_square`, `material_line_angle`) **moved to
`core.kinematics` in ch03** (done by the implementer; ch02 re-exports them). New Python candidate from ch03:
`show_eqs(text, EQ)` from `notebooks/build_ch03.py` → `tools/nbkit.py` (every builder needs it under the equation rule).

### ch03 candidates (listed, **not promoted**: the ch03 knowledge pass ran while the site-publisher was reading `viz/`)
Counts are explainer files that define the helper locally (`function name(` or `const name =` after the inlined
library), over ch01–ch03 (17 explainers). Same procedure as above; add `--chapter ch03 --quick` to the shot runs.
Ranked by payoff:

| Rank | Helper | Found in (file · local name) | Count | Proposed library name |
|---|---|---|---|---|
| 1 | fixed-decimal formatters (HTML minus / TeX minus, round-off → 0) | `fx` in ch01 `buckingham_pi_machine`, ch02 `cauchy_traction_principal_axes`, ch03 E1 E2 E3 E4 E5 E7; `tx` in ch02 cauchy, ch03 E1 E2 E3 E4 E5 E7; variants `fz`/`tz` in ch02 gauss, stokes, `f2z`/`tz` in ch02 strain; ch03 E6 `F`/`T` | **11 explainers** (fx 8, tx 7) | `Viz.fixed(v, d, {tex})` (+ `keepZeros`) |
| 2 | per-layout view-row hide (portrait / landscape / short portrait), Derivation tab exempt; preset strip and mode chips dropped on short screens | scoped `<style data-chapter>` in **all 7 ch03** explainers + ch02 E3, E4, E5 | **10 explainers** | engine: `views[i].hideOn: ['portrait', 'landscape', 'short']`, `rows[i].hideOn`, `presets.hideOn`, `modes.hideOnTabs` (keep the `:not([data-tab="derive"])` rule built in) |
| 3 | waterfall term bars with signed segments (`drawBars`) and hatched negatives (`hatchRect`, `hatch`, `hatchPattern`) | `drawBars` ch02 gauss, stokes; ch03 E2, E3, E7; hatch: ch01 `heat_work_paths` `hatchPattern`, ch03 E3 `hatch`, E7 `hatchRect` | **5 (bars) + 3 (hatch)** | `Viz.bars(v, items, {waterfall, measured, hatchNegative, arrowNegative})`, `Viz.hatch(ctx, rect, color)` |
| 4 | 2×2 symmetric principal values/angle | `principal2d` ch02 cauchy, strain; ch03 E4; `principal` ch03 E5 | **4** | `Viz.num.principal2d(S)` |
| 5 | closed-form 2×2 matrix exponential | `expm2` ch02 strain; ch03 E4, E5 | **3** | `Viz.num.expm2(M, t)` |
| 6 | 2×2 SVD (finite-time ellipse axes of e^{Gt}) | `svd2` ch03 E5 | 1 (pairs with `expm2`; Ch. 13 frontogenesis will want it) | `Viz.num.svd2(M)` |
| 7 | JS port of `show_eqs` (write the equation next to its first bare number) | `showEqs` ch03 E5 (Python twin `show_eqs` in `notebooks/build_ch03.py`) | 1 JS + 1 Python | `Viz.showEqs(text, EQ)` in the library and `nb.show_eqs(text, EQ)` in `tools/nbkit.py` (the rule is project-wide) |
| 8 | ticks drawn inside narrow plots | `inTicks` ch03 E4, E7 | **2** | `P.axes({ticksInside: 'auto'})` below ~90 px |
| 9 | inline legend that returns false when it does not fit (caller moves it into the title) | `legend` ch03 E7; E3 needs it (open Should-fix) | 1 (+1 needed) | `Viz.legend(v, P, items)` → fits |
| 10 | Gauss–Legendre nodes mirroring `core.integral_theorems.gauss_legendre_nodes` | `gl`/`glStd` ch03 E7 (+ `midpoint`) | 1 | `Viz.num.gaussLegendre(n, a, b)`, `Viz.num.midpoint(a, b, n)` (Ch. 4 CV budgets will reuse) |
| 11 | viewport width class (phone / mid / wide) for text variants | `widthClass` ch03 E4, E7 | **2** | `Viz.widthClass()` (or expose the engine's layout/density) |
| 12 | scientific-notation formatters with Unicode exponents | `fsci`/`tsci` ch03 E4, `dec` + `fs3`/`ts3` ch03 E7 | **2** | `Viz.fmt(v, {sci: true})` / `Viz.tnum` exponent option |
| 13 | filled polygon in plot coordinates | `fillPoly` ch02 cauchy; ch03 E4 (`polyFill` E7) | **3** | `P.fillPoly(pts, color, alpha)` |
| 14 | label past an arrow tip | `tipLabel` ch02 rotation, cauchy | 2 (carried from ch02) | `P.tipLabel(from, to, str, opt)` |

### ch04 candidates (listed, **not promoted**: the ch04 knowledge pass ran while the site-publisher was reading `viz/`)
Counts = explainer files (of all 26, ch01–ch04) that define the helper locally after the inlined library, counted by
script (`function name(` / `const name = (…) =>`); "(ch04 n)" = how many of them are ch04's. Same procedure as above,
plus `tools/shot.py --chapter ch04 --quick`. Ranked by payoff:

| Rank | Helper | Found in (file · local name) | Count | Proposed library name |
|---|---|---|---|---|
| 1 | **number formatters** (fixed decimals, sig figs, HTML vs TeX minus, round-off → 0, unit switch) | every explainer: `fx`/`tx` (9/8), `tn` (10), `fm` (7), `fz`/`tz`, `f2z`, `f4`/`t4`, `fp`/`tp`, `F`/`T`, `fsci`/`tsci`, `fs3`/`ts3`, `dec`; ch04 E7 also `fT` (K/mK/µK), 3 ch04 files `fT`, `fL` | **26 of 26** (ch04 9) | `Viz.fixed(v, d, {tex, keepZeros})`, `Viz.fmt(v, {sci, unitScale: 'K'\|'m'\|'s'})` — the single most duplicated code in the repo |
| 2 | waterfall / term bars with signed segments, hatched negatives, measured ◇ | `drawBars` ch02 gauss, stokes; ch03 E2, E3, E7; **ch04 E1, E4, E5, E8**; `hatchRect` ch03 E7, ch04 E1 | **9** (ch04 4) + hatch 3 | `Viz.bars(v, items, {waterfall, measured, hatchNegative, signedLabels})` |
| 3 | per-mode control/readout hiding (`.xxx-off` toggled in `applyMode`) | ch03 E7 `.rtt-off`; ch04 E1 `.cvb-off`, E2 `.sfs-off`, E3 `.nsl-off`, E5 `.rfc-off`, E6 `.wb-off`, E9 `.dsm-off` | **7** (ch04 6) | engine: `params[k].modes: ['wake', 'bore']`, `readouts[i].modes`, and `modes.onChange` built in |
| 4 | scoped per-layout row hide (portrait / landscape / short), Derivation tab exempt | `<style data-chapter>` in all 9 ch04 files (+ 10 earlier) | **19** | engine `views[i].hideOn`, `rows[i].hideOn` (carried from ch03 rank 2) |
| 5 | image of a scalar field inside a plot (off-screen canvas, colour ramp, NaN transparent) | ch02 `heatFor`, `curlImage`; ch03 E6 `shadeImage`; ch04 E8 `heatImage`, E2 `bandImage` | **5** (ch04 2) | `Viz.field.image(P, f, {vmin, vmax, cmap, diverging, bands})` (library `heatmap` works on grids, not callables) |
| 6 | 2 × 2 symmetric principal values/angle | `principal2d` ch02 ×2, ch03 E4; `principal` ch03 E5; `princ2` ch04 E3 | **5** (ch04 1) | `Viz.num.principal2d(S)` (carried from ch03 rank 4) |
| 7 | closed-form 2 × 2 matrix exponential | `expm2` ch02 strain, ch03 E4, E5, **ch04 E3** | **4** | `Viz.num.expm2(M, t)` |
| 8 | label past an arrow tip | `tipLabel` ch02 rotation, cauchy, **ch04 E5** | **3** | `P.tipLabel(from, to, str, {color, offset})` |
| 9 | Gauss–Legendre nodes mirroring `gauss_legendre_nodes` | `gl`/`glStd` ch03 E7, **ch04 E2, E9** | **3** | `Viz.num.gaussLegendre(n, a, b)` |
| 10 | filled polygon in plot coordinates | `fillPoly` ch02 cauchy, ch03 E4, **ch04 E3** | **3** | `P.fillPoly(pts, color, alpha)` |
| 11 | `clamp(x, a, b)` | ch02 cauchy, ch03 E2, E5, ch04 E3, E4, E7 | **6** | `Viz.num.clamp` (trivial, but six copies) |
| 12 | **double-precision erf / erfc** | `erfHP` ch04 E1, `erfcHi` ch04 E4 | **2** | fix `Viz.num.erf`/`erfc` to ~1e-15 (series below 2.5, Lentz continued fraction above) — a library defect, not only a helper |
| 13 | arrow in pixel coordinates (tail, head, width, head size) | `pxArrow` ch04 E1, `arrowPx` ch04 E4 (+ E5 `arrowScale`) | **2** (ch04 only) | `Viz.arrowPx(ctx, x0, y0, x1, y1, {color, width, head})` (`P.arrow` is plot-space only) |
| 14 | 3-vector maths (`add`, `cross`, `dot`, `norm`) | ch04 E5 (four functions), E2 `cross`; ch03 E4 `dot` | **3** | `Viz.vec3.{add, sub, scale, dot, cross, norm}` — Ch. 13 (3-D Coriolis, Ekman) will need it |
| 15 | log–log axes with decade ticks | `logAxes` ch04 E9; ch02 E4 `drawLimit`, E5 `drawBars` log-x; ch03 E6/E7 log panels | 1 named (5 hand-made) | `P.axes({xlog: true, ylog: true})` |
| 16 | streamline tracing for a callable field in plot space | `traceLine` ch04 E6 (library `Viz.field.streamline` needs a grid) | 1 | `Viz.field.traceCallable(u, x0, {h, n, both})` |
| 17 | sine-mode table for series profiles; profile curve x = f(y) drawn/filled | `sinTable`, `vline`, `hfill` ch04 E7 | 1 | `P.fnY(f, {y0, y1})`, `P.fillBetweenX` (Ch. 5, 8 profiles) |
| 18 | per-mode transport range and rate | ch04 E1 `TMAX`/`RATE` rewrites `a.params.t.max`, `a.cfg.transport.max` and the slider DOM | 1 (hack) | engine: `transport.max`/`rate` as functions of state, or `modes.options[i].transport` |

**Top five for the next library pass**: formatters (26), waterfall bars (9), per-mode control hiding (7, engine
flag), field image from a callable (5), principal2d (5) — then fix erf/erfc precision (a defect) and add `vec3`
before Ch. 13.

**Library and tool quirks reported by the ch04 builders** (none fixed yet; each has a chapter workaround):

| # | Quirk | Symptom | Workaround in use | Engine fix |
|---|---|---|---|---|
| Q1 | Code-tab `{{placeholder}}` names share one namespace across all code snippets | two snippets using the same placeholder name for different quantities show the same live value | give every snippet's placeholders distinct names | scope `live()` results per snippet id |
| Q2 | The status badge sits over the top-left corner of the first view | it covers whatever is drawn at the view's top-left | keep that corner free of labels | reserve the badge's height in the stage layout (or place it in the strip) |
| Q3 | Transport `max` and `rate` are fixed at app creation | per-mode time ranges need three writes (params, cfg, DOM slider) | E1 `applyMode` rewrites all three | `transport.max`/`rate` accept a function of state |
| Q4 | No per-mode control/readout hiding | Explore lists every mode's controls (10–15 pages on phones) | scoped `.xxx-off` class (rank 3 above) | `params[k].modes`, `readouts[i].modes` |
| Q5 | `onChange` runs after a walkthrough step's `set()` and can overwrite the values the step just set | a step that sets a mode and slider values can end with the mode's defaults | let `onChange` reset dependants only for the key that changed, and re-check the step's screenshot | apply `set()` values after `onChange`, or pass a `fromStep` flag |
| Q6 | Mode chips do not wrap | many chips overflow a phone strip | shortened chip labels, scoped `flex-wrap: wrap` (ch02 E4/E5 pattern) | `.viz-seg { flex-wrap: wrap }` by default (carried from ch02) |
| Q7 | `Viz.num.erf`/`erfc` precision ~1.2e-7 | parity rows at 1e-8 fail; second differences of erfc profiles are noise | local `erfHP`, `erfcHi` (rank 12) | double-precision implementations |
| Q8 | `tools/eq_refs.py` flags selftest row names (and `<meta name="viz:concept">`) that cite an equation number, although both are meant to be exempt | false hits during review (E6's meta is the only one left) | reword names/meta without numbers | exempt `selftest` `name:` strings and `<meta>` content in the scanner |

### ch05 candidates (listed, **not promoted**: the ch05 knowledge pass ran while the site-publisher was reading `viz/`)
Counts = explainer files that define the helper locally after the inlined library, counted by script over **ch04 + ch05
(18 files)**; "(ch05 n)" = how many are ch05's; "all" = over all 35 explainers ch01–ch05 where it matters. Same
procedure as above, plus `tools/shot.py --chapter ch05 --quick`. Ranked by payoff:

| Rank | Helper | Found in (file · local name) | Count ch04+ch05 | Proposed library name |
|---|---|---|---|---|
| 1 | **number formatters** (fixed decimals, sig figs, TeX vs HTML minus, round-off → 0) | ch05: `fz`/`tz` (E3, E5, E6, E8), `f3`/`t3`/`t4` (E2, E4, E5, E7), `f4` (E2, E4, E5), `fm`/`tn` (E1, E9), `tp` (E3, E6, E8); ch04 all nine | **18 of 18** (all 35 of 35) | `Viz.fixed(v, d, {tex, keepZeros, roundoff})` — still the most duplicated code |
| 2 | **per-mode control/readout hiding** (`.xxx-off` class toggled in `applyMode`) | ch05 all nine (`.vtc-off`, `.vpf-off`, `.kml-off`, `.bt-off`, `.vst-off`, `.bsf-off`, `.ver-off`, `.pvl-off`, `.vsr-off`); ch04 E1, E2, E3, E5, E6, E9 | **15** (ch05 9; 16 with ch03 E7) | engine: `params[k].modes`, `readouts[i].modes` (quirk Q4) |
| 3 | **arrow in pixel coordinates** | `arrowPx` ch04 E4; ch05 E2, E4, E5, E6, E7, E8; variants `pxArrow` ch04 E1, `arrowScale` ch04 E5 | **7** (+2 variants = 9; ch05 6) | `Viz.arrowPx(ctx, x0, y0, x1, y1, {color, width, head, minLen})` (`P.arrow` is plot-space only; a `minLen` closes the ch04 E5 "half arrow" Should-fix) |
| 4 | **3-vector maths** (`add`, `sub`, `dot`, `cross`, `norm`, `scale`, `vec`) | ch04 E5; ch05 E2 (`norm`), E5, E6, E8 (`add`) | **5** (ch05 4) | `Viz.vec3.{add, sub, scale, dot, cross, norm, unit}` — Ch. 13 (3-D Coriolis, Ekman) and Ch. 14 filaments need it |
| 5 | **log axes with decade ticks** | `logAxes` ch04 E9, ch05 E4, E6; `logTicks` ch05 E3, E6, E9 | **5** (ch05 4) | `P.axes({xlog: true, ylog: true})` |
| 6 | waterfall / term bars | `drawBars` ch04 E1, E4, E5, E8; ch05 E5, E7 | **6** (11 overall) | `Viz.bars(v, items, {waterfall, measured, collapseZeros, minRowPx: 25})` — add the E7 zero-row collapse |
| 7 | image of a scalar field from a callable | `washImage` ch05 E2, E4; `heatImage` ch04 E8, `bandImage` ch04 E2 | **4** (7 overall) | `Viz.field.image(P, f, {vmin, vmax, cmap, diverging, bands})` |
| 8 | Gauss–Legendre nodes mirroring `gauss_legendre_nodes` | `glNodes` ch05 E3, E6; `glStd` ch04 E2, `gl` ch04 E9 | **4** (5 with ch03 E7) | `Viz.num.gaussLegendre(n, a, b)` |
| 9 | vector RK4 step / integrator for many points | `rk4` ch05 E1, E3, E8 (library `rk4Step` works on one state; these carry arrays of points) | **3** | `Viz.num.rk4Many(f, X, t, dt)` or document `rk4Step` on flat arrays |
| 10 | light 3-D canvas projector (orbit angles → 2-D, depth sort) without three.js | `projector` ch05 E3, E5; `proj` ch05 E6 | **3** | `Viz.proj3(view, {theta, phi, scale})` with `.p(x)`, drag-to-orbit — cheaper than `Viz.three` for line drawings |
| 11 | `app.set` wrapper that re-applies a step's values after `onChange` (quirk Q5 workaround) | ch05 E5, E8 | **2** | engine fix Q5: apply `set()` after `onChange` (or pass `fromStep`) |
| 12 | lazy trajectory cache keyed by the parameters | `cache` ch05 E3, E6 | **2** | `Viz.memo(fn, keyOf)` |
| 13 | length in everyday units | `fmtLen` ch05 E2, E4; ch01 E1 `lenLabel` | **2** (3 overall) | `Viz.fmtLength(m)` |
| 14 | polygon in pixel coordinates | `polyPx` ch05 E6, E8; `fillPoly` (plot space) ch04 E3 | **3** | `P.fillPoly` + `Viz.polyPx(ctx, pts, opt)` |
| 15 | `clamp(x, a, b)` | ch04 E3, E4, E7; ch05 E9 | **4** (7 overall) | `Viz.num.clamp` |
| 16 | complete elliptic integrals K(m), E(m) by AGM with the parameter m | `ellipKE` ch05 E6 | 1 | `Viz.num.ellipKE(m)` — Ch. 14 rings/filaments will reuse; document m = k² |
| 17 | round-off-aware term value ("≈ 0 (round-off)") | `termRO` ch05 E3 | 1 (E3's Helmholtz title needs it too) | option `roundoff: true` on `Viz.tnum` / term items |
| 18 | general 3 × 3 matrix exponential (scaling and squaring) | `expm` ch05 E5; `expm2` ch04 E3 (4 overall) | **2** | `Viz.num.expm(M, t)` (2 × 2 closed form as a fast path) |
| 19 | double-precision erf/erfc | `erfHP` ch04 E1, `erfcHi` ch04 E4; none in ch05 | 2 | fix `Viz.num.erf/erfc` (quirk Q7) |
| 20 | label past an arrow tip | `tipLabel` ch04 E5 (ch02 ×2); none in ch05 (ch05 labels arrows inline) | 1 (3 overall) | `P.tipLabel(from, to, str, opt)` |

**Top five for the next library pass** (unchanged in spirit, re-ranked with ch05): formatters (35/35), per-mode hiding
(16, engine flag), `arrowPx` (9 arrow helpers), `vec3` (5) and log axes (5); then the narrow-view y-title strip as an
engine default (three phone Must-fixes in ch01 and ch05 came from it), `Viz.bars` with zero-row collapse, field image,
Gauss–Legendre, and the Q5 fix that removes the `app.set` wrappers. No ch05 builder reported a new engine quirk beyond
Q1–Q8; Q4 (per-mode controls) was hit in all nine ch05 explainers and Q5 (`onChange` after `set`) in two.

### ch07 candidates (listed, **not promoted**: the ch07 knowledge pass ran while the site-publisher was reading `viz/`)
Counts = explainer files (all 53, ch01–ch07) that define the helper locally outside the inlined library, counted by
script (`function name(` or `const name = … =>`, plus chapter CSS rules); "(ch07 n)" = how many of the nine ch07 files.
Ranked by payoff for the chapters ahead (Ch. 8 laminar, Ch. 10 CFD, Ch. 11 instability, **Ch. 13 GFD**, Ch. 15):

| Rank | Helper | Found in (file · local name) | Count (ch07) | Proposed library name | Recommendation |
|---|---|---|---|---|---|
| 1 | **pager label on one line** | chapter CSS `.viz-pager-label { white-space: nowrap; min-width: 52px; }` in `hydraulic_jump`, `internal_wave_beams`, `particle_orbits`, `two_layer_modes`, `wave_rays_refraction` | 5 (5) | one rule in `assets/viz_base.css` | **promote first** — CSS only, backwards compatible, removes five local copies; re-run `--quick` shots of all chapters (the header height changes by at most one line) |
| 2 | **pixel-space arrow and dot** | `arrowPx` ch04–ch07 (ch07: capillary, packets, internal, orbits, two_layer, rays); `dotPx` ch06 E6, ch07 capillary, packets, internal, two_layer; `pxArrow` ch04 E1, ch07 `hydraulic_jump`, `seiche_standing_waves`; `linePx` packets, internal | 21 + 5 + 3 + 2 (6 + 4 + 2 + 2) | `Viz.arrowPx(ctx, x0, y0, x1, y1, {color, width, head, minLen})`, `Viz.dotPx(ctx, x, y, r, fill, stroke)` | **promote** — the top candidate by count since ch05; three local names for one arrow |
| 3 | **narrow-view plots: y-title gutter and square data windows** | `gutterPlot` ch05–ch07 (ch07: capillary, packets, internal, orbits, rays); `squarePlot(v, R, small)` `internal_wave_beams` (square ±R window centred in the view with its own y-title strip) | 12 + 1 (5 + 1) | engine default below ~420 px (or `v.plot({ylabelGutter: true})`) and `v.plot({square: true})` | **promote** — the phone-legibility Must-fix of ch01, ch05, ch06 and ch07 E9; the square K plane is Ch. 13's wavenumber-space view |
| 4 | **number, length and speed formatters** | `f3` 28 files (ch07 9/9), `t3` 21 (9/9), `t4` 21 (7), `f4` 14 (1); `fmtLen`/`texLen` (`dispersion_relation`, `seiche_standing_waves`, ch05 ×2), `lenStr`/`lenTex` (capillary, packets, two_layer); `spd`/`spdTex` (packets, internal, two_layer: m/s ↔ cm/s ↔ mm/s) | 28 / 21 / 5 / 3 (9 / 9 / 5 / 3) | `Viz.fixed(v, d, {tex, keepZeros, roundoff})`, `Viz.fmtLen(m, {tex})` (µm … km), `Viz.fmtSpeed(v, {sig})` | **promote** — unchanged top candidate since ch04; ch07 adds length and speed units (4 s.f. option for rule-9 numbers) |
| 5 | **wave numerics mirroring `core.waves`** | `tanhKH`, `xOverSinh` (capillary, dispersion, packets, rays); `coshSinh`/`sinhSinh` (dispersion, seiche); `coshOverSinh`/`sinhOverSinh` (orbits, overflow-safe); `omegaG` (orbits, seiche, rays); `omegaCG`, `cgv` (capillary, packets); `cph` (capillary, dispersion, packets); `kFromT` (dispersion, Newton inverse of (7.28)); `kFromOmega` (rays) | ≈ 8 names × 2–4 files (ch07 only) | `Viz.waves = {tanhKH, xOverSinh, coshOverSinh, sinhOverSinh, coshOverCosh, omega(k, H, g, sigma, rho), c, cg (signed), kFromOmega(omega, H, …)}` with selftest parity against `core.waves` | **promote before Ch. 13** (Poincaré, Kelvin, Rossby explainers will need ω(k), c, c_g and the inverse again); keep the overflow-safe forms exactly as in `core.waves` |
| 6 | **per-mode control and readout hiding** | `markMode` (packets, internal, two_layer: data attribute + CSS); `applyMode` 13 earlier files; `applyVisibility` 3 | 3 + 13 + 3 (3) | engine `params[k].modes` / `readouts[k].modes` (quirk Q4) | **engine fix** — 19 local workarounds; also the cure for 10–19 page Explore/Explain on phones |
| 7 | **layout test** | `narrow()` in 9/9 ch07 files (window width < 560 or height < 460/500) | 9 (9) | read `g.app.layout` / density from the engine instead | **replace, don't promote** — a window test disagrees with the engine's layout inside a notebook iframe |
| 8 | **log axes with decade labels** | `logAxes` ch04 E9, ch05 ×2, ch07 `dispersion_relation`, `two_layer_modes`; `logTicks` 4 earlier | 5 + 4 (2) | `P.axes({xlog, ylog})` | promote with rank 3 (engine plot work) |
| 9 | complex arithmetic | `csqrt` ch06 E5, E6, ch07 `group_velocity_packets`; ch06 `cmul`/`cdiv`/`cabs` | 3 (1) | part of `Viz.cx` (ch06 rank 1) | promote with `Viz.cx` (Ch. 11 normal modes, Ch. 14) |
| 10 | Gauss–Legendre quadrature | `glInt` ch07 `wave_rays_refraction`; `glNodes` ch05 ×2, ch06 E9 | 1 + 3 (1) | `Viz.num.gaussLegendre(n, a, b)` | promote (four files now) |
| 11 | batch particle integration with dense lookup | `integrateMany`, `posAt` ch07 `particle_orbits`; `posAt` ch05 `point_vortex_lab` | 2 (1) | `Viz.num.odeintMany(f, X0, t)` + interpolated `posAt` | wait for a third user (Ch. 13 floats) |
| 12 | wave-specific drawing | `crestLines` (rays), `nearestReal` (dispersion, jump, seiche: nearest row of a real-world table by log distance) | 1 / 3 | `nearestReal` → a `current: 'nearest'` option of `Viz.work.table`; `crestLines` stays local | small; with rank 4 |

**Recommendation for the next library pass (orchestrator, when nothing reads `viz/`)**: ranks 1–4 first (pure additions:
CSS rule, `Viz.arrowPx`/`dotPx`, gutter + square plots, formatters incl. `fmtLen`/`fmtSpeed`), then `Viz.waves` (rank 5)
and per-mode hiding (rank 6) **before Ch. 13's explainers**, then log axes, `Viz.cx`, Gauss–Legendre; replace `narrow()`
by the engine's layout; keep `crestLines`, `integrateMany` local. After promotion: `tools/viz_inline.py --all`, then
`tools/shot.py --chapter ch01 … ch07 --quick` and `tools/shot.py templates/viz_example.html --quick` must all PASS.

### ch12 candidates (listed, **not promoted**: the ch12 knowledge pass ran while the site-publisher was reading `viz/` and `assets/`)
From `reports/ch12_viz.md` ruling l (14 library and tool findings) and a grep for named definitions in the ten chapter
files (E1 `reynolds_averaging_window`, E2 `correlation_and_spectrum`, E3 `reynolds_stress_parcels`, E4
`energy_cascade_spectrum`, E5 `turbulent_energy_budget`, E6 `turbulent_jet_similarity`, E7 `law_of_the_wall`, E8
`mixing_length_closure`, E9 `stratified_surface_layer`, E10 `taylor_dispersion`). Where the reviewer's list is wider than
the grep, the helper is an inline expression there.

| Rank | Helper | Files carrying a local copy (ch12) | Proposed library name | Why |
|---|---|---|---|---|
| 1 | **log axes for `v.plot()`** (`logFrame`, `L10`, log ticks, semi-log and log–log frames) | **eight files**: `logFrame` in E1, E4, E6, E10; `L10` in E1, E4, E5, E7, E8, E9, E10; `logTicks`-type helpers in E7, E8 (+ the ch10 list: E1, E2, E4) | `v.plot({xlog: true, ylog: true})` with decade ticks and minor ticks | every spectrum, every wall profile, every dispersion curve; Ch. 13 needs it on day one (spectra, Ekman depth vs K, Rossby dispersion on log k); third chapter asking |
| 2 | **bracketed powers of formatted numbers** (`tsq`, `tpow`: `\big(a\times10^{b}\big)^{n}`) | `tsq` in E1, E4, E5, E8, E9; `tpow` in E1, E10 | `Viz.tpow(x, n, sig)` (and `Viz.tsq`) | the round-1 KaTeX double-superscript failures (E1, E4) were this; one tested helper removes the class |
| 3 | **`fitText`** (shrink-to-fit a canvas label, never below 12 px, else wrap or drop) | E5, E8, E9 | `Viz.fitText(ctx, str, maxWidth, {min: 12})` or a `maxWidth` mode of `Viz.text` | labels on bars and narrow views; the truncated-title finding (11) is the same need |
| 4 | **`erf` to 1e-14** (series + continued fraction) | E6, E10 (both also call `Viz.num.erf`, good to ≈ 1e-7 only — their rows sit at 1e-6) | fix `Viz.num.erf`/`erfc` in place (asked since ch08) | parity rows at 1e-10 instead of 1e-6; Gaussian profiles and dispersion in Ch. 13 |
| 5 | **Gauss–Legendre quadrature `gl8`** (8-point rule on a panel, composite) | E5, E8 (local `quad`-type helpers also in E2, E3, E6) | `Viz.num.gauss(f, a, b, {n: 8, panels})` beside `trapz` | wall-profile integrals, spectra and budget integrals to 1e-10 with few evaluations |
| 6 | **seeded Gaussian and Ornstein–Uhlenbeck draws** (Box–Muller on `Viz.rng`; exact OU update e^{−Δt/τ_c}) | E1, E6, E10 (ch01 E1 has `normal`/`poisson` samplers) | `Viz.rng(seed).normal()`, `Viz.num.ou(rng, n, dt, sigma, tau_c)` | every stochastic explainer; must match the fluidpy recursion for statistical parity rows |
| 7 | **log-paced transport** (a 0…1 clock mapped to time on a log scale by hand) | E10 | `transport: { param: 't', min, max, scale: 'log' }` | dispersion, spin-up and adjustment span decades (Ch. 13) |
| 8 | **two-line status** (E9 injects `<br><span class="ssl-conv">` into the status text) | E9 | `status: s => ({ text, sub, tone })` | the two-convention badge; any verdict that needs its criterion underneath |
| 9 | **per-layout view hiding that frees the space** (`hidePortrait` keeps the flex share; E8 and E9 hide the whole row by scoped CSS and set `flex-grow` by hand) | E8, E9 (+ ch10 list `views[i].hideOn`) | `views[i].hideOn: ['portrait', 'landscape']` removing the view from the flex row | three views cannot stack on a 360×640 phone |
| 10 | **`.viz-der-chain.viz-sliced { align-items: flex-start; }`** (a sliced derivation line was re-centred in its slice) | E8, E9 | rule in `assets/viz_base.css` | CSS fix, no API |
| 11 | `sub { font-size: max(12px, .8em) }` in titles and status | E2, E7, E8 (grep also E10) | rule in `assets/viz_base.css` | 12 px floor |
| 12 | hide the preset strip on phones outside Explore (the same six CSS lines) | at least seven files (reviewer's count) | an app option `presets: { phones: 'explore' }` | removes a copy-pasted block |
| 13 | `arrowPx` (fifth chapter) | E3, E5, E6, E8, E9 (+ 24 earlier files) | `Viz.arrowPx(ctx, x0, y0, x1, y1, opts)` | unchanged from the ch09–ch11 lists |

Order for the library pass (when nothing reads `viz/` or `assets/`): ch12 ranks 1–4 and 10–11 first (they remove audited
failures), then 5–9; with the ch11 list (complex numbers, `bisect`/`goldenMax`, table interpolation, cubic solver) and the
ch10 list (`Viz.num.linalg`, `gutterPlot`). Then `tools/viz_inline.py --all`, `tools/shot.py --chapter ch01 … ch12 --quick`
and `tools/shot.py templates/viz_example.html --quick`; revert on any failure.

### ch12 library and tool findings (the reviewer's 14, verbatim in substance; nothing was edited)
1. `.viz-der-chain` is `align-items: center`: a line the pager cuts over two pages is re-centred in its slice (rank 10).
2. A `hidePortrait` view keeps its flex share of the row (rank 9).
3. `<sub>` in titles and status falls below 12 px (rank 11).
4. No log axes in `v.plot()`: local helpers in eight files (rank 1).
5. `Viz.num.erf` is good to about 1e-7 only (rank 4).
6. No log-paced transport (rank 7).
7. No two-line status (rank 8).
8. `Viz.tnum(x) + '^2'` is a KaTeX double superscript when x prints as a×10^{b} (rank 2).
9. **`tools/shot.py` does not fail on `.katex-error` nodes and audits Explain only in the default and step states.**
10. A quoted derivation step's "All steps →" button can be cut by a page slice on phones (E2).
11. View titles are truncated with "…" instead of wrapping (E5, E9).
12. At least seven files hide the preset strip on phones outside Explore with the same six CSS lines (rank 12).
13. `Viz.fmt`/`Viz.tnum` print values below 1e-12 as 0 unless `keepTiny` is passed.
14. Two builders' probe scripts collided by file name in the shared scratch folder.

### ch12 open explainer items (Should, non-blocking; from `reports/ch12_viz.md` round 2)
- **E3, E9, E6: label slips #16, #17, #18** in the house wording (the numbers did not exist when the files were built;
  the corrected forms are already shown). E9: move the −$\overline{uw}$ = u_*² correction from the *watch* line of D25
  step 1 into its *why*; E6: D13 says "trap".
- E1: "0.241 = 0.280 + -0.040" → "− 0.040" with a sum that closes on the printed digits; "(no wave in these runs)" on the
  N → ∞ line; **no fluidpy twin for the OU window variance** (`ch12.time_average_variance_ou` at the next implementer pass).
- E2: "All steps →" cut on `phone-tall` step 4; say the rose/blue exception in Explain section 0.
- E5: D10 step 4 takes 9 card pages at 360×640 (hide the step dots there); title truncated at 1000×700; "with damping" in
  the undamped note.
- E6: D16 step 1 *why* is 36 words.
- E7: step 2 should write $u_*^2\equiv\tau_0/\rho$ beside (12.81); Explore 10 pages on a phone (cut readouts to five);
  "equal shares at 9.9 (Spalding's curve; the mixing-length channel gives 10.5)"; "B −17.55" with a hyphen; model fit
  0.382 / 3.94 here vs 0.383 / 3.99 in the notebook.
- E8: DNS intercept "about 4.29" vs 4.28 in E7 and the notebook; Explore 10 pages on a phone.
- E9: the Code tab prints fluidpy's "Γ < Γa: 6.5 < 9.8 K/km" two lines under "dT/dz > Γa: −6.5 > −9.8 K/km" (add a comment
  line, or change the library string); ∂ beside (12.106); label the readout "Ri = Pr_T·Rf (neutral shear)"; step 1 over 4
  card pages on a phone.
- E10: a lone full stop after u_rms in the Explore intro on a phone.
- Backup B1 `k_epsilon_decay` not built (C13 has no explainer; neither has C05).

### ch12 machinery / tooling TODOs (for the orchestrator; not done in this pass)
- **`tools/shot.py`: fail on any `.katex-error` node, on every page of every pager, and walk the presets (and both ends of
  every slider) with the Explain tab open** — third chapter; the reviewer's DOM probe and string probe exist only as
  session scratch files.
- `tools/check_public.py` now reads `_forbidden_public_regex` from the private `tests/book_values_chNN.json` and exits on an
  unparsable private file; **still tracked files only** — document "run after `git add`, or pass the files explicitly".
  Seed the regex list for a chapter **before** scripts and design are written.
- `tools/nbkit.py`: move `self_check_ctrl`, `self_check_eaten` (TeX command tails without a backslash) and
  `self_check_near` out of `build_chNN.py` into the kit; generic comment phrases still generated on 113 lines (fourth
  chapter); stray `<matplotlib.legend.Legend …>` outputs (end the cell with a semicolon).
- `tools/eq_refs.py`: `ref:` labels of the form "Eq. (12.54), exponent corrected" on cards that show the TeX are flagged
  (4 hits accepted by hand) — allow a suffix after the number on `ref:` labels.
- A contract audit tool: compare `inspect.signature` of every design Part C name with the module (done by hand on resume).
- Per-slug scratch subfolders (collisions in ch10, ch11, ch12).
- `viz_lint` control-character rule (seven chapters).

### ch12 lesson candidates for the skills (not promoted in this pass: the brief forbade library and skill edits)
- `interactive-viz` Lessons: (ch12) **never append `^n` to a formatted number — bracket it (`tsq`/`tpow`)**, and run a
  `.katex-error` probe over presets × slider ends × tabs × pager pages before review · **compute displayed constants from
  the page's own functions; no typed results in reader strings** · a derivation step shows its own equation, never a
  placeholder · run every Code block as displayed · a status badge shows the bare verdict word, then the criterion in the
  chosen convention, then the other convention — never a raw library string, never a regime label as a stability verdict ·
  a fitted constant is shown with its window, and the window is draggable · term bars of a budget sum to zero and the
  remainder is labelled "remainder" · two views at a time on a 360×640 phone, swapped by step · pass `keepTiny` for
  physically tiny values · "printed vs corrected" toggles carry parity rows on both.
- `math-to-python` §7: (ch12) **an argument that carries a convention or a physical assumption is a required keyword**
  (`Gamma_a`) · **empirical constants are never silent defaults** (κ, B, Π, z₀) and illustrative ones are visibly not the
  book's · dictionary keys carry their sign (`minus_uv_plus`) · state one-sided vs two-sided beside every spectral constant ·
  "to the first zero" is not an integral scale when the correlation has a negative lobe · a fitted constant is a function of
  its window · outside its range a formula returns NaN · classifier boundaries are compared on the ratio, inclusive, and
  mirrored the same way in JS · synthetic random-phase fields are kinematic; 2-D and 3-D isotropy differ.
- `verify-implementation`: (ch12) **a test that compares a constant with its cited value is "constant consistency", not
  V5** · **a benchmark whose source nobody read first-hand is not a benchmark** (flag `verified_first_hand: false`, test the
  flag) · models against DNS are "approximate" with the band in the test · an identity that holds for any input cannot test
  which input was used (R16) · sampled assertions: seed + 5 standard errors + a measured scatter exponent · a self-audit test
  that every contract name is called by a test · run the mutant campaign on `@slow` tests too.
- `teaching-style` Lessons: (ch12) **a number may stand only beside the equation actually printed under it**; never start
  a sentence with a symbol · show an artefact beside the clean result with its predicted size, and **test an explanation by
  removing the alleged cause** · "Where this holds" under an order-of-magnitude estimate · a printed gap before a
  comparative sentence · a slip box says what is printed, why it fails, what is right, and checks both — and never
  over-claims ("a valid family", not "the valid one") · a library string in another notation is translated at every
  appearance · "theorem vs observation" for look-alike criteria · the ch12 derivation moves in `concept_map.md`.
- `colab-notebook` / nbkit: (ch12) **never write LaTeX through a shell heredoc**; the builder scans its own output for
  control characters and command tails without a backslash, and the scanner is tested on planted damage · design "expect"
  values are hypotheses (report, do not edit) · runtime is load-dependent — measure on an idle machine; trim animation frames
  before anything else · the committed notebook should be the executed one.

### ch11 candidates (listed, **not promoted**: the ch11 knowledge pass ran while the site-publisher was reading `viz/` and `assets/`)
From `reports/ch11_viz.md` ("Helpers worth promoting") and a grep for named function definitions in the nine chapter scripts
(the reviewer's file lists are wider than the grep where a helper is an inline expression or a `const`).

| Rank | Helper | Found in (ch11) | Proposed library name | Why |
|---|---|---|---|---|
| 1 | **complex-pair formatting and complex square root** (`ctex`, `ctxt`, `cpy`, `csqrtReal`: a ± bi in TeX, in text, as a Python literal for `py:` rows; √ of a real that may be negative) | reviewer: E1–E4, E9; named definitions in `normal_mode_growth.html`, `salt_fingers.html` | `Viz.cx.fmt/tex/py`, `Viz.cx.sqrtReal` (with the ch09 `Viz.cx` request: add, mul, div, exp, cos, cosh — E3 computes a complex 3 × 3 determinant) | every stability explainer prints eigenvalue pairs; Ch. 13 will too |
| 2 | **`bisect` and `goldenMax`** (bracketed root and bracketed maximum) | `normal_mode_growth.html`, `inviscid_shear_criteria.html` (`bisect`); `normal_mode_growth.html` (`goldenMax`) | `Viz.num.bisect`, `Viz.num.goldenMax` (beside `brentq`) | band edges and fastest-growing modes on every σ(k) curve |
| 3 | **interpolation of a table on a rectangular grid** (bilinear and cubic; optional snap to a known zero curve; log-spaced axes) | reviewer: E5–E8 (`cubW` in `orr_sommerfeld_neutral_curve.html`; local variants in E3, E5, E6, E7) | `Viz.num.table2d(xs, ys, values, {kind: 'bilinear'\|'cubic', logx, logy})` | the "table of fluidpy values" pattern needs one tested interpolator instead of four |
| 4 | **real-coefficient cubic solver** (three roots, real or a conjugate pair) | `salt_fingers.html` (`cubicCoeffs`, `cubicRoots`), `lorenz_attractor.html` (`cubicC`, `solveCubic`) | `Viz.num.cubicRoots(a2, a1, a0)` | characteristic cubics (double diffusion, Lorenz; Ch. 13 dispersion relations) |
| 5 | **`arrowPx`** (pixel-space arrow; asked since ch09) | `normal_mode_growth.html`, `salt_fingers.html`, `richardson_shear_instability.html`, `lorenz_attractor.html` (+ 20 earlier files) | `Viz.arrowPx(ctx, x0, y0, x1, y1, opts)` | fourth chapter defining it locally |
| 6 | **text-checksum row helper** for exact-string parity with a fluidpy text; `pyG` (Python `%g`-style formatting for mirrored strings) | reviewer: E4, E6; `pyG` in `salt_fingers.html`, `richardson_shear_instability.html` | `Viz.selftest.textRow(name, jsString, pyExpr)`, `Viz.pyG` (beside ch10's `pyRound`) | pins mirrored verdict strings without ndarray methods |
| 7 | **orthographic 3-D projector with drag-to-orbit** on a 2-D canvas | `lorenz_attractor.html` (`projector`) | the offline path of `Viz.three` (`Viz.three.canvas(view, {theta, phi})`) | 3-D line data with no CDN; Ch. 13 Ekman spirals |
| 8 | **`gutterPlot`** (fourth pass asking) | `kelvin_helmholtz_boundary.html` (+ 26 earlier files) | `Plot` option `gutter: true` | unchanged from the ch10 list |

`Viz.num.linalg` (ch10 rank 1) was **not** needed after all: no ch11 explainer solves an eigenproblem in the browser — E3
evaluates a 3 × 3 complex determinant, the others read parity-checked tables. Order for the library pass: ch11 ranks 1–4 with
the ch10 list; then `tools/viz_inline.py --all`, `tools/shot.py --chapter ch01 … ch11 --quick` and
`templates/viz_example.html --quick`. **Already done by the orchestrator in the ch11 viz phase** (not candidates any more):
L1 coloured `<b class="c-…">` in view titles and L2 `.viz-stage-badge:empty { display: none }` are in `assets/viz_base.css`
and inlined in every explainer and template; the five local copies of the badge rule and E6's legend workaround can be deleted
at the next edit of those files.

### ch11 table generators (`tools/viz_tables/ch11/`)
`benard_neutral_curve_sigma_table.py` (E3 growth table; paths relative to the script), `taylor_couette_onset_sigma_table.py`
(E5), `inviscid_shear_criteria_tables.py` (E7), `orr_sommerfeld_neutral_curve_table.py` (E8; continuation of one mode family,
Bickley on a 36 × 37 log-k grid; writes to `outputs/ch11/viz_tables/`). **To tidy**: the E5 and E7 scripts hard-code the
absolute repository path as `ROOT`, and the E3/E7 scripts write their JSON/JS next to themselves — make all four
`Path(__file__)`-relative with one output folder; **the E8 table is spliced into the HTML by hand** (a `tools/viz_embed_table.py`
that replaces a marked `const OS = …;` block would make it repeatable). Every embedded table is ours, labelled, and pinned by
parity rows against the fluidpy function.

### ch11 open explainer items (Should, non-blocking; from `reports/ch11_viz.md` round 2)
- E8 `orr_sommerfeld_neutral_curve`: rotated y title "wavenumber kL (log)" clipped by one character at 360×640 (→ "k L
  (log)"); quiz 4 answer "Its lower edge below k = 0.02, not resolved." (add "is"); status "grows already at Re = 500" for
  the jet far above Re_c (drop "already").
- E4 `salt_fingers`: 3 of 15 code lines wrap at 1366×768; the EQTAG split-string workaround stays until L6.
- E3 `benard_neutral_curve`: 1 of 41 code lines wraps.
- E1 `normal_mode_growth`: thin Bénard strip at 360×640 (draw fewer wavelengths).
- Backup B1 `period_doubling_route` not built.

### ch11 machinery / tooling TODOs (for the orchestrator; not done in this pass)
- **L3 `tools/shot.py`: fail on any `.katex-error` node and on control characters in rendered text, on every page of every
  pager** (E2 passed with two errors; a builder's probe script and the reviewer's sweep exist only in session scratchpads).
- **L4 deep links**: apply `#<param>=…` after `goStep(0)` (all nine explainers override the hash with step 1's `set`).
- **L5 code-line length**: lint `code[].src` lines longer than ≈ 44 characters, or shrink the code font one step.
- **L6 `tools/eq_refs.py`**: an allow-marker for strings that mirror a fluidpy text byte for byte.
- `tools/check_public.py`: it sees only tracked/staged files — document "run after `git add`" in the tool's help and the
  pre-push hook text (a selftest label matching a private book string was committed locally before being reworded).
- One private scratch subfolder per agent (collisions in ch10 and ch11); never pipe an audit through `tail` when its exit
  code matters.
- `viz_lint` control-character rule (now six chapters).

### ch11 lesson candidates for the skills (not promoted in this pass: the brief forbade library and skill edits)
- `interactive-viz` Lessons: (ch11) **embed a table of the tested fluidpy function instead of a rough in-page estimate**, and
  snap it to exact zero on the neutral curve the page computes live · interpolate so that growth ends on the exact neutral
  point; never bilinear across a neutral curve · follow one mode family by continuation when tabulating eigenvalues · pin
  mirrored fluidpy strings with checksum rows · derive every "the others decay" phrase from the computed signs · a status may
  be computed from the run · code lines ≤ 44 characters at 1366×768 · TeX strings are raw templates; sweep `.katex-error`
  before review · when a page disagrees with a fluidpy table near a boundary, report the solver — do not add "≈".
- `math-to-python` §7: (ch11) **an unbounded profile needs a box that scales with 1/k** (≥ 12 e-folds), proved by box and N
  doubling · **return the leading mode or raise; a 0 in a table has a stated meaning** · near a neutral point use a complex
  path round the critical layer and keep an independent route (shooting, determinant, compound matrix) · least-damped ≠ the
  mode to interpolate; neutral curves may have several bands; never join across NaN · eliminate boundary unknowns rather than
  replacing rows · spectral convergence is digits and has a round-off rise at large N · a new solver option needs its own
  test before it is documented (`semi_infinite`).
- `verify-implementation`: (ch11) **classifier and verdict functions: test both sides of every boundary and every degenerate
  input, and scan against an independent computation on a grid** · **pin sets and physics, not counts** (299 positives pinned a
  defect) and re-derive pinned numbers when the numerics change · design-document "expect" values are claims to test, not
  references · mutants for every "old behaviour" after a fix (61 planted, 1 equivalent survivor recorded with its reason).
- `teaching-style` Lessons: (ch11) **one-way theorems are taught with the word "necessary" or "sufficient" and a
  counter-example** · "trust what does not move with N" overlays · say what a table does not resolve next to its figure · a ★★★
  check cell tests the derived result itself · a numerical disagreement can be the lesson (RK4 vs DOP853 past the
  predictability time) · the ch11 derivation moves in `concept_map.md` (conjugate-multiply and take the imaginary part;
  substitute ψ̂/(U − c)ⁿ to make the equation self-adjoint; map a new problem onto a solved one; Galerkin truncation).
- `colab-notebook` / nbkit: (ch11) animations dominate the runtime (126–239 s by machine load; dpi 60 and 18–22 frames) —
  cache animation frames or cap frames in the builder; echo comments (`# draw <y> against <x>`) are still generated — fix the
  template; a live widget must catch `ValueError` from solvers that now raise.

### ch10 candidates (listed, **not promoted**: the ch10 knowledge pass ran while the site-publisher was reading `viz/`)
Found by the builders' notes and a script over the eight ch10 chapter scripts (the part after the inlined library); "files" =
ch10 files defining the helper locally (earlier-chapter counts from the ch09 list where the helper was already known).
Ranked by payoff for **Ch. 11 instability (eigenproblems, growth curves), Ch. 12 turbulence, Ch. 13 GFD (live shallow-water /
C-grid models)** and by risk (pure additions first):

| Rank | Helper | Found in (file · local name) | Files | Proposed library name | Recommendation |
|---|---|---|---|---|---|
| 1 | **dense and structured linear algebra** | E6 `luFactor`, `luSolve`, `pinnedLU`, `luSolveRefined`; E8 `luFactor`, `luSolve`, `solveLU`, `lstsqMinNorm`, `cholesky`, `lowerSolve`, `upperTSolve`, `jacobiEig`; E2 `implicitLU`, `luSolve`; E7 `bandLU`, `bandSolve`; E4 `gtsv`; E5 `thomas` | 6 (E2, E4, E5, E6, E7, E8) | `Viz.num.lu(A)` → {solve}, `Viz.num.lstsq(A, b)` (min-norm), `Viz.num.banded(A, kl, ku)` → {solve}, `Viz.num.thomas(a, b, c, d)` (alias `gtsv`), `Viz.num.cholesky`, `Viz.num.eigSym(A, B?)` (Jacobi; generalised via Cholesky) | **promote first, before Ch. 11** — every stability explainer needs eigenvalues and every live implicit/Poisson solve needs a factorisation; parity rows already exist in E2, E6, E7, E8 to move with them |
| 2 | **log axes and the narrow-view gutter** (third chapter asking) | `logAxes` E1, `logPlot` E4, `drawLogPanel` E7 (+ ch09 E5 log–log C_D); `gutterPlot` E2, E3, E4, E5 (+ 22 earlier files) | 3 + 4 (26 files with `gutterPlot` overall) | `P.axes({xlog, ylog})` with decade ticks and minor ticks; engine y-title gutter below ~420 px | **promote** as one engine-plot pass (ch07 ranks 3/8, ch08 rank 6, ch09 rank 5, now ch10) |
| 3 | **filled polygon in data coordinates** | `polyFill` E3, E4, E5 | 3 | `P.polygon(xs, ys, {fill, stroke, alpha})` (or `P.fill(poly)`) | promote with rank 2 (pure addition) |
| 4 | **matrix view** (heat-coloured cells, highlighted block, cell pick) | E5 `drawMatrix`/`cellColor`/`matTitle` (assembly), E6 Poisson-matrix and stencil views | 2 | `Viz.matrix(v, M, {cellColor, highlight: [[r0, c0, r1, c1]], labels, onPick})` | promote — Ch. 11 (discretised operators), Ch. 13 (C-grid operators) will draw matrices again |
| 5 | **linear colour fill of P1 triangles** | E8 `fillLinear` (+ `pColor`, `mix`, `rgbOf`) | 1 | `Viz.field.triFill(v, P, tris, values, cmap)` (barycentric / gradient fill) | promote when a second FE explainer appears (Ch. 14 low-Re bodies, Ch. 16) |
| 6 | **Python-compatible numbers for Code `live` and parity** | `pyRound` E3 (round-half-even), `rint` E3; `pn` (Python literal) E2, E3, E4, E5, E6, E7, E8 (+ 18 earlier) | 1 / 7 (25 overall) | `Viz.num.pyRound(x, nd)`, `Viz.pyNum(v, sig)` | promote with the ch09 formatter pass (rank 4 there) |
| 7 | **a time march with checkpoints** (scrub back without re-running from zero) | E3 `class Run` (`getRun`, `exactAt`); E7 `newRun`, `advance`, `stateAt`, `simulate` | 2 | `Viz.run(stepFn, state0, {dt, checkpointEvery, maxSteps})` → {at(t), advance(t)} | promote — every live PDE explainer (Ch. 11 saturation, Ch. 13 shallow water) scrubs a march |
| 8 | **background pump for heavy live solvers** (advance a solve in slices between frames, status "computing …") | E7 `pump`, `want`, `ensure` (setTimeout slices) | 1 | `Viz.pump(task, {budgetMs, onProgress})` + an engine status hook | promote with rank 7; fixes E7's disagreeing clocks by design |
| 9 | **derivation result page "key lines only"** (builders: the result page repeats the whole chain; long chains page 3–4× on phones) | E3 D18 (14 steps), E4 D16 (12), E5 D11 (13), E6 D20 (12) | — | `derivations[i].result.chain: 'key' \| 'all'` with `steps[k].key: true` | engine option (ch08 asked for a result-chain option too) |
| 10 | **per-layout view hiding and empty-row collapse** | `hidePortrait` used in 7 of 8 files; builders: hiding the only view of a row left an empty row (scoped CSS to collapse it); E6 wants n capped on portrait | 7 (+ 19 earlier) | engine `views[i].hideOn: ['portrait', 'landscape', 'short']`, collapse rows with no visible view, `params[k].maxOn` | **promote** (ch09 rank 2; the most requested engine feature since ch03) |
| 11 | **exact fractions** | E1 `class Fr`, `gcd`, `frPow` (+ ch01 E5 `Frac`) | 1 (+ 1) | `Viz.Frac` | promote — now two chapters (the ch01 rule "on second use") |
| 12 | **complex formatting and arithmetic** | E2 `ctex`, `ctxt`, `imZero`, `sg` (+ ch09 E6 `cmul`/`cdiv`/`cadd`/`cexp`, ch06–ch07) | 1 (+ 4) | `Viz.cx` with `Viz.cx.tex(z)` (drops ±0 i) | promote before Ch. 11 (ch09 rank 11) |
| 13 | **layout probes** | `phone()` and `lay()` in all 8 ch10 files (+ 16/14 earlier) | 8 | expose `g.app.layout` / density | **replace, don't promote** (ch09 rank 6) |

**Recommendation for the next library pass (orchestrator, when nothing reads `viz/`)**: the ch09 rank 1 lint rule (tool)
first; then ch10 ranks 1–3 (linear algebra, log axes + gutter, polygon) and 10 (`hideOn` + empty-row collapse) **before Ch. 11's
explainers**; then ranks 6–8 (formatters, `Viz.run`, `Viz.pump`) and 12 (`Viz.cx`) before Ch. 13; rank 9 with the ch08 pager
"now"-box slicing. After promotion: `tools/viz_inline.py --all`, then `tools/shot.py --chapter ch01 … ch10 --quick` and
`tools/shot.py templates/viz_example.html --quick` must all PASS (serially).

### ch10 open explainer items (Should, non-blocking; from `reports/ch10_viz.md` rounds 1–2)
- `fd_stencil_order`: quiz Q2 (x² vs x³) not answerable with the offered functions — add a polynomial chip or reword; D01
  step 6 *why* ends "…the third line of (10.6)" without the formula; orange stencil nodes vs amber time term too close in
  FTCS mode.
- `von_neumann_amplification`: Explain §3 for upwind could add the two-row table (C ≤ 1, C + 2β ≤ 1) with the failing edge
  lit. (The "±0 i" at θ = π and the β = 100 Equations clip were fixed in the orchestrator's round-2 resolution.)
- `upwind_cfl_advection`: a negative D_num at C > 1 (Code comment −2.632e-4 m²/s) needs "(negative: anti-diffusion, C > 1)".
  (The 1 px phone-land "Right now" clip was fixed: note shortened, three clean runs.)
- `cell_peclet_wiggles`: derivation header "EQS. (10.90)–(10.93)" is a bare range (the goal page could name the two schemes'
  equations).
- `fem_hat_assembly`: on phones the matrix hides its numbers for n ≥ 4, so tour step 4 ("click any matrix cell") relies on
  the inspector — keep the row-of-node-A numbers in the title on every step.
- `mac_projection_staggered`: 16² illegible on portrait phones — cap n at 12 there (storyboard said so); tour step 1 uses the
  default shear field, not the 8² preset (note in the design).
- `lid_driven_cavity`: MacCormack ◆ missing from the convergence legend and Explain §0; transport, view title and readout
  clocks disagree on phones while a live grid spins up.
- Backup B1 `operator_splitting_theta` not built.

### ch10 machinery / tooling TODOs (for the orchestrator; not done in this pass)
- `tools/viz_lint.py`: fail on any byte < 0x20 other than `\n` in an explainer and on `(?<!\\)\\(f|t|r|b|v)[a-z]` inside plain
  JS strings (**five chapters now**; ch10 lost `\approx` via a heredoc).
- `tools/eq_refs.py`: accept en-dash ranges on `ref:` badges ("Eqs. (10.90)–(10.93)" next to shown TeX) and scan template
  strings sentence by sentence (it missed E8's bare "(10.184)").
- `tests/test_machinery.py`: a slider-drag test (move a transport or parameter slider by keyboard with the Explain tab open,
  assert the Explain HTML changed) — the ch10 reviewer did it by hand in 25 cases.
- `tools/check_public.py`: `--untracked` to scan new files before they are added (the viz reviewer grepped untracked `viz/ch10/`
  by hand), and **fail — not skip — on an invalid `tests/book_values_*.json`** (today an invalid file silently disables the
  `_forbidden_public` scan).
- `tools/shot.py`: give `py:` rows `np.min`/`np.max`/`np.add.reduce` and ndarray methods (or document the restriction); reset
  presets before each tab audit or add a preset sweep (state-dependent clips); `--repeat 3` for flake detection under load.
- `tools/nbkit.py`: move `build_ch10.py`'s `self_check_ctrl` (control characters) and `self_check_near` (bare equation numbers)
  into nbkit for every builder; lint templated code comments ("# show the numbers computed above").
- `tools/embed_check.py`: the Derivation-tab vs `--dump` diff after the lesson review (ch09 TODO; ch10 matched by hand).
- Process: run the merge-gate pytest alone (≈ 10 min; > 1 h in parallel with notebook execution and shot audits); one
  scratchpad subfolder per agent.

### ch10 lesson candidates for the skills (not promoted in this pass: the brief forbade library and skill edits)
- `interactive-viz` Lessons: (ch10) **pin every mirrored fluidpy string by an exact-text parity row** (`stability_verdict`
  reasons went stale) · **audit tabs at the extreme presets** — clips can be state-dependent (β = 100 on the Equations tab) ·
  intermittent 1–26 px phone clips under load: leave a line of slack and run shot 3× · `py:` rows: `np.min(x)`, not `x.min()` ·
  Code-tab `{{placeholders}}` share one namespace across blocks — unique names · a live solver in the page needs a factorise-once
  matrix and a background pump, and its status must say when the picture lags the transport · per-layout parameter caps (n ≤ 12
  on portrait) · never write explainers through shell heredocs.
- `math-to-python` §7: (ch10) **a table written by the same function is not a reference** — test fitted quantities against an
  exact field (paraboloid) and must-hold inequalities (fitted min ≤ grid min), then regenerate every cached table and its
  consumers · **series stop on the term envelope**, never on the actual term (a term vanishes where its sine does) · **every
  explicit solver takes an explicit Δt rule, computes the full limit (all rates) and raises above it**; growth caps, never NaN ·
  build 1-D operators link by link from real neighbours and test null vectors on 1-cell and 2-cell grids · evaluate floating-point
  references by integer shifts · guard the degenerate case of every division (GCI p = 0) · scheme-dependent quantities take
  required arguments (`numerical_diffusivity(C=)`); boundary-like defaults (`lid=1.0`) default to None.
- `verify-implementation`: (ch10) observed orders only in the asymptotic range — state the grids (BTCS 1.80 → 1.98, upwind 0.76
  → 0.947); a benchmark saturates at its own error (MAC 128² vs Ghia) — use self-convergence as the convergence evidence;
  "benchmark" only inside the stated tolerance (MacCormack 64² 1.38 % is not); "converged" only with an asserted order
  (removed from six functions); mutants in code reached only through cached heavy runs are invisible — keep one slow test that
  recomputes a coarse cached run to 6 s.f.; derive stability regions in closed form and confirm by a brute-force scan.
- `teaching-style` Lessons: (ch10) **pre-asymptotic honesty** (a failed convergence study taught as failed, p = 0.43) · a
  **conventions table up front** for overloaded letters (i, α, β, θ, g, D, M, K, S, n, R, δ, F, G, p) and enforce it in every later
  cell · name printed slips "#k" so they never read as recap IDs · **DEVIATION boxes** with both numbers computed · drawings
  must draw what the reading note says · a named test case (Taylor–Green) gets a primer before first use · derive the discrete
  solution, not only the continuous one (rʲ) · the ch10 derivation moves in `concept_map.md` (Taylor-cancel stencils, modified
  equation via the PDE, one Fourier mode ⇒ G, s = sin²(θ/2) endpoint test, CFL as the characteristic's foot, test function +
  parts, unit coefficient vectors, parent element + scatter-add, geometric trial rʲ, divergence of the correction ⇒ Poisson,
  cell flux budget without a pressure BC, planted (−1)^{i+j}, three-grid Richardson).
- `colab-notebook` / nbkit: (ch10) `self_check_ctrl` + `self_check_near` for every builder; templated code comments are a
  review failure — the builder template should generate meaning + unit comments; notebook 72 s with eight explainers, four
  animations, seven plotly figures and cached heavy runs (the < 5 min budget held because 128² and unsteady FE runs are
  script-only).

### ch09 candidates (listed, **not promoted**: the ch09 knowledge pass ran while the site-publisher was reading `viz/`)
Counts = explainer files (all 71, ch01–ch09) defining the helper locally outside the inlined library, counted by script
(`function name(` / `const name = (…) =>` in the chapter script); "(ch09 n)" = how many of the nine ch09 files. The one
library change of ch09 (app.set rebuilds Explain on user-driven transport changes) is already in `assets/viz_lib.js`.
Ranked by payoff for **Ch. 10 CFD, Ch. 11 stability, Ch. 12 turbulence, Ch. 13 GFD** and by risk (pure additions first):

| Rank | Helper | Found in (file · local name) | Count (ch09) | Proposed library name | Recommendation |
|---|---|---|---|---|---|
| 1 | **control-character guard** (tool, not library) | four chapters shipped `\f`/`\t`/`\r` inside TeX (ch04 E7, ch06 E2, ch09 E3, E5); ch09's E5 bytes were in a `String.raw` string, so the source text itself was corrupted | — | `tools/viz_lint.py`: fail on any byte < 0x20 other than `\n` in an explainer, and on `(?<!\\)\\(f\|t\|r\|b\|v)[a-z]` inside plain JS string literals | **do first** (tool change; zero risk for existing files — a scan of all 71 finds none today) |
| 2 | **per-layout view hiding** (viz-builder request: `views[i].hideOn: ['short']`) | ch09 builders asked for it (E3's 80 px views at 360×640 need the `curve` view hidden on the smallest portrait only); `hidePortrait` + scoped `<style data-chapter>` rules in every ch09 file; ch03 rank 2, ch04 rank 4 (19 files then) | 9 (9) + 19 earlier | engine `views[i].hideOn: ['portrait', 'landscape', 'short']`, `rows[i].hideOn`, `presets.hideOn` (Derivation tab exempt, as the scoped CSS does) | **promote** — the most requested engine feature since ch03; removes chapter CSS |
| 3 | **pixel-space arrow** | `arrowPx` E1, E2, E4, E5, E6, E7, E8, E9 (+ 27 earlier) | 35 (8) | `Viz.arrowPx(ctx, x0, y0, x1, y1, {color, width, head, minLen})` | **promote** (ch07 rank 2, ch08 rank 3; now 35 files) |
| 4 | **number and length formatters** | `f3` 44 (7), `t3` 36 (6), `t4` 33 (4), `f4` 23 (4), `pn` 18 (7, Python literal for Code `live`), `tmm` 3 (3) + `mm` 5 (3) (mm/µm switch) | 44 / 36 / 33 / 23 / 18 / 8 | `Viz.fixed(v, d, {tex})`, `Viz.pyNum(v, sig)`, `Viz.fmtLen(m)` | **promote** (ch08 rank 4) |
| 5 | **narrow-view plot gutter and log axes** | `gutterPlot` E1, E2, E4, E5, E6, E7, E8 (+ 15 earlier); E5's log–log C_D(Re) and log-Re dial | 22 (7) | engine y-title gutter below ~420 px, `P.axes({xlog, ylog})` | **promote** as one engine-plot pass (ch07 ranks 3/8, ch08 rank 6) |
| 6 | **layout probes** | `phone()` 7 ch09 (16 total), `lay()` 7 ch09 (14 total) | 16 + 14 (7 + 7) | expose `g.app.layout` / density (`Viz.isPortrait(app)`) | **replace, don't promote** — a window test disagrees with the engine inside a notebook iframe |
| 7 | **similarity-profile panel** (profiles at several stations, raw vs rescaled, legend in the title) | `drawProf` E2, E3, E7, E8, E9 (5, **all ch09**), `profTitle` E2, E7, E8, `drawJet`/`jetTitle` E7, E8; ch08 E5/E6 stage the same idea | 5 (5) | `Viz.profiles(v, {curves, scale: 'raw'\|'rescaled', marks})` | **promote after Ch. 12's jets/wakes** confirm the API (ch08 said "wait for Ch. 9's Blasius explainer" — five files now; the shape is clear) |
| 8 | **label with a background box** | `lab` E1, E2, E4, E5, E6, E7, E8 | 22 (7) | `Viz.label = (ctx, s, x, y, o) => Viz.text(ctx, s, x, y, {bg: true, size: 12, …o})` | small; with rank 4 |
| 9 | **moving hatched wall** | `hatch` E1, E2, E4 (+ 3 earlier; ch08 rank 7 counted 8 variants) | 6 (3) | `Viz.hatch(ctx, x, y, w, h, {step, shift, color, angle})` | **promote** — every wall-bounded flow draws one |
| 10 | **table interpolation and quadrature** | `hermite` E4, E8 (cubic Hermite on embedded fluidpy tables), `simpson` E7 (+ ch08 E4), `sech2` E6, E7, E8 | 2 / 2 / 4 (2 / 1 / 3) | `Viz.num.hermite(xs, ys, dys)` or `Viz.num.pchip`, `Viz.num.simpson`; `Viz.num.sech2` trivial | promote with the ch08 numerics pass (erfc precision, gamma, golden) |
| 11 | **complex arithmetic** | `cmul`, `cdiv`, `cadd`, `cexp` E6 (+ ch06 E5/E6, ch07 E5) | 4 / 4 / 3 / 3 (1) | `Viz.cx` (ch06 rank 1) | promote before **Ch. 11** (normal modes, Orr–Sommerfeld spectra) |
| 12 | **pager slicing of a derivation quote in a walkthrough card** | E2 phone step 2: "All steps →" cut by "continued on the next page" | — | treat the quote card's footer as keep-with-previous (library bug 1) | library fix with the ch08 rank-12 "now"-box slicing |

**Recommendation for the next library pass (orchestrator, when nothing reads `viz/`)**: rank 1 (lint rule) now; then the
ch08 ranks 1–2 (pager-label CSS, erfc precision — a defect) and ranks 2–5 here (per-layout hiding, arrow, formatters, plot
gutter + log axes); then `Viz.cx` before Ch. 11 and `Viz.waves` before Ch. 13 (ch07 list); `Viz.profiles` after Ch. 12.
After promotion: `tools/viz_inline.py --all`, then `tools/shot.py --chapter ch01 … ch09 --quick` and
`tools/shot.py templates/viz_example.html --quick` must all PASS (serially if memory is short).

### ch09 open explainer items (Should, non-blocking; carried from `reports/ch09_viz.md` round 2, plus two found in this pass)
- **Derivation-tab drift (found by the knowledge pass; the ch08 rule makes it Must-level for the next viz touch):**
  `blasius_similarity_collapse` D06 = 10 steps vs the notebook's 14 (lesson round 1 split step 6 into five moves: f f″ =
  f(f′ − 1)′ → integrate (9.27) → by parts → solve for f″(0) → η → ∞; the explainer also has an extra first step "Translate
  the conditions"); `karman_street_stability` D14 = 15 vs 16 (notebook step 10 "Conjugate the lattice sums" missing). Rebuild
  both from `build_ch09.py --dump`, renumber walkthrough `derive: {id, step}` links. **RESOLVED before publish:** D06 14/14,
  D14 16/16, links renumbered, shot PASS; karman_street_stability should-fixes (lost `\Gamma`, 31 vs 61 wavenumbers, ±Γ
  convention stated) also done. Note `build_ch09.py --dump` truncates — builders rebuilt the reference from the `PF` data.
- `bl_scaling_thicknesses`: `2\,tau_0\,x` lacks a backslash (Explain §5); amber δ* vs orange θ too similar (use `#d9a300` or
  hatch one area); phone u/U title clipped (drop "δ*" on phones).
- `blasius_similarity_collapse`: rescaled collapse drawn from one function, not from each station's data; phone step 2
  slices the derivation quote (shorten or drop the quote on phones).
- `falkner_skan_family`: at 360×640 the two views are ~80 px and the profile's y title clips (hide `curve` on the smallest
  portrait — rank 2 above — or one row of preset chips).
- `thwaites_marching`: design note says 60 closure rows, the file (correctly) embeds 64 (text only).
- `cylinder_drag_crisis`: Explain §0 not mode-aware (sphere mode); "φ_s = 125°" label overlaps the U arrow in the notebook
  frame.
- `karman_street_stability`: inspector `'lower row, $+\Gamma$'` in a plain string renders "+Gamma" (`'$+\\Gamma$'`); Explain
  §0 says 31 wavenumbers, §6 says 61 (61 is right); +Γ on the upper row in D14 step 1 vs −Γ in the picture and code.
- `free_jet_similarity`: "inflow v = 59.3 mm/s" collides with "y stretched ×40" (laptop/notebook); `x^{-1/3}` wraps between
  base and exponent in tour step 3.
- `wall_jet_invariant`: "peak f′ = 0.07875 at η = 8.11" overprinted by the "f∞ = 1 (dotted)" legend; phone bars clip "1.7".
- `teacup_secondary_flow`: blue u(z) curve without a scale and not named in Explain §0; phone cup labels overprint the
  swirl symbols (hide below ~200 px); unlabelled grey dashed curves; plane (6.2) cited for the axisymmetric loop.
- Backup B1 `ball_swing_magnus` not built.

### ch09 machinery TODOs (for the orchestrator; not done in this pass)
- `tools/viz_lint.py`: control-character rule (rank 1 above).
- `tests/test_machinery.py`: open the template (or a ch09 explainer), move the transport slider by keyboard, assert the
  Explain HTML changed — the only guard against a return of R2-M1 (`shot.py` sets states and cannot see it).
- `tools/embed_check.py` (or the notebook phase gate): diff `app.cfg.derivations` step counts/titles against
  `build_chNN.py --dump` **after the lesson review**, not only during the viz review (ch09 D06, D14 drifted after it).
- Tests: cap solver iterations (`solve_bvp(max_nodes=…)`) in functions a mutation run exercises — a planted FS sign mutant
  hung for 900 s instead of failing.
- Whole-tree `pytest -q` stalls near 96 % when the browser machinery test runs after the long ch09 script test: run
  `tests/test_machinery.py` separately or mark the script test `slow`.
- `tests/book_values_ch09.json` was invalid JSON (`2.0/3.0` literals) until the verifier fixed it: add a JSON-validity check
  for `tests/book_values_*.json` to `tools/check_public.py` or a machinery test.
- Docstring: `karman_street_spectrum` calls the eigenvalue σ (notebook λ; D14's σ is sinh/cosh²).

### ch09 lesson candidates for the skills (not promoted in this pass: the brief forbade library and skill edits)
- `interactive-viz` Lessons: (ch09) **the transport parameter is often physical** (x, n, Re, a cursor) — the engine now
  rebuilds Explain on every user `set()`; an explainer must never assume "transport = time" · **lost backslashes → control
  characters** (fourth chapter): `String.raw` or doubled backslashes, and scan the written file for bytes < 0x20 ·
  **book value first, exact value second, both labelled** in status, Explain, chips and presets (Example 9.2: λ = −0.090 →
  0.158, then −0.068148 → 0.126) · re-dump the notebook's derivations after every lesson round (D06, D14 drifted after the viz
  review) · tighten parity rows to the achieved agreement (1e-9, not 1e-4, for a copied table value) · defaults whose numbers
  mean something (U = 1 m/s: an audible street) · one sign convention across derivation, picture and code · break the
  invariant on purpose (wrong-exponent toggle) to show why the exponents are forced.
- `math-to-python` §7: (ch09) **a test on |x| cannot see the sign of x** — pin every signed transcription with a manufactured
  field ((9.8) was negated and passed an |δp| ∝ 1/Re test) · a fold is mathematics: parametrise by the quantity that passes
  through it (f″(0)), not by the parameter that turns back · prefer a scaling symmetry (Töpfer) or `solve_bvp` with
  continuation to long shooting (sensitivity e^{η²/4}); give shooting brackets that scale with the solution · truncation η_max
  per member (thick layers near a fold; the wall jet's free scale) · never type a constant a function can compute
  (`LAMBDA_SEP_FS`), keep the book's as a separate named constant · defaults encode physics (`cp_base=None` gave 2.59) · sum
  lattice series in ±n pairs and assert the truncation order · cap solver iterations so wrong variants fail, not hang.
- `verify-implementation`: (ch09) a mutant that hangs is not "killed" — record it; the reviewer re-ran 6 of 25 mutants and
  found the (9.8) sign survivor: re-run a sample of the verifier's mutants independently · "slip" is a claim about the page:
  read the sentence before calling it one (R5 was a trap), and read nearby sentences for unrecorded slips (R17) · labels:
  "approximate" for few-percent agreement, "secondary" for encyclopedia correlations, a course handout is not V5 ·
  report a bracket, not a number, where the result is inlet-sensitive (marched separation near the fold).
- `teaching-style` Lessons: (ch09) **two criteria, the book's first**, both computed, in the same order everywhere · **every
  derived number in a reading note is computed by a cell** (a new note said the C_p curves cross at 110°; it is 132°, from
  1 − 4 sin²φ = −1.2) · a tool first used inside a derivation gets its primer before it (integration by parts, P218a) and a
  step that needs four moves is four steps (D06) · a shaded gap between two curves is a force only with its weight (cos φ) ·
  a sign-checking tiny example for every transcribed equation whose sign matters · keep "slip" (printed-vs-correct box +
  failing code option) and "trap" (⚠️ callout) apart · the ch09 derivation moves in `concept_map.md`.
- `colab-notebook` / nbkit: (ch09) a builder self-check that numbers in reading notes (angles, percentages) appear in a
  printed cell output; a coverage check that a derivation's *Tools* line names only primers placed before it; notebook
  ~94 s with nine explainers and six animations.

### ch08 candidates (listed, **not promoted**: the ch08 knowledge pass ran while the site-publisher was reading `viz/`)
Counts = explainer files (all 62, ch01–ch08) defining the helper locally outside the inlined library, counted by script
(`function name(` / `const name =` in the `Viz.app` script, plus `<style data-chapter>` rules); "(ch08 n)" = how many of
the nine ch08 files. Ranked by payoff for the chapters ahead (**Ch. 9 boundary layers** — similarity profiles, erf/erfc,
log axes; Ch. 10 CFD; Ch. 11; **Ch. 13 GFD**) and by risk (pure additions first):

| Rank | Helper | Found in (file · local name) | Count (ch08) | Proposed library name | Recommendation |
|---|---|---|---|---|---|
| 1 | **pager label on one line** | chapter CSS `.viz-pager-label { white-space: nowrap; … }` — ch08 E1, E3–E9 (8), ch07 ×5 (rank 1 there too) | 12 (8) | one rule in `assets/viz_base.css` | **promote first** (CSS only, backwards compatible; two chapters now carry it) |
| 2 | **double-precision erfc / erfcinv** (continued fraction above 3, series below, Newton inverse; `erfHP`, `erfinvHP` beside them) | `stokes_first_problem` `erfcHP`, `erfcinvHP`; `similarity_exponents` `erfcHP`, `erfHP`, `erfinvHP`; ch04 `navier_stokes_term_balance` `erfcHi` | 3 (2) | fix `Viz.num.erf`/`erfc` to ~1e-15 and add `Viz.num.erfinv`, `erfcinv` | **promote** — `Viz.num.erfc`'s ~1.2e-7 is a *library defect* (too coarse for 1e-12 parity rows); Ch. 9 (Blasius vs erfc, temporal layers) and Ch. 13 (spin-up) will need it again |
| 3 | **pixel-space arrow** | `arrowPx` ch08 E2, E3, E5, E7, E8, E9 (+ 21 earlier); `pxArrow` ch08 E1, E4 (+ 3 earlier); `dotPx` 5 earlier | 27 + 5 (6 + 2) | `Viz.arrowPx(ctx, x0, y0, x1, y1, {color, width, head, minLen})`, `Viz.dotPx` | **promote** (ch07 rank 2, now 32 files) |
| 4 | **number formatters incl. a Python-literal formatter and a length switch** | `f3` 37 (9/9), `t3` 30 (9/9), `t4` 29 (8), `f4` 19 (5); `pn(v, k)` Python literal for Code-tab `live` (E3, E5, E7, E8, E9, E4 — three variants of the exponent rule; `pyNum` 8 earlier files); `fmtLen`/`fmtLenM`/`texLen` µm/mm/cm/m switch (E2, E4, E6; E5 inline) | 37 / 30 / 29 / 13 + 8 / 6 (9 / 9 / 8 / 6 / 3) | `Viz.fixed(v, d, {tex})` (document `sig` shorthands), **`Viz.pyNum(v, sig)`** (one exponent rule: `1e-06` style), **`Viz.fmtLen(m)`** beside `Viz.fmtTime` | **promote** — pure additions; `pyNum` removes three inconsistent Python-literal rules |
| 5 | **layout probes** | `phone()` 9/9 ch08 (+1), `lay()` 7 ch08; ch07 `narrow()` 9/9 | 10 + 7 (+ 34 `narrow`) | read `g.app.layout` / density, or add `Viz.isPortrait(app)` | **replace, don't promote** — a window test disagrees with the engine inside a notebook iframe; expose the engine's layout |
| 6 | **narrow-view and log–log plots** | `gutterPlot` ch08 E3, E5, E7 (+ 12 earlier); `logPlot` E8, `logFrame` E9, `logAxes`/`logTicks` 9 earlier | 15 + 11 (3 + 2) | engine y-title gutter below ~420 px (`v.plot({ylabelGutter: true})`); `P.axes({xlog, ylog})` with decade ticks | **promote** as one engine-plot pass (ch07 ranks 3 and 8) |
| 7 | **moving hatched wall** | `hatch` E1, E5; `hatchShift` E7 (moving with the wall); inline in E3; earlier `hatchPattern` (ch01 E3), `hatch` (ch03 E3), `hatchRect` (ch03 E7, ch04 E1) | 8 (4) | `Viz.hatch(ctx, x, y, w, h, {step, shift, color, angle})` | **promote** — every wall-bounded flow of Ch. 9–11 draws one |
| 8 | label with a background box | `lab` 31 files (9/9) — `Viz.text(…, {bg: true, size: 12})` wrapped | 31 (9) | document the idiom (`Viz.label = (ctx, s, x, y, o) => Viz.text(ctx, s, x, y, {bg: true, size: 12, …o})`) | small; with rank 4 |
| 9 | Γ function | `gammaFn` (Lanczos) E4, E6 (Huppert η_N) | 2 (2) | `Viz.num.gamma`, `lgamma` (ch01 E1 has `lgamma`) | promote with rank 2 (numerics pass) |
| 10 | golden-section maximiser | E3 (slider optimum, parity 1e-6 vs Brent) | 1 (1) | `Viz.num.golden(f, a, b, {max, tol})` beside `brentq` | promote with rank 2 — optima recur (Ch. 11 most unstable k, Ch. 13 maximum growth) |
| 11 | **derivation-result "chain" option** | E1, E3, E5 result pages repeat every line and page on phones | — | `derivations[i].result.chain: 'short' \| 'full'` (show the first and last lines on phones) | library feature; shortens result pages |
| 12 | **pager unit-box placement** | E3 D14 s14, E8 D29 s5: a tall "now" box placed as a unit leaves ~450 px blank on phone-tall | — | slice the "now" box at `aligned` rows like other atoms, or pack "now" first when "we had" is long | **library fix** (known, not counted by the reviewer) |

**Recommendation for the next library pass (orchestrator, when nothing reads `viz/`)**: rank 1 (CSS) and rank 2 (erfc
precision — a defect) first, then ranks 3–4 (arrow, formatters incl. `pyNum`, `fmtLen`), then the engine work (rank 5
layout exposure, rank 6 gutter + log axes, rank 12 pager slicing, rank 11 chain option), then hatch, gamma, golden.
Together with the ch07 list (`Viz.waves` and per-mode hiding **before Ch. 13's explainers**). After promotion:
`tools/viz_inline.py --all`, then `tools/shot.py --chapter ch01 … ch08 --quick` and
`tools/shot.py templates/viz_example.html --quick` must all PASS (run shot serially if memory is short).
JS stages the curator flagged as possible helpers (`Viz.profiles`: profiles at several times with a raw/rescaled mode, E5
and E6; `Viz.channel`: a channel with profile arrows and tracers, E1 and E3) stayed chapter-specific — wait for Ch. 9's
Blasius explainer before abstracting.

### ch08 machinery TODOs (for the orchestrator; not done in this pass)
- `tools/nbkit.py` / builders: **automatic equation insertion (`show_eqs`) must skip mentions of the book's printed form**
  ("printed", "book prints", "the text") — it showed the corrected equation where the text named the slip (lesson M2).
- **Post-processing regexes must never run inside maths they just created** (`tidy_raw_tex` turned `$e^{…}$` into
  `$e^(…)$` in 20 places, lesson M1); add a build-time guard (no `$…^(…$`) to `self_check_prose` — move both into nbkit.
- `_title_no_numbers` deleted "(8.29)", "(8.39)" from derivation headings instead of writing the equation (lesson M8).
- ch08 scripts need `--no-show` (plt.show() blocks headless runs); `tests/test_ch08.py::test_scripts_V1_*` passes it.
- `tools/shot.py`: run serially when memory is short; the viz reviewer ran each file alone and confirmed flaky clips with
  two consecutive runs.
- A scripted comparison of `app.cfg.derivations` step titles/counts with `build_chNN.py --dump` (the ch08 reviewer's
  Playwright script) belongs in `tools/embed_check.py` or `shot.py`.
- API naming: `slider_bearing_state["inlet_backflow"]` → add `wide_end_backflow` (alias) in a machinery pass.

### ch08 lesson candidates for the skills (not promoted in this pass: the brief forbade library and skill edits)
- `interactive-viz` Lessons: (ch08) **the notebook is the reference for derivation step counts** — build Derivation tabs
  from `build_chNN.py --dump`, not from design Part F, and re-dump after every lesson round (four explainers were one step
  short) · live rows use text-size integrals (display-size ∫ caused an intermittent clip) · parity rows cover every regime
  the sliders allow (α < −½ was reachable but unpinned) · compute every preset's regime numbers before building (8 Pa/m
  already past zero net flow; ice 10¹³ Pa s; Re = 2 not creeping) · draw the rejected branch of a solution (the growing
  root) so a "because bounded" step is seen · an independent JS brute-force of a regime flag is a second route that can
  catch fluidpy · run shot serially under memory pressure and confirm clips twice.
- `math-to-python` §7: (ch08) test every parameter at a non-zero value with both signs (the extremum sign hid behind
  U = 0) · derive a regime flag from the field and brute-force it over the full parameter range (slider recirculation at
  either end) · every asymmetric solution gets a direction test (Oseen wake) · compare magnitudes of vector terms, not
  one component (far-field prefactor ½, not ⅛) · keep the book's sign in the argument name (`dpdx`) with the earlier
  chapter's as a keyword alias (`G=`) · code the corrected form and keep the printed one as a named option
  (`form="book"`, `model="book"`, `printed_8_13b=True`) · Crank–Nicolson needs backward-Euler start-up steps after a jump
  (pure CN stalls) · degenerate nonlinear diffusion: implicit Newton + precursor film + conservative faces.
- `verify-implementation`: (ch08) a verifier "correction" is itself a claim — the implementer's independent magnitude
  computation overturned the loop-1 "⅛" prefactor; record both and pin the settled value with sympy and a numeric
  asymptote · write the finite-difference round-off budget into the test (ε r/h_rel² tripled by Richardson) before
  asserting · "converged" needs an asserted order; Gauss–Legendre exact from n = 1 is analytic; a correlation bracket is
  benchmark evidence for the bracket only · plant the old code of every review fix as a wrong variant (41/41).
- `teaching-style` Lessons: (ch08) a slips table up front (printed | corrected | where) plus ghosts and failing code
  options · three-way parity cells across chapters (ch03 developing profile far downstream = ch04 exact solution = ch08
  derived) · **measure every "~" claim** (prefactor, crossover, direction) · flag inserted derivation steps ("the book
  skips this move; we add it") · an approximation that is not required by the theory is labelled as extra, with its cost
  (the slider's linear load +15.6 % at α = 0.1 is not lubrication's assumption) · the ch08 derivation moves in
  `concept_map.md`.
- `colab-notebook` / nbkit: (ch08) `show_eqs` skips "printed" mentions; regex tidies run outside `$…$` only; headings keep
  their equations; notebook ~196 s with nine explainers.

### ch07 machinery TODOs (for the orchestrator; not done in this pass)
- `tools/shot.py` `py:` namespace: the expressions have no builtins and no `np`, so a lambda inside a row must bind numpy
  through default arguments (`H_fn=lambda X, th=ch07.np.tanh, hy=ch07.np.hypot: 20.0*th(…)` in E6's island rows). Either
  add `np` to the evaluation namespace or document the idiom in the `interactive-viz` skill.
- Pager: re-pack after KaTeX has rendered (or reserve one line of slack) — E7's intermittent 3 px clip came from late
  typesetting in a card packed to within a few pixels.
- Rule 9: a scripted scan of explainer strings and notebook markdown against the private `tests/book_values_chNN.json`
  (the ch07 viz- and lesson-reviewers did it by hand and still found the minimum group speed in three places) — e.g.
  `tools/check_public.py --book-values` run locally before review.
- `notebooks/build_ch07.py` has a `tidy_raw_tex` pass plus a build-time guard for raw `e^{…}` in plain *why* lines —
  move it into `tools/nbkit.py` so every builder gets it.
- `tools/eq_refs.py` still reports `ref:` badges on the Derivation header and Code tab (3 accepted hits in E9) — exempt
  them like the other `ref:` labels.
- `ch04.linear_wave_surface` overflows for kH ≳ 710 (verification O8): re-point it to `core.waves.cosh_over_sinh`.

### ch07 lesson candidates for the skills (not promoted in this pass: the brief forbade library and skill edits)
- `interactive-viz` Lessons: (ch07) **never type a number a function can print** — prose, checks and legends take their
  numbers from the mirrored function, and a computed value that equals an exercise answer at 3 s.f. is printed at 4 s.f.
  (rule-9 near-miss twice in E3/E5) · physical space + wavenumber space on portrait phones = two square views side by
  side at full height (E9) · one Code block per geometry, traced off-screen, so live comments never depend on the
  displayed mode (E6) · per-mode hiding via a data attribute + CSS until the engine has it (`markMode`) · `py:` parity
  lambdas bind numpy by default arguments · cards packed near the page height need one line of slack (KaTeX renders
  late) · a status that counts ("rays on the lee side") uses a strict geometric criterion · sign-safe presets for
  formulas the book prints for k > 0, with a parity row per quadrant.
- `math-to-python` §7: (ch07) complex-step derivatives need complex-safe analytic continuations of every |·|
  (|k| → −k where Re k < 0; `np.abs` silently falls back and loses precision) · keep speed and signed derivative apart
  (`phase_speed` |c|, `group_velocity` sgn(k)·c_g) · code the gradient form (c_g = ∇_K ω of N|k|/K) when a printed
  formula assumes a sign, and test all quadrants · overflow-safe hyperbolic ratios (e^{kz}(1 + e^{−2k(z+H)})/(1 −
  e^{−2kH}), x/sinh x via expm1) and `np.errstate(over="ignore")` around numpy's complex tanh · invert a monotonic
  dispersion relation with brentq between its two asymptotic roots and assert the residual · state asymptotic orders as
  relative or absolute (O(ka) relative = slope 2 absolute at fixed k) · two published "g′" differ by ρ₂/ρ₁: name the
  reference density in the keyword (`ref="lower"`).
- `verify-implementation`: (ch07) quote a benchmark's error bound from the **primary** paper at its printed precision
  (Fenton & McKee 1.7 %, Guo 0.75 % for β = 2.4908), not from a later note (1.5 %, 0.7 %) · a literal exercise set-up
  can give a different coefficient than the published result (Stokes γ = 3/8 vs 1): check the expansion's consistency
  order before attributing · plant sign-convention variants at k < 0 (six `_abs_k`/`group_velocity` variants, all caught
  by two tests) · long drift integrations: DOP853, rtol 1e-10, ≥ 10 periods, release at the mean depth · an independent
  numerical route for a closed form (FD sloshing eigenproblem → (7.65) at order 1.98) settles a book-value discrepancy
  (Exercise 7.6, +0.72 %).
- `teaching-style` Lessons: (ch07) make a linearisation measurable (residual scan on log axes; relative vs absolute
  slope) · show a printed slip as a ghost that visibly fails (frozen envelope, failing residuals) · "take it on trust
  for now (derived in C03–C04)" is an honest forward use · size the dropped terms from the forcing (the surface rises at
  aω) · compute crossovers and thresholds instead of quoting them (Ursell 4π²/9) · captions count what they claim
  (rays reaching the lee printed by the cell) · climate numbers carried through (38 kW/m swell, 198 m/s tsunami, 0.99 vs
  22 m/s thermocline vs surface wave) · idiom primers whose two-line demo output *is* the explanation (P179–P184) · the
  ch07 derivation moves in `concept_map.md`.
- `colab-notebook` / nbkit: (ch07) a `tidy_raw_tex` pass + build-time guard for raw TeX in plain text; the full run takes
  ~189 s (ch06 ~170 s) — set timeouts from a loaded run.

### ch06 candidates (listed, **not promoted**: the ch06 knowledge pass ran while the site-publisher was reading `viz/`)
Counts = explainer files (all 44, ch01–ch06) that define the helper locally outside the inlined library, counted by
script (`function name(` or `const name = (…) =>`); "(ch06 n)" = how many of the nine ch06 files. Same procedure as the
ch02 list, plus `tools/shot.py --chapter ch06 --quick`. Ranked by payoff for the chapters ahead (Ch. 7 waves, Ch. 9
boundary layers, Ch. 10 CFD, Ch. 13 GFD, Ch. 14 aerodynamics):

| Rank | Helper | Found in (file · local name) | Count (ch06) | Proposed library name |
|---|---|---|---|---|
| 1 | **complex arithmetic with numpy's branch conventions** (principal `sqrt` with signed zero on the cut, principal `log`, rotatable cuts) | `cx`, `cadd`, `csub`, `cscale`, `cmul`, `cdiv`, `cabs`, `cexp`, `csqrt`, `carg` in E4 `complex_potential_corners`, E6 `conformal_joukowski`; `cmul`/`cabs`/`csqrt` in E5; `cdiv` in E1; `logBranch`, `powerBranch` in E4 | 4 (4) | `Viz.cx = {c, add, sub, scale, mul, div, abs, arg, exp, log, sqrt, pow, conj}` + `Viz.cx.logBranch(z, cut)`, `powerBranch(z, p, cut)` mirroring `core.potential.log_branch/power_branch` (parity rows); **Ch. 7 (complex wave amplitudes), Ch. 11 (complex growth rates), Ch. 14 (airfoil maps)** |
| 2 | **arrow in pixel coordinates** | `arrowPx` ch04 E4; ch05 E2, E4, E5, E6, E7, E8; ch06 E1, E2, E4, E5, E6, E8, E9, E3 | **15** (8) | `Viz.arrowPx(ctx, x0, y0, x1, y1, {color, width, head, minLen})` — top by count since ch05 |
| 3 | **narrow-view plot with a y-title gutter** | `gutterPlot` ch05 E2; ch06 E2, E3, E4, E6, E8, E9 | **7** (6) | engine default below ~420 px, or `v.plot({ylabelGutter: true})` — the ch01/ch05 phone Must-fixes, now copied six times |
| 4 | **number formatters** | `f3` ch06 9/9 (19 overall), `t3`/`f4`/`t4` ch06 8/9, `fz` 4, `tp` 4, `fsci` 2 | ch06 **9/9** (all 44) | `Viz.fixed(v, d, {tex, keepZeros, roundoff})`, `Viz.sci` (unchanged top candidate since ch04) |
| 5 | **contour lines that skip branch cuts and singular rows** | `maxJump` option of E1's contour wrapper (7 uses); E4 masks the cut; E8's ψ = 0 hook at the axis (open) | 3 (3) | `Viz.field.contour(f, {maxJump, mask, clip})` — drop a segment whose end values jump by more than `maxJump` (a cut) or that touches a masked cell |
| 6 | **per-mode control/readout hiding** | `applyMode` ch06 E3, E4, E7, E9 (13 overall); `applyVisibility` E1, E8 | 15 (6) | engine: `params[k].modes` (quirk Q4); also the cure for E1/E5's 15–16 Explore pages |
| 7 | **scalar-field image from a callable** | `washImage` ch05 E2, E4, ch06 E3; `putImageData` wash in E1, E3; `heatImage` ch04 E8 | 5 (2) | `Viz.field.image(P, f, {vmin, vmax, cmap, diverging, alpha})` |
| 8 | **log axes with decade labels** | `logTicks` ch05 E3, E6, E9, ch06 E8; `logAxes` ch04 E9, ch05 E4, E6; E7's inline "decade grid + labels"; E4's log–log tip-speed view | 8 (3) | `P.axes({xlog, ylog})` with decade ticks 10ⁿ |
| 9 | term / waterfall bars | `drawBars` ch06 E5, E8 (13 overall); `barsTitle` 2 | 13 (2) | `Viz.bars(v, items, {measured, faint, collapseZeros, minRowPx})` — add E5's faint/solid pair per item |
| 10 | Gauss–Legendre nodes and a force quadrature on a sphere | `glNodes` ch05 E3, E6, ch06 E9; `forceQuad` E9 | 3 (1) | `Viz.num.gaussLegendre(n, a, b)` |
| 11 | small dense linear algebra | `luSolve`, `cond1` ch06 E8 | 1 (1) | `Viz.num.lu(A)`, `.solve(b)`, `Viz.num.cond1(A)` — **Ch. 10 and Ch. 14 panel methods** will need it |
| 12 | per-state transport range | `syncRange` ch06 E7 | 1 (1) | engine: `transport.range(s)` callback (quirk "transport range fixed per app") |
| 13 | memoised state / solver run cache | `MEMO` ch06 E9; `cache` ch05 E3, E6; E7 checkpoints of sweep states | 4 (2) | `Viz.memo(fn, keyOf)`; solver runs keyed by parameters with checkpoints every k sweeps for scrubbing |
| 14 | segment stroking and C_p drawing | `strokeSegs` ch06 E2, E6, E8; `drawCp` E1, E2, E8; `drawFlow` E1, E2 | 3 (3) | `P.segments(list, opt)`; C_p stays chapter code |
| 15 | engine fixes proposed by the builders (no helper) | Explain refreshed on pause and scrub (E9 re-renders its Explain by hand); `onChange(s, key, source)` telling user from walkthrough changes (E7 re-applies modes); views in a row need `flex-grow ≥ 1` (a `flex: 0.4` view left a gap) | — | engine: refresh Explain live values on `pause`/`scrub`; pass a `source` ('user', 'tour' or 'preset') to `onChange` (with quirk Q5); clamp row flex ≥ 1 or document it |

**Top five for the next library pass** (re-ranked with ch06): `Viz.cx` complex module (new, needed by Ch. 7, 11, 14),
`arrowPx` (15 files), the y-title gutter as an engine default (7 local copies, three phone Must-fixes since ch01),
formatters (all 44), per-mode hiding (Q4). Then contour `maxJump`/mask, field image, log axes, bars with faint/solid
pairs, Gauss–Legendre, LU/cond1, and the engine fixes of rank 15.

### ch06 machinery TODOs (for the orchestrator; not done in this pass)
- `tools/viz_lint.py`: fail single-backslash TeX in JS string literals — `(?<!\\)\\(f|t|r|b|v|n)(rac|heta|ho|au|eta|u|abla)`
  and, more generally, any control character (`\x08 \x09 \x0b \x0c \x0d`) inside a `'…'`/`` `…` `` string that also
  contains `$` (E2's round-1 Must; the reviewer's scan regex is in `reports/ch06_viz.md`).
- `tools/run_notebook.py`: refuse to run (or unset) `MPLBACKEND=Agg` — with Agg set, inline figures are not captured and
  `coverage_check` fails "no visual output" (orchestrator observation in the ch06 notebook phase). Scripts are the
  opposite case: `scripts/ch06_*.py` call `plt.show()` unless `--no-show`, which blocks a headless run; ch01's scripts
  switch to Agg themselves when no display is present — make that the script template.
- `notebooks/build_ch06.py` / `tools/nbkit.py`: a shared converter from Part F TeX (`\lvert x\rvert`) to prose `|x|` that
  keeps the following space (lesson round 2 N1: ~15 glued "|d|fixed"), strips `*…*` inside `\text{}` (M2), and a scan
  for literal TeX commands outside `$…$` in `coverage_check` (M1: 98 hits passed the checker).
- `pf_sub` patches in `build_ch06.py` (D23 steps 7–8, traps) repeat the ch05 lesson: fix the design (Part F), then drop
  the patch.
- `eq_refs.py` still flags `<meta>` content (E6) — exempt `viz:concept` or keep metas number-free (quirk Q8).
- A `knowledge/` check that tables keep their column count (an unescaped `|` in `|z|` or `"lex"|"book"` breaks a row).

### ch06 lesson candidates for the skills (not promoted in this pass: the brief forbade library and skill edits)
- `interactive-viz` Lessons: (ch06) double every TeX backslash in JS strings — `\f`, `\t`, `\r`, `\b`, `\v` are
  silent control characters, and only a screenshot of the derivation page shows the damage (E2 round-1 Must) · every
  check question's premise and answer are computed with the fluidpy function before shipping (E5 cut-cylinder, E6
  plate) · a walkthrough step's preset must make its claim visible (E4 z² showed nothing) · when a chapter has two sign
  conventions, every status shows both with units at every size · draw tracers and contours inside `P.clip`, and mask
  branch cuts and axis rows in contours · re-read explainer text after every fluidpy change it mirrors (E8 odd N).
- `verify-implementation`: (ch06) **before calling an equation a book slip, compare the whole printed form with the
  reference form times every plausible factor** ((6.82) = r × App. B was wrongly flagged) · a flux sum that telescopes to
  the boundary values is V1, not V4 — use a genuine invariant (discrete cell circulation) · public tests pin a book
  example's geometry and BCs (a planted non-uniform outlet passed until S3) · Newton-based searches are tested over the
  whole slider range (m from 1e-4 to 200) — the default case hid F1 · report singular-corner convergence orders with
  their reason (4/3 for a 270° corner) instead of the smooth-case 2 · symmetric collocation at odd N is exactly singular.
- `math-to-python` §7: (ch06) seed Newton at the problem's own length scale (rings at ℓ·(¼…4) round each singularity,
  ℓ = m/2πU) plus deflation; a double root is found only to √tol · complex inverse maps take the physical root explicitly
  (Zhukhovsky √(z − 2b)√(z + 2b), never numpy's principal √(z² − 4b²)) · put branch cuts inside bodies (`cut_angle`) ·
  stop iterative solvers on a residual far below the wanted accuracy (residual and change ≈ (1 − ρ) × error when
  ρ → 1) · one ρ default per medium across sibling functions (1.2 air / 1000 water), keyword-only when ambiguous ·
  functions following a book with a reversed sign convention take convention-named keywords (`Gamma_cw=`, `Gamma_ccw=`)
  · return components *and* the invariant magnitude in rotated problems (F_perp, F_par, stream angle).
- `teaching-style` Lessons: (ch06) **one numeric case carried through conventions table → tiny example → code →
  explainer** (stagnation −9.16°/−170.84°, L = 24 N/m) · "state what converges and what doesn't" (4/3 with its reason,
  alternating strengths vs converging moments, exact circle panels) · a suspected slip is retracted in writing when it
  was only a normalisation · every figure caption is checked against the plotted levels (M6: "on three teal lines" that
  were not drawn) · derivation *why* texts are checked for direction/sign claims on the whole body (D17 step 8:
  parallel **or antiparallel**) and for limit orders in substitutions (D26) · the ch06 derivation moves in
  `concept_map.md`.
- `colab-notebook` / nbkit: (ch06) execute notebooks without `MPLBACKEND=Agg` (inline figures vanish) · convert design
  TeX to prose with spaces preserved · strip markdown emphasis before `\text{}` · runtime ~170 s alone (ch05 67 s) —
  set timeouts from a loaded run.

### ch05 lesson candidates for the skills (not promoted in this pass: the brief forbade library and skill edits)
- `interactive-viz` Lessons: (ch05) **phone legibility checklist** — all three round-1 Must-fixes were phone legibility at
  390×844: a rotated y title over tick minus signs (third time since ch01: make the strip an engine default), two
  labels drawn over each other, seven term rows in ~120 px; check every view at 360×640 and 390×844 for overlapping
  labels, ≥ ~25 px per term row (else collapse zero rows or hide a schematic row), and only the 0 tick + the value in
  views < 60 px · a derivation's result page is the boxed result + 2–5 key lines + "all steps" (25 → 5 pages) · print
  round-off as "≈ 0 (round-off)", never as a value · no combining diacritics in canvas labels · every "sum so far"
  from one counter · parity tolerances near the achieved agreement · **builders compute, not copy, design claims**
  (three design errors — θ convention, 2a vs 3a, "short sides cancel" — were caught this way) · guard interactive
  singularities (min separation, clamping, halt distance, sub-stepping) before shipping click-to-place.
- `verify-implementation`: (ch05) **probe every piecewise function exactly at its boundaries** (r = a, 0.999a, 1.001a) —
  a central difference straddling the Rankine/cylinder kink gave half the torque and a 1/h "force" (F1) · check a
  velocity field against its own stream function (Hill's exterior u_R sign, F2) · a convention test needs a field where
  the wrong variant differs (a symmetric G hid the planetary-term transpose; add a tilting field) · published closed
  forms reproduced identically are V1 form cross-checks, not V5 (ring speed, García–Haziot) · 40 wrong variants + 5
  re-planted fixes all caught.
- `math-to-python` §7: (ch05) `scipy.special.ellipk/ellipe` take m = k² · when the book's sign slips cancel, code the
  corrected intermediate and keep the printed sign as an option with a test that it reverses the physics ((5.14)) ·
  thin-core or other limited models return `stop_time`/`valid` instead of running silently out of range · user-reachable
  singular inputs (a vortex at the circle centre) are handled, not NaN-propagated · stencil-residual tolerances relative
  to the cancelling terms (|u||ω|/σ), not to the residual.
- `teaching-style` Lessons: (ch05) **keep the term the book drops** in the ★★★ sympy check (D15 carries u_{j,j}(ω_n + 2Ω_n)
  for a generic u) · teach a book slip by computing both versions as curves or numbers (printed −1/(4π) as a dashed
  trace, fluid at G at −1.70 m/s) · every number in prose is the printed number (four Must-fixes; numpy print options
  hid 2.9e-9 as 0 — use format strings) · physical analogies are checked like equations (rotating tank ≠ geostrophy:
  no Coriolis force at rest) · "what would change if" predictions are computed (viscous cellular flow: Γ ∝ e^{−2νt}
  whatever the loop shape) · a named flow (Taylor–Green) or tool (polar vector Laplacian, moment transfer) needs its
  gloss before first use · the ch05 derivation moves in `concept_map.md`.
- `colab-notebook` / nbkit: (ch05) fix the design when a derivation text is wrong; builder `pf_sub` text patches break
  as soon as the design is corrected upstream · notebook runtime depends on machine load (67 s alone, 345 s during
  parallel browser audits) — set execution timeouts from a loaded run · a derivation "check" must describe exactly what
  its cell runs (two ★★★ checks again in round 1).

### ch04 lesson candidates for the skills (not promoted in this pass)
- `interactive-viz` Lessons: (ch04) a sign-convention toggle must flip labels as well as values — pin the label text
  with exact-text rows · every scenario of a multi-scene budget reports its full balance per direction, not just the
  totals · never write explainer JS through a shell heredoc; scan for single-backslash TeX (`\r`, `\t`, `\p` eat
  characters) · walkthrough cards: display-wide equations go to `eq:`, text ≈ 30 words (5–6-page cards at 360×640 in 7
  of 9 ch04 explainers) · a forward-pointer mode is labelled as a preview everywhere it appears · parity rows call the
  explainer's own functions, never literals · per-mode transport ranges and control lists (Q3, Q4) until the engine
  has them.
- `verify-implementation`: (ch04) **test every scenario's full momentum balance (each direction), not only the budget
  totals** — the oblique-jet 50/50 split and the 2 × 2 mean pressure (missing τ₃₃) passed 147 tests and were caught
  only by the derivation review; afterwards four discriminating tests + planted pre-fix variants (26 wrong variants in
  total caught) · a function that accepts 2-D input must either handle the third dimension's physics or raise.
- `math-to-python` §7: (ch04) a truncated series must not stop at the first zero coefficient (symmetric Couette modes
  vanish for even n) and must not alias its projection with a fixed-node quadrature — compute b_n in closed form ·
  plane (2 × 2) stress needs τ₃₃ for p̄ · a scalar where a vector is expected must raise, not become a z-component ·
  two default g values (core 9.80665 vs chapter 9.81) must be passed explicitly when mixed · book typos: code the
  corrected physics and keep a test that proves the printed form wrong ((4.15) "= 0", (4.51) dA, (4.74) gauge sign).
- `teaching-style` Lessons: (ch04) the equations-vs-unknowns ledger threaded through a chapter (0 → 6/13 → 4/5 → 5/5 →
  7/7) · a "four primes / reused symbols" callout up front plus one code name per meaning (`u_rot`, `p_pert`) · a
  derivation check must cite a cell that actually runs (5 of 30 did not in round 1) · explainer IDs ("E1") are jargon
  in the notebook — name the explainer · every reading note of a vector figure is checked against printed numbers (the
  Lamb–Oseen "arrows forward outside" note was false) · the ch04 derivation moves in `concept_map.md`.
- `colab-notebook` / nbkit: (ch04) write builder and design files with the Write tool — a heredoc corrupted the design's
  backslashes · keep one lake/scenario per name across question, tiny example and code.

## Open library bugs (for the next library pass; `assets/viz_lib.js`, found by the ch03 viz-reviewer)
Each is worked around in chapter CSS/JS today; fix in the engine, then drop the workarounds and re-run
`tools/shot.py --chapter ch01/ch02/ch03 --quick`.

| # | Bug | Symptom | Chapter workaround in use | Engine fix |
|---|---|---|---|---|
| 1 | **Pager has no keep-with-next** | a heading, a `W.step` title or a derivation *move* line can end a page alone, its body on the next page (also ch02 stokes Code title, S4) | keep headings short and put the first sentence in the same item (step title carries the equation; E7 D21/D22 short titles; the "now" live line first on a step) | mark headings/move lines `data-keep-next`; `Pager.layout` never ends a page on one |
| 2 | **Page packed at a different height than shown at 360×640** | a page packed for 145 px is displayed at 135 px because layout runs before the final `fit()`, so its last line is clipped or the density changes after packing | drop the preset strip and the mode chips on short/portrait screens in chapter CSS (E1, E4, E5, E6, E7) so the stage height is stable before packing; leave slack on text-heavy pages | re-pack every pager after the final `fit()` (and whenever the stage height changes) |
| 3 | **Tall `\int` glyph top clipped when a slice break falls just above it** | the slice boundary ignores the KaTeX `.vlist` overhang of display-style integrals | E7 writes derivation rows as `\textstyle` integrals inside `aligned` and uses display style only on desktop (`widthClass() === 'wide'`) | measure each row's ink box (including `.vlist` overhang) when choosing a cut |
| 4 | **`optional` controls stay hidden on phones even when their value drives the picture** | E6 real-vortex modes set Γ and σ on log sliders that phones never show, so the reader could not see the values in use | E6's status line leads with the values ("🌀 cyclone Γ 1.8×10⁷ m²/s, σ 36000 m") at every size, pinned by a selftest row; E7 hides other modes' controls instead of marking them optional | an `optional` control whose value differs from its default (or that a preset/mode set) is shown, or its value is echoed as a chip |
| 5 | **Explain stale after a hand-dragged transport parameter** (ch09 R2-M1) — **FIXED in aa34f17** | `app.set` refreshed only `explain.live` when the transport parameter alone changed; six ch09 explainers showed the previous station | none needed now | user-driven `set()` → `updateDynamic(true, true)`; add a slider-drag machinery test |

### Lesson candidates for the skills (not promoted in this pass; one line each for the next pass that may edit skills)
- (ch03 additions, not promoted) `interactive-viz` Lessons: a measured ◇ on the formula's bar/curve is the strongest
  check an explainer can show · a preset that switches to real-world units must move the sliders and echo the values in
  the status · a slider whose meaning depends on the mode gets a neutral name ("Rate k") and a meaning line · `\textstyle`
  integrals in derivation rows until library bug 3 is fixed · the TeX next to an equation number must be that equation
  (retrofit R2) and a heading must not repeat an equation shown below it (R1) · run `tools/eq_refs.py` before review.
- (ch03) `math-to-python` §7: pin every reused letter with a discriminating test (γ = 2S₁₂ vs Γ ≡ S₁₂ vs Γ circulation) ·
  a book typo is exposed by a dimension test and implemented corrected ((3.6), Ex. 3.2) · signed swept volume, never
  |b·n| · extremum in x = r²/σ² returns r = σ√x · substitute generic polynomials before `simplify` when sympy produces
  `Subs(Derivative(...))` · a time-step default is never the space step · trace over axes (0, 1) for (d, d, N) arrays.
- (ch03) `verify-implementation`: form cross-checks against an encyclopedia are V1, not V5; a label says "converged" only
  if an order is asserted · 13 wrong variants patched in, all caught.
- (ch03) `teaching-style` Lessons: explain a measured deviation with a tiny example (the ellipse turns at ≈ γ/4, not
  γ/2) · the observer-dependence thread (steady/unsteady, local/advective, relative/absolute vorticity) · check every
  "Dxx step N" pointer by script after a split (again, round 2) · the ch03 derivation moves in `concept_map.md`.
- (ch03) `colab-notebook` / nbkit: builder strings with LaTeX must be raw strings (coverage check 9 caught a form feed
  from `\frac`) · promote `show_eqs` into nbkit.
- `interactive-viz` Lessons: (ch02) `autoplay: false` + `play: false` on every step that quotes numbers ·
  `\begin{aligned}` for long live lines · scoped `<style data-chapter>` to drop a view row on phones with
  `:not([data-tab="derive"])` · `py:` parity expressions have no builtins · parallel builders use private scratch
  subfolders · a regime-dependent label must give the regime's reason (singular vs quadrature).
- `math-to-python` §7: (ch02) a test asserting f(B) == g(B) between two code functions is not evidence for the book's
  convention; pin to a field with a known answer · build sympy expressions from the library's own `coordinates`
  (assumptions make different symbols) · normalise symbolic vectors with `b/sqrt(b·b)`, not `.norm()` · theorems whose
  hypotheses can fail return a flag instead of asserting · preset names are physics: test them by their defining
  property.
- `verify-implementation`: (ch02) patch every pinned convention with its wrong variant and show a test fails
  (8/8 caught in ch02) · report `hypothesis_ok` rather than asserting through a singularity.
- `teaching-style` Lessons: (ch02) reorder the book when a derivation needs it, and say so ((2.15) before (2.12)) ·
  counter-example residual tables for "is it a vector/tensor?" · assert every symbolic twin against the numeric
  result in the same cell · look at every animation's last frame (clipped readouts appear when numbers grow) · demo a
  Python idiom on an object where a wrong statement fails · the reusable derivation moves listed in
  `knowledge/concept_map.md` → "Reusable derivation moves".
- `colab-notebook` / nbkit: (ch02) `nb.derivation` must not prefix "Eq." to labels like "Exercise 2.8"; design Part E
  reminder rows need a marker so `tools/coverage_check.py` stops raising 28 false name-matching warnings.

## Ideas for later chapters
(seeded from `book.yaml → viz_seeds`; the curator decides)
- **Ch. 2 tensors** (done): the ch01 E5 matrix layout became ch02 E1's clickable direction-cosine matrix; the
  arrow-on-a-plane idea became ch02 E2's element + Mohr stage.
- **Ch. 4 (done — what came of the plans below)**: the CV idea became E1's five-scene budget (face flux bars + measured
  ◇); Cauchy/stress became E3 (G → S, R → τ grid → traction on a rotatable plane, cube mode); the rotating frame became
  E5's two observers + side toggle; Bernoulli became E6's hypothesis decision table with probes along and across
  streamlines; similarity became E9 (prototype/model pair, paired group bars, sphere collapse); Boussinesq became E8
  (rising blob + validity chart). `Viz.Frac` was not promoted (no second user yet).
- **Ch. 5 vorticity (done — what came of the plans above)**: the term-bar idea became E7's budget mode (seven terms,
  zero rows collapsed) and E5's stretching/tilting split; Kelvin became E3 (six flows, ✓/✗ hypothesis table, measured ◇
  that misses for a fixed loop); absolute vs relative vorticity became E7's column and ring modes with the conserved
  quantity flat in amber; new: E1 tubes with a "broken" field, E4 torque ◇ with an R² gap panel, E6 "build the sum",
  E8 click-to-place point vortices, E9 sheet + roll-up. Backup B1 `vortex_rings` (ring toward a wall, leap-frogging)
  not built — reuse E8's stage with `ring_dynamics` (stop at `stop_time`) if Ch. 14 wants it.
- **Ch. 6 potential flow (from ch05)**: E8's wall/bucket images and click-to-place vortices for images of cylinders and
  a cylinder with circulation; E3's hypothesis table for "why the flow stays irrotational"; E6's sum of pieces for
  superposition of elements.
- **Ch. 11 instability (from ch05)**: E9's sheet + roll-up with the linear-theory ghost is the Kelvin–Helmholtz start.
- **Ch. 12 turbulence (from ch05)**: E5's stretching/tilting arrows and Burgers balance for the cascade and fine scales.
- **Ch. 13 GFD (from ch05)**: E7's column over topography ((ζ + f)/h flat), ring moved poleward and the f-table are
  the seeds of PV conservation, Taylor columns and Rossby waves; E4's disc with tilted isopycnals for thermal wind and
  fronts; E2's four-vortex section for gradient-wind (cyclostrophic vs geostrophic) balance.
- **Ch. 14 aerodynamics (from ch05)**: E6's filament stage (segments, square, ring, "build the sum") for the horseshoe
  vortex, lifting line and downwash; E8's wall image for ground effect.
- **Ch. 6 potential flow**: E2's ψ contours at equal Δψ + spacing = speed + draggable gate, E6's B uniform check, the
  ch03 cylinder flow; superposition toggles per element.
- **Ch. 7 waves**: the C14 surface-particle animation as an explainer (backup `kinematic_free_surface` storyboard is in
  `analysis/ch04_curation.md`), E6's unsteady Bernoulli U-tube as the dynamic free-surface condition, E8's
  Boussinesq blob for internal waves.
- **Ch. 8 laminar flows**: E4's exact-solution term balance (already has Couette, Poiseuille, Stokes' first problem) and
  E7's Couette heating; reuse `exact_solution` and `ns_terms_preset`.
- **Ch. 13 GFD**: E5 (two observers, signed Coriolis bars, f-plane highs and lows) + E9's Ro map (a preview now) +
  E8's validity chart for the Boussinesq ocean/atmosphere; needs `Viz.vec3` (rank 14) and a hodograph view for Ekman.
- **Ch. 4 conservation laws (planned before the chapter)**: integral mass/momentum/energy of a moving CV = ch03 E7 (swept band by sign of
  b·n, budget waterfall + measured ◇, three geometries) with a momentum-flux bar per face; Cauchy's equation = ch02 E2
  element + ch03 E4 measured strain rates → Newtonian stress 2μS; rotating frame §4.7 = ch03 E3 observer slider plus a
  Coriolis bar and ch03 E5's co-rotating observer; Bernoulli = E2's probe-vs-float along a streamline.
- **Ch. 3 kinematics** (done — what was planned): streamline/pathline/streakline stage = template field stage + ch01 E2's "ghost + trail" idea;
  for linear flows use ch02 E3's `expm2` closed form for exact pathlines. **Deformation of a fluid element (§3.4)** =
  ch02 E3's "one clock, three squares" with the book's R = G − Gᵀ and ω = ∇×u. Keep the factor-2 readout
  ("spin ½ω₃") and the parity rows "vorticity = vector(R)", "R = 2A". Principal strain rates = ch02 E2 in strain mode.
  Reynolds transport theorem = ch02 E4's box + waterfall bars with a moving boundary.
- **Ch. 4 conservation laws**: dynamic similarity explainer = E5 with presets for Re, Fr, Ro and a model-vs-prototype
  table (promote `Viz.Frac` then); Boussinesq buoyancy reuses E4's column + profile and `brunt_vaisala_sq`.
- **Ch. 4 stress and Cauchy's equation**: ch02 E2 (element + traction split + Mohr), now with the symmetry of τ proved;
  control volumes reuse ch02 E4's face bookkeeping.
- **Ch. 5 circulation and Kelvin's theorem / Ch. 6 lift**: ch02 E5 (curl heat image, paddle wheel, unrolled u·t loop,
  "Stokes failing on purpose" for a line vortex).
- **Ch. 5 vorticity diffusion / Ch. 8 Stokes' first and second problems**: E2's gap + profile + wall-stress stage with
  a transport and a steady (or periodic) ghost; E2's modes show "same PDE, different D".
- **Ch. 7 internal waves**: E4's column + profile + time series on one clock, ω ≤ N marker, atmosphere/ocean modes and
  the two-convention badge where lapse rates appear.
- **Ch. 11 instability (Richardson number)**: E4's badge pattern for Ri < ¼ with the exact criterion text from fluidpy.
- **Ch. 13 GFD**: E4 for stratification (θ view, inversion preset), E5 for Ro/Ri/Rossby-radius scaling.
- **Ch. 15 compressible**: E3's linked state diagrams with term bars for Fanno/Rayleigh lines and shock entropy rise.
- **Ch. 6 ideal flow (done — what came of the plans above)**: superposition became E1's element kits with ψ-bars;
  lift became E2 (term bars of (6.39)) and E5 (Laurent bars, contour-size view); images became E3 with pressure on the
  wall (linked to ch05 `point_vortex_lab` for dynamics); conformal mapping became E6 (two planes, wrong-branch rose
  lines); relaxation became E7; the inverse method E8; added mass E9. Backup `flow_net_sources_vortices` not built.
- **Ch. 7 waves**: reuse E3's sensor trace (tide gauge / bottom-pressure recorder under a passing wave: p decays as
  cosh k(z + H)/cosh kH), E9's speed-part + acceleration-part split (Morison force on a pile: drag + inertia with the
  added-mass coefficient), E5's per-mode bars for Fourier modes of a packet; `Viz.cx` for complex amplitudes.
- **Ch. 9 boundary layers**: E1's element kit as an outer-flow picker (cylinder, half-body, wedge Azⁿ) feeding U_e(x)
  and dp/dx into a Falkner–Skan / Thwaites stage; mark separation where dp/dx turns adverse; the qualitative
  separated-cylinder band of E2 becomes measured data if a dataset is fetched.
- **Ch. 10 CFD**: E7 is the seed — stencil stepping, matrix view with b − Aψ, residual plot with ghosts, predicted vs
  measured ρ, honest refinement order; add a lid-driven cavity ω–ψ mode and multigrid.
- **Ch. 14 aerodynamics**: E6 + E5 → Zhukhovsky airfoil with the Kutta condition (offset circle, trailing-edge
  stagnation), lift vs angle with the Laurent bars; E8's panel mode → vortex panels with the Kutta condition.
- **Ch. 7 waves (done — what came of the plans above)**: the dispersion idea became E1 (tank + c(λ) with ghost limits
  and regime bands + profile); the surface-particle animation became E2 (linear vs exact paths, Stokes drift); E3's
  sensor-trace idea became the pressure profile and sea-bed banner of E1; Fourier-mode bars became E5's chord/tangent
  and pond; new: E4 seiches, E6 in-page rays, E7 jump budgets, E8 two-layer modes, E9 internal-wave beams with a square
  K plane. Morison wave forces (ch06 E9 split) were not built. Backup B1 `linearised_free_surface` (residual scan as an
  explainer) not built — a good small explainer for Ch. 13's rigid-lid approximation.
- **Ch. 11 instability (from ch07)**: E8's interface stage + a shear U₁ − U₂ → Kelvin–Helmholtz growth rate vs k with
  E3's surface-tension term bars as the short-wave cut-off; ρ₁ > ρ₂ as a Rayleigh–Taylor mode (`interface_omega` already
  returns NaN + warning); E9's K plane for Richardson-number criteria.
- **Ch. 13 GFD (from ch07)**: E1's dispersion diagram with f and β sliders (Poincaré ω² = f² + gHk², Kelvin, Rossby
  ω = −βk/(k² + l²) — `group_velocity_numeric` already handles it); E9's square K plane with a second circle for
  inertia–gravity waves (f < ω < N); E6's in-page rays for WKB internal waves in N(z) and topographic Rossby rays; E8's
  two magnifications + "which g′" for reduced-gravity models and the Rossby radius √(g′H)/f; E5's pond step as
  geostrophic adjustment (`linear_evolve` with f); E2's drift panel for Stokes–Coriolis drift.
- **Ch. 15 compressible (from ch07)**: E7's momentum/energy bars and "forbidden" preset → normal shock (an expansion
  shock is the forbidden case); the simple-wave steepening and x–t characteristics of E6's D22 page → Riemann invariants.
- **Ch. 8 laminar (done — what came of the plans above)**: the ch04 term-balance idea became E1's line + parabola = sum
  with a floor tangent at onset; the rotating-cylinder idea stayed the backup (B1 `rotating_cylinders_couette`, not built —
  C04 taught with a slider figure); new: E2 log term bars with a "dropped" band, E3 pad-frame slider with a printed ghost,
  E4 conservative thin-film march onto t^{1/5}, E5 raw/rescaled Stokes-first collapse on a log clock, E6 exponent plane,
  E7 Stokes layer with the rejected root drawn, E8 Stokes/Oseen/ideal on one control set, E9 traction sweep.
- **Ch. 9 boundary layers (from ch08)**: E5's raw/rescaled stage with a log clock → Blasius profiles at several x
  collapsing on η = y√(U/νx) (the temporal layer of E6's vortex-sheet mode as a ghost, C_f√Re_x 1.128 vs 0.664); E6's
  exponent plane + bracket chips → Falkner–Skan m and β; E2's log term bars with a "dropped" band → the boundary-layer
  ordering (δ/L = Re^{−1/2}, ∂p/∂y ≈ 0); E1's floor tangent at onset → separation where the wall shear crosses zero;
  E8's measured crossover circle → inner/outer matching.
- **Ch. 10 CFD (from ch08)**: E5's CN dots on the exact curve + E7's "start from rest" march → scheme stability and order
  panels (CN vs BE vs FTCS on one clock, the CN oscillation after a jump with and without BE start-up); E4's conservative
  finite-volume march with volume in the status.
- **Ch. 11 instability (from ch08)**: E1's channel as the Couette/Poiseuille base state picker; the circular-Couette
  backup as a Taylor–Couette stage with the Rayleigh line Ω₂/Ω₁ = (R₁/R₂)² drawn in the (Ω₁, Ω₂) plane.
- **Ch. 13 GFD (from ch08)**: E7 → the Ekman layer (ω → f, add the hodograph spiral and the transport arrow; the
  "Ekman look-alike" preset already exists); E5's √(νt) clock → spin-up times (E^{−1/2}/f vs L²/ν); E4 → viscous gravity
  currents and lava/ice (μ up to 10¹³ Pa s preset); E9's settling table → droplets, aerosols, sediment (10 µm droplet
  1.21 cm/s).
- **Ch. 9 boundary layers (done — what came of the plans above)**: ch08 E5's raw/rescaled stage became E2 (Blasius collapse
  + a shooting overlay); ch08 E6's exponent plane was not reused — Falkner–Skan became E3's one-dial-through-a-fold with
  equal wall term bars; ch08 E2's log term bars became E1's scaling story (with streamlines lifted by δ*); ch06 E1's outer-flow
  kit became E4's outer-flow picker feeding Thwaites; ch08 E8's measured circle was not needed. New: E5 log-Re dial (cartoon
  + C_p + C_D + regime table), E6 configuration + growth + spectrum, E7/E8 invariant bars (one falls, one stays flat), E9
  force bars that add with height. Backup B1 `ball_swing_magnus` not built.
- **Ch. 10 CFD (from ch09)**: `march_boundary_layer` as an explainer — the von Mises grid with the σ² mapping, BDF2 vs backward
  Euler orders on one log–log panel, the wall-shear error vs Δσ (plate and wedge), and a "march into the fold" preset that
  stops at τ₀ = 0; E2's shooting overlay as the BVP-methods panel (shoot vs `solve_bvp` vs Töpfer).
- **Ch. 11 instability (from ch09)**: **E6's configuration + growth-rate curve + complex-plane spectrum** is the template for
  every normal-mode problem (Kelvin–Helmholtz via a vortex sheet, Rayleigh–Taylor, Orr–Sommerfeld on Blasius and the
  Falkner–Skan profiles of E3 — the inflection point marked by `profile_inflection`, Rayleigh's criterion as a status);
  E7's sech² jet as the Bickley-jet stability case; `Viz.cx` first.
- **Ch. 12 turbulence (from ch09)**: E7/E8's invariant bars and wrong-exponent toggle for turbulent jets and wakes (spreading
  ∝ x with an eddy viscosity; the same J); E1's shape chips with a 1/7-power and log-law profile (H ≈ 1.3 vs 2.59); E4's
  integral march with a turbulent closure; E5's plate drag curve (laminar, transition, turbulent) — cite a primary source
  for 0.074 first.
- **Ch. 13 GFD (from ch09)**: **E9 teacup → the Ekman layer and spin-down** (add f, the hodograph spiral and the pumping
  arrow; the force-bar subtraction becomes pressure − Coriolis); E5's regime dial and E6's street for island wakes (Kármán
  streets in stratocumulus behind islands, with a stratification slider); E7's entrainment arrows for plumes and gravity
  currents.
- **Ch. 14 aerodynamics (from ch09)**: E4 with a panel-method airfoil as the outer flow (laminar separation and stall), E1's
  displacement-thickness lift as the viscous–inviscid correction (N15), E5's form-drag model for bluff vs streamlined bodies.
- **Ch. 10 CFD (done — what came of the plans above)**: ch08's "CN vs BE vs FTCS on one clock" became E2's von Neumann lab (G(θ)
  plane + Noye region + live march; BTCS and CN as schemes) and E3's CFL ring; ch06's stencil stepping with honest order became
  E1 (log–log error vs h with a round-off floor); ch08 E4's conservative march with the volume in the status reappeared as E6's
  "max\|∇·u\| after projection" status and E7's live solver; ch09's `march_boundary_layer` explainer and the shooting-vs-BVP panel
  were **not** built (C09/C14 took the slots). New: E4 discrete-root sign strip, E5 assembly as a transport with a matrix
  inspector, E6 algorithm stages + checkerboard mode, E7 live solver + benchmark + Richardson, E8 counting → LBB. Backup B1
  `operator_splitting_theta` not built.
- **Ch. 11 instability (from ch10)**: E2's complex-plane G(θ) + parameter-plane region + live kick is the discrete twin of
  every dispersion-relation explainer — reuse the three-view layout with σ(k) (growth rate) in place of \|G(θ)\|, the neutral
  curve in the (k, Re) or (k, Ri) plane with a brute-force dot scan, and a kicked perturbation growing at the predicted rate;
  E8's generalised eigen-solver (needs `Viz.num.eigSym`, rank 1) for small Orr–Sommerfeld/Rayleigh matrices; E7's live MAC
  solver for the nonlinear saturation of a shear-layer (KH) or Rayleigh–Bénard instability; ch09 E6's spectrum view for the
  eigenvalues.
- **Ch. 12 turbulence (from ch10)**: E4's modified-equation ghost and D_num bars → numerical vs eddy viscosity (resolution of the
  Kolmogorov scale); E7's verification panel → a turbulent channel compared with DNS profiles (validation, P253); E1's
  log–log slope reading → energy spectra (−5/3).
- **Ch. 13 GFD (from ch10)**: **E6's staggered grid is the Arakawa C-grid** — a shallow-water explainer with η at centres and
  u, v on faces, Coriolis by four-point averaging, a geostrophic-adjustment transport and a Kelvin wave along a wall; E3's x–t
  stencil with √(gH) as the speed (the gravity-wave CFL of a real model: 126 s at 25 km for the deep ocean); E4 → Ekman and
  thermocline layers, upwind tracer advection's numerical diffusivity vs model κ; E7's pump for any live model run.
- **Ch. 15 compressible (from ch10)**: E3's MacCormack/Lax–Wendroff ring → a shock tube (ripples behind the shock = the D18
  dispersion), the x–t stencil with two characteristic families (u ± c), CFL on \|u\| + c.
- **Ch. 11 (done — what came of the plans above)**: the three-view layout became E1/E2 (σ(k) or the (k, ΔU) boundary + the
  eigenvalue plane + the wave); no live eigen-solver was needed (E3 a 3 × 3 complex determinant, E5–E8 parity-checked tables,
  E9 RK4); the MAC saturation run was not built (roll-up shown by the ch05 sheet in animation A2). Backup B1
  `period_doubling_route` not built.
- **Ch. 12 turbulence (from ch11)**: E8's production/dissipation bars → the turbulent kinetic-energy budget across a channel
  (production −⟨uv⟩U′, dissipation, transport) with DNS profiles; E9's two runs on one clock → sensitivity of a turbulent
  signal and "statistics are predictable" (time-mean bars that agree while the states do not); E1's σ(k) → an energy spectrum
  with a cursor; E6's three-way status → "laminar / transitional / turbulent" verdicts that separate *allowed* from *happening*.
- **Ch. 13 GFD (from ch11)**: E7 with a β slider → Rayleigh–Kuo (the criterion becomes "β − U″ changes sign"; the same
  `rayleigh_eigs_contour` table generator with `Upp − β`); E3's parameter plane with a neutral curve and a live zero → the Eady
  growth-rate curve and its short-wave cut-off; E5's ring swap with energy bars → inertial instability (parcel exchange with
  absolute angular momentum, f(f + ζ) < 0); E6 → Ri in mixing schemes and critical levels of internal waves; E4 → thermohaline
  staircases and a (R_ρ, regime) map; E9's canvas projector → the Ekman spiral without three.js; E8's "click a point → mode +
  energy bars" → baroclinic energy conversion.
- **Ch. 15 compressible (from ch11)**: E2's c-plane where two real roots collide → the compressible vortex sheet (stabilised
  above a Mach number); E1's mode switcher → acoustic vs vorticity vs entropy modes.
- **Ch. 12 (done — what came of the plans above)**: ch11 E8's production/dissipation bars became E5's paired budgets with
  the mirrored production term and DNS-era channel numbers; the "turbulent channel against DNS profiles" became E7 (five
  Lee & Moser profiles, draggable fit window) and E8 (closure constants); the −5/3 slope reading became E4 (ladder + model
  spectrum + units inspector); ch01's two-convention badge became E9's two-line status; ch09's raw ↔ rescaled collapse
  became E6 with an invariant line as the judge; ch11 E9's "statistics are predictable" became E1 (two estimators of one
  mean). Not built: numerical vs eddy viscosity (ch10 idea), a MAC 2-D decaying field (it would show the opposite
  cascade), backup B1 `k_epsilon_decay`.
- **Ch. 13 GFD (from ch12)**: **E9's two-convention badge, H-as-transport and ✕-at-the-edge-of-validity are the template
  for every stratified explainer** (Ri-based mixing, stable boundary layer, convective adjustment); **E5's paired budgets
  with a mirrored term → geostrophic and Ekman force balances with bars that sum to zero ("turn off Coriolis" presets) and
  the Lorenz energy cycle (baroclinic conversion g α$\overline{wT'}$ as the shared bar)**; E7's draggable fit window →
  fitting a spectral slope (k⁻³ vs k^{−5/3}) or an Ekman e-folding depth on data; E2's lag-and-spectrum stage → a Rossby
  or inertial peak in a spectrum; E3's click-one-sample-then-thousands → an eddy heat flux $\overline{v'T'}$ from correlated
  v′ and T′; E4's ladder → the two cascades of 2-D / geostrophic turbulence (energy upscale, enstrophy downscale) with the
  Rossby radius marked; E10's log-paced transport and puff ↔ plume modes → spin-up over many rotation periods, tracer
  spreading on isopycnals; E8's "κ turns it, damping slides it" → an Ekman layer with constant K vs K(z) (what each
  closure constant does to the spiral); E6's invariant line → potential vorticity as the invariant that turns rose when a
  wrong scaling is tried; E1's noise · drift · leak error bars → "how long a record for a climate mean?".

## Reference explainers (the depth to match)
See skill `interactive-viz` §4–§5: Shammunul's preferred MIT-mathlet re-implementations (forced damped vibrations — the
explanation panel; amplitude and phase, second order II — a numbered live derivation; angular frequency explorer —
linked views and modes) and the fast.ai labs (FID formula lab, forward noising lab, pixels as parameters, random copy,
overfitting curves, stride & padding playground, 3-D U-Net and ResNet). Passing in-repo templates that use every engine
feature: `templates/viz_example.html`, `templates/viz_example_field.html`, `templates/viz_example_3d.html`. In-chapter
models after ch01: E4 `parcel_stability` (linked views, conventions, 4 derivations) and E3 `heat_work_paths` (term bars,
modes, linked diagrams). After ch02: `cauchy_traction_principal_axes` (the reference-depth Explain panel of the
chapter: 9 sections, a ★★★ derivation with live numbers) and `strain_vs_rotation_split` (the cleanest one-clock story).
After ch03: `reynolds_transport_cv` (★★★ D22 in 12 steps, three geometries, measured ◇, dropped-term panel) and
`spin_and_principal_axes` (5 derivations, co-rotating observer, mode-neutral "Rate k"); best screenshots
`reports/viz/ch03/flow_lines_unsteady/desktop__tour-step1.png`, `reports/viz/ch03/reynolds_transport_cv/desktop__tour-step5.png`.
After ch04 (reviewer "excellent"): `newtonian_stress_lab` (★★★ D09 in 13 steps with colour-switching τ grid, cube mode),
`rotating_frame_coriolis` (5 derivations, ★★★ D15 with the split Coriolis arrow, signed side toggle, 31 rows) and
`which_bernoulli` (hypothesis decision table, 4 derivations, a book slip with live numbers); best screenshot
`reports/viz/ch04/rotating_frame_coriolis/desktop__tour-step5.png` (signed side panel, force side).
After ch05 (reviewer "excellent"): `biot_savart_filament` (four derivations incl. two ★★★ with both book slips shown,
"build the sum" transport, 23 rows), `kelvin_material_loop` (three derivations, six flows, hypothesis table, deliberate
miss for a fixed loop), `baroclinic_torque` (torque ◇ with an R² gap panel, real-case table) and `point_vortex_lab`
(click-to-place with guards, invariants, a caption correction with a number); best screenshot
`reports/viz/ch05/biot_savart_filament/desktop__tour-step4.png` (pieces tip to tail, unrolled integrand, live code).

After ch06 (reviewer "excellent"): `blasius_kutta_contour` (three derivations incl. two ★★★, Laurent bars "size vs what
survives", force flat vs contour radius with a rose "cuts the body" band, tilted-stream mode, 30 rows) — the template
for contour/residue arguments; `laplace_relaxation` (stencil stepping + matrix + residual ghosts + honest 4/3 order) and
`conformal_joukowski` (two planes, same particles, wrong-branch streamlines) — seeds for Ch. 10 and Ch. 14.

After ch07 (reviewer "very good"): `hydraulic_jump` (momentum and energy term bars, a "forbidden" preset, bore and
solitary modes with KdV invariants, 24 rows) and `two_layer_modes` (two ★★★ derivations incl. one checked with sympy,
separate magnifications for surface and interface, both roots with limits, 27 rows); `wave_rays_refraction` (in-page
Hamiltonian rays with parity against `ray_trace`, one Code block per geometry) and `internal_wave_beams` (sign-safe
k < 0, square K plane beside the tank on phones) — seeds for Ch. 13; best screenshots
`reports/viz/ch07/internal_wave_beams/phone__explore.png` (two square views at 360×640) and
`reports/viz/ch07/group_velocity_packets/phone-tall__tour-step7.png` (pond step with the live calm-disc callout).

After ch09 (reviewer "strong" for eight of nine): `wall_jet_invariant` (three ★★★ derivations, the two-invariant bars — one
falls, one stays flat — printed-ODE toggle, gauge preset), `karman_street_stability` (★★★ D14 with live γ, σ, P, Q; spectrum
reaching the axis at 0.2805 — the Ch. 11 stability template) and `thwaites_marching` (four derivations, outer-flow picker,
inspector re-doing (9.50) at a cursor, book/exact criteria side by side); best screenshots
`reports/viz/ch09/wall_jet_invariant/desktop__explain.png` and `reports/viz/ch09/thwaites_marching/desktop__tour-step5.png`
(book criterion first after round 2).
