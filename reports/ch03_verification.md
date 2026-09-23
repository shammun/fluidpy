# Chapter 3 verification — Kinematics                     2026-09-23, commit 266a1ca (working tree: ch03 phase-4 code)

Verifier: math-verifier. Suite: `tests/test_ch03.py` (120 test functions = 120 collected items), figures/metrics:
`tests/ch03_verify_figures.py` → `outputs/ch03/verify/` (git-ignored), cited data: `reference/ch03/` (`make_refs.py`,
`benchmarks.json`, `SOURCES.md`). Nothing in `fluidpy/` or `scripts/` was edited by the verifier.

## Environment
Python 3.11.5 · numpy 2.4.6 · scipy 1.17.1 · sympy 1.14.0 · pint 0.25.3 · matplotlib 3.11.2 (Windows, `.venv`).

Runs: `pytest tests/test_ch03.py -q -p no:cacheprovider` → **120 passed, 0 failed in 28 s** (with the private book file;
without it **116 passed, 4 skipped** — the four V6 tests). Full suite `pytest -q` → **326 passed, 1 failed in 105 s**: the
one failure is `tests/test_machinery.py::test_viz_library_inlined_and_template_lints` (14 "STALE … run
tools/viz_inline.py" lines for ch01/ch02 explainers and templates). Cause: the working tree carries uncommitted edits
to `assets/viz_lib.js`, `assets/viz_base.css` and `tools/shot.py` made outside this phase, so every inlined copy is
stale. It is not a ch03 file and not in the verifier's scope (Open item O5); every ch01 (107) and ch02 (92) physics test
passes. All 10 `scripts/ch03_*.py` exit 0 (2.9–6.4 s each). `tools/check_public.py` → OK (204 files).

Collected items per evidence tag: **V1 54 · V2 32 · V3 14 · V4 4 · V5 5 · V6 4 · V7 7** (= 120).

**Discrimination proofs** (skill lesson ch01): thirteen wrong variants were patched into the live modules in a
subprocess (scratchpad `discriminate.py`; nothing on disk changed) and the whole ch03 suite re-run — every variant fails
at least one test:

| wrong variant patched in | tests that fail |
|---|---|
| R = ½(G − Gᵀ) (vorticity halved — the ch02 factor-2 convention) | 11 |
| Galilean map with the wrong sign, x = x′ − Ut′ | 2 |
| shear rate γ instead of S₁₂ (2 n₁·S·n₂; ch02's Γ ≡ S₁₂ vs ch03's γ = 2S₁₂) | 3 |
| rotating frame ω′ = ω − Ω instead of ω − 2Ω | 2 |
| Leibniz lower-limit term added instead of subtracted | 4 |
| Gaussian r* = σx* (1.2564σ) instead of σ√x* | 3 |
| (3.6) as printed: ∂F/∂s without the factor \|u\| | 1 |
| (3.19) sign slip du_rot = −½ω × dx | 1 |
| unsigned RTT sliver \|b·n\| | 5 |
| (3.23) with ∂u_θ/∂r (r outside the derivative) | 3 |
| Rankine core ω = Γ/2πσ² (ω₀ instead of 2ω₀) | 4 |
| material derivative local − advective | 7 |
| velocity gradient transposed (G[j, i]) | 8 |

## Validation table — A items (CORE, ≥ 2 independent levels, one of V1/V2/V3/V5)
| Concept / Eq. | Tier | fluidpy target | Evidence (V-levels) | Numbers (error, order, residual) | Label | Notes |
|---|---|---|---|---|---|---|
| C01 Lagrangian ↔ Eulerian bridge (3.2); N07, N08 (3.1); D01 | CORE | `lagrangian_map_example`, `lagrangian_velocity_acceleration(_sym)`, `lagrangian_to_eulerian` | V1 worked number X = 2, α = 0.5, t = 1 (3.297 m, 1.649 m/s, 0.824 m/s²), stencils = αx, α²x on 21 (X, t) pairs; V2 stretching, rotation, shear maps: Eulerian u recovered and Du/Dt − d²r/dt² ≡ 0 (sympy); V7 compatibility (3.2) along the rotation map; V3 stencil order; V2 D01 | 1e-14; stencil 1e-8; orders 2.00 / 2.00 | analytic, symbolic, converged | |
| C02 material derivative (3.3)–(3.5); N09–N12; D02 | CORE | `material_derivative(_terms, _sym)`, `streamwise_derivative`, `thermal_front(_terms)`, `expand_indices` | V1 stencils = d/dt F(r(t), t) along an independent `solve_ivp` path line (4th-order difference) on 6 paths; V2 sym = exact on 200 points, vector form = acceleration; V3 order 2.000; V7 steady F → local 0, u = 0 → advective 0, u ⟂ ∇F → 0; V1 climate number (−1e-4 K/s, +0.36 K/h, DT/Dt = 0, "warm advection"), exact split = stencil split on linear and tanh fronts; N12 streamwise = advective term (1e-7), raises at u = 0; V2 dimensions show the printed (3.6) is 1/s, not K/s; N11 index form 3 terms; V2 D02 (chain rule with generic F, r(t); concrete flow with known paths; vector step 9) | 1e-7; order 2.000 | analytic, symbolic, converged | book typo (3.6) tested corrected |
| C03 streamlines (3.7); N14, N15; D03 | CORE | `streamline`, `streamline_slope`, `flux_through_disc`, `cylinder_streamfunction` | V1 Ex. 3.1 line y = x tan ωt′ (cross-product form 1e-12), solid-body circle (1e-10 at rtol 1e-12); V4 ψ constant along computed streamlines (cylinder, 3 seeds; stagnation xy = 1) < 1e-8; V7 stagnation point stops the curve and raises at the seed; V1 slope 2 for u = (1, 2); V4 stream tube: equal flux through two sections (1e-12) and rim streamline stays on R²z = const; V2 D03 | 1e-8–1e-12 | analytic, conserved | |
| C04 path lines (3.8), streak lines, Ex. 3.1; N13, N16–N19; D04, D05 | CORE | `pathline`, `streakline`, `example_3_1`, `unsteady_flow_preset`, `preset_field` | V1 five linear presets + a random 3-D G vs expm (1e-9), Ex. 3.1 path circle (radius 1e-9, parametric form), streak mirror circle (181 releases, 1e-7), closed-form centres, slopes, tangency; V3 error decreases monotonically with rtol (1e-4 → 1e-10); V7 steady flows: streak = path = streamline; presets (rotating strain traceless with axis at Ωt, trochoid for ex31 + current); V2 D04, D05; V6 | radius error 5e-11 (figure) | analytic, symbolic, converged, book-value | |
| C05 Galilean invariance (3.9); N02, N04, N20–N24; D06 | CORE | `galilean_transform`, `acceleration`, `cylinder_flow`, `cylinder_velocity_field`, `frame_acceleration_terms` | V1 a equal in both frames at 200 points of a random unsteady field (1e-7) while the split differs (> 0.1); wrong-sign map and rotating observer rejected; sine-wave worked number (−0.707 + 0.832 = 0.125; wave frame 0 + 0.125); frame terms: local_y = −16/27 = −0.592593, body frame steady, total independent of U_frame (5 values, 1e-8), t ≠ 0 at the same physical point; N22 λ, λ² scaling; cylinder: u·n = 0 on r = a, far field, fluid = body − U, NaN inside; V2 cylinder u, v from ψ, div = curl = 0; V2 D06 (see Open item O1) | 1e-7–1e-14 | analytic, symbolic | |
| C06 relative velocity (3.10); R01–R03 (3.11)–(3.13); D07 | CORE | `velocity_gradient_at`, `relative_velocity`, `strain_rate_tensor`, `rotation_tensor` | V1 stencil G = exact Jacobian (2-D and 3-D, 50 points, 1e-7), index order pinned; V3 stencil order 2.000 and Taylor remainder slope 2.0 (D07); V1 worked number (0.05, −0.02); R = G − Gᵀ (≠ ½(G − Gᵀ)), S + ½R = G, R = −ε·ω | 1e-7 / 1e-15 | analytic, converged | |
| C07 linear strain rate n·S·n; D08 | CORE | `linear_strain_rate`, `measured_strain_rates`, `linear_flow_map` | V1 S₁₁, = n·G·n, rigid rotation 0, worked number e^{0.02} = 1.0202 cm, 0 at 45°; V3 tracked segment (ruler) → n·S·n at order 1.00; V7 rotation of axes invariant (20); V2 D08 (limit, n·A·n = 0) | 1e-13–1e-15 | analytic, converged, symbolic | |
| C08 shear strain rate S₁₂ = ½D(α+β)/Dt; N27; D09, D10 | CORE | `shear_strain_rate`, `rigid_body_velocity`, `material_line_angle` | V1 = ½(G₁₂ + G₂₁), non-perpendicular raises, **γ = 2S₁₂ pinned (S₁₂ = γ/2, not γ)**, worked number dα = atan(0.01), S₁₂ = 0.5; V3 tracked corner closing → 2n₁·S·n₂ order 1.00; V1 rigid motion: S = 0, ω = 2Ω, spin = Ω (20 random); V7 S unchanged by adding a rigid motion, ω changed; V2 D09 (general pair via rotated components), D10 | 1e-9 / 1e-15 | analytic, converged, symbolic | |
| C09 volumetric strain rate (3.14); D11 | CORE | `volumetric_strain_rate`, `material_volume_ratio`, `measured_strain_rates` | V1 trace; Jacobi det e^{Gt} = e^{t tr G} (20 random, 1e-12); worked number e^{0.03}; shear keeps volume; V3 tracked area → tr G order 1.00; V7 50 rotations; V2 D11 (product rule, det(I + G dt) = 1 + tr G dt + O(dt²), trace invariance, Jacobi symbolically); V4 via `material_volume_rate` (C15) | 1e-12–1e-15 | analytic, converged, symbolic | |
| C10 spin = ½ω; R04–R07 (3.15)–(3.18); N28, N29, N34; D12, D13, D14 | CORE | `element_rotation_rate`, `material_line_rotation_rate`, `perpendicular_pair_rotation_rate`, `vorticity(_from_gradient)`, `vorticity_in_rotating_frame`, `parallel_shear_kinematics`, `circulation_circle`, `potential_velocity`, `velocity_potential_2d` | V1 pair average = ½ω₃ for 36 angles × 10 G (1e-14); shear θ̇ = −γ sin²θ (−0.25 at 30°); solid body: lines turn at ω₀, spin ω₀, ω = 2ω₀; V3 tracked pair → ½ω₃ order 1.00; V1 vorticity = exact curl (50 points, 1e-7), ch02 b × x → 2b; V1 ω′ = ω − 2Ω from the numeric curl of u − Ω × x (1e-7), ω − Ω rejected; V1 §3.5 shear kinematics (ω₃ = −γ, AB −γ, BC 0, S, S̄, 45°), aligned element shears / 45° element stretches, both spin −γ/2; circulation 2πr²ω₀ centred and off-centre, 2πB, 0 excluding the axis, Rankine core loop = ω × area; V2 curl ∇φ ≡ 0; V1 φ reconstructed (1e-4, trapezoid), line-vortex path dependence = 2πB (3e-6 rel.); V2 D12, D13, D14; V6 | 1e-7–1e-15 | analytic, converged, symbolic, book-value | |
| C11 (3.19) deformation + rigid rotation; N30, N33; D15 | CORE | `relative_velocity_split` | V1 2-D and 3-D: du = G·dx, strain + rotation = du, rotation ⟂ dx, = ½ω × dx (1e-14), worked number (0.5, 0) + (0.5, 0); V2 D15 for a general symbolic 3 × 3 G (steps 1–6, the ε-swap sign, the −½ω × dx slip rejected) and the index machine (`-1/2 eps_ijk w_k dx_j` = ½ω × dx) | 1e-14 / exact | analytic, symbolic | book cites (3.14), the step uses (3.15) — noted in the docstring |
| C12 principal strain axes; sphere → ellipsoid; N31, N32, N35; D16 | CORE | `principal_strain_rates`, `strain_velocity_principal`, `strain_ellipse_axes`, `deform_circle`, `deform_sphere`, `deform_square` | V1 SC = C diag(λ), orthonormal, det +1, sign convention, (3.21) dū = λ dx̄ (1e-13); shear t = 0.1: exact 1.05125/0.95125 at 43.57°, first order 1.05/0.95, strain only e^{±0.05}; product of semi-axes = det; deformed circle/sphere on the exact ellipse/ellipsoid (1e-12); V3 first-order vs exact axes order 2.01 (pure strain); Fig. 3.14 squares; V2 D16 (sphere → ellipsoid with semi-axes ε(1 + S̄dt), volume 1 + S_ii dt) | 1e-12–1e-14 | analytic, converged, symbolic | |
| C13 polar vorticity (3.23); (3.22), (3.25); N36–N41; D17 | CORE | `polar_vorticity_z(_sym)`, `annular_sector_circulation`, `mean_vorticity_in_disc`, `vortex_velocity_field`, `solid_body_rotation`, `line_vortex` | V1 = Cartesian curl of a non-axisymmetric field (30 points, 1e-7), 2ω₀ / 0 / Rankine / Gaussian through the kind names; V2 sym (3.22) → 2ω₀, (3.25) → 0, non-axisymmetric pair = Cartesian curl; V3 stencil order 2.000; V1 sector: line vortex 0 (1e-12), legs ±Bdθ, Rankine core ω × area, radial legs carry u_r; V3 Γ/area → ω_z order 2.000; δ-core 2B/r² slope −2 (1e-12); Cartesian twins = ch02 fields, numeric curls; V2 D17 (exact four-leg circulation of a general polynomial field → (3.23), steps 4 and 7 separately, leg orientation of the code checked); V6 | 1e-7–1e-14 | analytic, symbolic, converged, book-value | |
| C14 Rankine (3.28), Gaussian (3.29); N42–N44; D18, D19, D20 | CORE | `rankine_vortex`, `gaussian_vortex`, `gaussian_vortex_max_radius`, `vortex_profile` | V1 Rankine peak at σ (1 m/s), continuity, jump in ω, Γ(r ≥ σ) = Γ, core fraction; Gaussian Γ(r) = Γ(1 − e^{−r²/σ²}) (1e-13), no cancellation at r = 1e-9σ (1e-14 vs series), line-vortex tail; r* three routes (brentq, Lambert W, numerical maximisation) 1e-13 / 1e-6, u_max = 0.6382 Γ/2πσ, r* ≠ σx*; V5 SWIRL α = 1.256 (x* inside α's rounding interval; profiles identical when α = x*), Wikipedia Rankine and Lamb–Oseen (σ² = 4νt) forms (1e-12); V2 D18, D19, D20; V6 | see benchmarks | analytic, symbolic, benchmark, book-value | |
| C15 RTT (3.30)–(3.35); N45–N55; D21–D24 | CORE | `reynolds_transport`, `rtt_check`, `volume_integral(_rate_fd)`, `swept_volume_integral`, `swept_terms(_sphere)`, `surface_flux_term`, `volume_rate_term`, `material_volume_rate`, `leibniz_terms/check`, `leibniz_example`, `example_3_2`, `rtt_field`, `rtt_ellipse_2d`, CV shapes | V1 RTT total = FD of the definition (3.31) for a moving/deforming box, sphere, cylinder, cone and 2-D ellipse (rel 1e-7); V2 growing-sphere closed form (1e-13) and the C15 worked number 4.189 + 1.257 = 5.445; V3 FD agreement order 1.998, T4 2.005, sliver error 2.002, Taylor residual 2.0, one-sided quotient 1.008; V1 fixed volume (surface term exactly 0), F = const → F dV*/dt for every shape, signed shrinking sliver; V4 closed surfaces ∮n dA = 0 and ∮x·n dA = dV (1e-12); V1 slab RTT = Leibniz (1e-12); V1 Leibniz cases = sympy rates (1e-10), D21 numbers 18.667 + 128 − 8 = 138.667; V3 Leibniz FD order 2.000; V5 Wikipedia RTT and Leibniz signs; V4 material volume: flux = ∫∇·u (Gauss) and → ∇·u at order 2 in radius, incompressible 0; V1 Ex. 3.2 three routes (1e-12), 0.1047 m³/s, **b·n = 0 on the base while \|b\| > 0**; V2 cone integral (the "[?]" is z); `carried` field: RTT terms cancel (1e-13); V2 D21, D22 ★★★, D23, D24; V6 | see tables | analytic, symbolic, converged, conserved, benchmark, book-value | |

Every computable A item has ≥ 2 independent levels with at least one of V1/V2/V3/V5.

## Validation table — coded B/C items (NOTE, ≥ 1 level) and recaps
| Concept | fluidpy target | Evidence (test) | Label |
|---|---|---|---|
| N02 steady vs unsteady (frame-dependent) | `frame_acceleration_terms` | V1 body frame local = 0, fluid frame local ≠ 0 (`test_frame_terms_V1_…`) | analytic |
| N03 1-D section average | `cross_section_average`, `pipe_profile` | V1 Poiseuille U/2, uniform U, developing profile keeps its mean at 4 stations, u_max → 2ū (`test_section_average_V1_…`) | analytic |
| N04, N20 cylinder in two frames | `cylinder_flow`, `cylinder_streamfunction`, `cylinder_velocity_field` | V1 + V2 (`test_cylinder_flow_V1_…`, `…_V2_…`) | analytic, symbolic |
| N05 coordinate systems | `core.coords` (all 11 callables) | V1 round trips on 10⁴ points (1e-12), θ ∈ [0, π] from +z, φ ∈ (−π, π], axis finite; orthonormal right-handed bases; projections keep \|u\|, inverse; solid body u_φ = ΩR (`test_coordinates_V1_…` ×2) | analytic |
| N08 (3.1) | `lagrangian_velocity_acceleration(_sym)` | V1, V2, V3 (C01 tests) | analytic, symbolic, converged |
| N09 (3.3), N11 (3.5) index form | `material_derivative_sym`, `expand_indices` | V2 D02, `test_material_derivative_V2_index_form` | symbolic |
| N10 thermal advection | `thermal_front(_terms)` | V1 worked number + stencil agreement | analytic |
| N12 (3.6) corrected | `streamwise_derivative` | V1 + V2 dimensions | analytic, symbolic |
| N13 lines coincide in steady flow | `streakline`, `pathline`, `streamline` | V7 | analytic |
| N14 (3.7) slope | `streamline_slope` | V1 | analytic |
| N15 stream tube | `flux_through_disc` | V4 | conserved |
| N17, N18, N19 streak lines, Ex. 3.1, Fig. 3.7 | `streakline`, `example_3_1`, `scripts/ch03_fig3_7_flow_lines.py` | V1, V2, V6; figure verdict below | analytic, symbolic |
| N21–N23 Galilean transformation, nonlinearity, split | `galilean_transform`, `acceleration`, `frame_acceleration_terms` | V1 | analytic |
| N27 rigid motion | `rigid_body_velocity` | V1, V7, V2 D10 | analytic, symbolic |
| N28 rotating frame | `vorticity_in_rotating_frame` | V1, V2 D13 | analytic, symbolic |
| N29 (3.17) potential flow | `potential_velocity`, `velocity_potential_2d` | V2, V1 (path dependence 2πB reported and checked) | symbolic, analytic |
| N31, N32 (3.20)–(3.21) | `principal_strain_rates`, `strain_velocity_principal` | V1 | analytic |
| N34, N35 shear flow, Fig. 3.14 | `parallel_shear_kinematics`, `deform_square` | V1, V6, figure | analytic, book-value |
| N36–N41 (3.22)–(3.27), Γ_ABCD | `core.vortices`, `annular_sector_circulation` | V1, V2, V3 | analytic, symbolic, converged |
| N43 Lamb–Oseen pointer | `gaussian_vortex` | V5 | benchmark |
| N44 max radius | `gaussian_vortex_max_radius` | V1, V5, V2 D20, V6 | analytic, benchmark, book-value |
| N46, N47 Leibniz (3.30), Fig. 3.17 | `leibniz_terms`, `leibniz_check`, `leibniz_example`, `scripts/ch03_fig3_17_leibniz.py` | V1, V3, V5, V2 D21, smoke | analytic, converged, benchmark |
| N48 control volumes | `ControlVolume` shapes, `MovingEllipse2D` | V4 closed-surface identities, exact volume and volume rate | conserved |
| N49–N52 (3.31)–(3.34) | `volume_integral_rate_fd`, `swept_terms(_sphere)`, `swept_volume_integral` | V3 orders | converged |
| N53 interpretations | `material_volume_rate` | V4, V2 D23, D24 | conserved, symbolic |
| N54 Ex. 3.2 | `example_3_2`, `GrowingCone` | V1, V2, V6 | analytic, symbolic, book-value |
| N55 Fig. 3.18 | `rtt_ellipse_2d`, `rtt_field`, `scripts/ch03_fig3_18_rtt.py` | V1 (FD agreement, carried pattern), smoke | analytic |
| R01–R03 (3.11)–(3.13) | `core.tensors` | V1 `test_strain_rotation_split_V1_recap_R_is_G_minus_GT` (the factor-2 convention discriminated) | analytic |
| R04–R06 (3.15)–(3.16) | `vorticity_from_gradient`, `vector_from_antisymmetric` | V1 `test_vorticity_V1_equals_curl_of_field` | analytic |
| R07 (3.18) circulation | `circulation`, `circulation_circle` | V1 Stokes cross-check on the Rankine core | analytic |

Conceptual NOTE items without computable output (N01, N06, N07, N16, N24, N25, N26, N30, N33, N42, N45) are words only
(curation §2); SKIP S01, S02 not coded.

## Derivations (curation §4b, design Part F) — every ★★ and the ★★★ re-derived with sympy
| D | ★ | test | what is re-derived (intermediate lines checked) | finding |
|---|---|---|---|---|
| D01 | ★ | `test_lagrangian_map_V2_derivation` | steps 1–6, invertibility, steady field, forward check Du/Dt = a | ✓ |
| D02 | ★★ | `test_material_derivative_V2_derivation` | chain rule (3.3) with generic F, r(t); (3.5) on a flow with known paths; step 9 vector form | ✓ |
| D03 | ★ | `test_streamline_V2_derivation` | ds = λu ⇒ u × ds = 0, (3.7), \|dx\| = ds, circles for solid body | ✓ |
| D04 | ★★ | `test_pathline_V2_derivation` | steps 2–7 (constants, parametric form, circle), satisfies (3.8), counterclockwise | ✓ |
| D05 | ★★ | `test_streakline_V2_derivation` | steps 1–8 (release-time label, snapshot, circle, one-period fill, both implicit slopes = tan ωt′) | ✓ |
| D06 | ★★ | `test_galilean_V2_derivation` | (3.9) for a generic cubic primed field; step 3 (−U·∇′u′ term; the slip ∂u/∂t = ∂u′/∂t′ rejected); step 4 | **O1**: the design's check_src does not print [0, 0] in sympy 1.14 |
| D07 | ★ | `test_velocity_gradient_V3_…` | Taylor remainder O(\|dx\|²), slope 2.0 | ✓ |
| D08 | ★ | `test_linear_strain_rate_V2_derivation` | limit (ℓ − 1)/dt = n·S·n, n·A·n = 0, step 4 | ✓ |
| D09 | ★★ | `test_shear_strain_rate_V2_derivation` | α, β from the tracked corner, limits G₁₂, G₂₁, ½(α+β)/dt = S₁₂, step 8 in rotated axes = n₁·S·n₂ | ✓ |
| D10 | ★ | `test_rigid_motion_V2_derivation` | S = 0, R = 2G, ω = 2Ω for symbolic U, Ω | ✓ |
| D11 | ★★ | `test_volumetric_strain_rate_V2_derivation` | product rule, det(I + G dt) linear term = tr G, shear exact, trace invariance, Jacobi symbolically | ✓ |
| D12 | ★★ | `test_spin_V2_derivation` | −α + β limit = −R₁₂/2 = R₂₁/2 = ω₃/2 (ε₂₁₃ = −1); α + β gives a different quantity | ✓ |
| D13 | ★★ | `test_rotating_frame_V2_derivation` | curl(u − Ω × x) = curl u − 2Ω for generic u; co-rotating Ω = ω₃/2 | ✓ |
| D14 | ★★ | `test_shear_pair_spin_V2_derivation` | G·e, −γ sin²θ, AB/BC values, partner −γcos²θ, average −γ/2, −1/4 at 30° | ✓ |
| D15 | ★★ | `test_relative_velocity_split_V2_derivation` | steps 1–6 for a general 3 × 3 G, the ε swap, sign slip rejected, index machine | ✓ |
| D16 | ★★ | `test_principal_axes_V2_derivation` | (3.21), moved sphere on the ellipsoid ε(1 + S̄dt), volume 1 + S_ii dt | ✓ (43.6° drift number verified: ½ atan 20 = 43.57°) |
| D17 | ★★ | `test_polar_vorticity_V2_derivation` | exact four-leg circulation of a polynomial (u_r, u_θ) → (3.23); step 4 and step 7 separately; step 10 | ✓ |
| D18 | ★ | `test_rankine_V2_derivation` | ω per piece by (3.23), continuity, Γ, kink maximum | ✓ |
| D19 | ★★ | `test_gaussian_V2_derivation` | Γ(r), u_θ, both limits, (3.23) closes it | **O2**: check number "u_θ(5) = 0.19999998" is wrong |
| D20 | ★★ | `test_gaussian_max_V2_derivation` | f′, ×2x^{3/2}, ×eˣ, trivial root, g(0.5) = 0.351, g(3) = −13.09, nsolve, Lambert W, √x* = 1.12091, f(x*) = 0.6382 | ✓ |
| D21 | ★★ | `test_leibniz_V2_derivation` | Φ, I = Φ(b) − Φ(a), ∂Φ/∂t under the integral, step 6, (3.30) for arbitrary a(t), b(t); wrong lower sign rejected | ✓ |
| D22 | ★★★ | `test_rtt_V2_derivation` | the design's cell (growing sphere, F = t x₁²): (3.35), T4 starts at Δt², lim T3/Δt = surface term; plus (3.32) exact for this F, (3.33) limit, step 8 (tangential motion sweeps nothing); step 12 numerically (`test_rtt_V1_one_dimensional_reduction_is_leibniz`) | ✓ (the cell runs as written) |
| D23 | ★★ | `test_material_volume_V2_derivation` | 2-D material disc, polynomial u: ∮u·n = ∫∇·u, limit = ∇·u(0) | ✓; **O3** on the E7 number |
| D24 | ★★ | `test_material_derivative_from_rtt_V2_derivation` | per-area RTT limit = DF/Dt + F∇·u, ≠ DF/Dt (the hidden term), product rule for ∇·(Fu) | ✓ |

## Functions used by the notebook and explainers (design Part C) — test name each
`test_contract_V1_every_part_c_name_exists_with_signature` checks **every** Part C name (C.0 reused, C.1–C.5 and the
design's "extra callables", 118 entries) exists on `ch03` with the contracted leading parameters, and that the names are
re-exports of the core modules (not copies; `ch02.velocity_gradient_preset is core.kinematics.velocity_gradient_preset`
— the promotion). `test_contract_V1_every_part_c_name_is_exercised_in_this_file` asserts that every contract name is
called as `ch03.<name>` somewhere in the suite (the reused ch02 names are exercised in
`test_contract_V1_reused_names_smoke`, each with a physical assertion: A = ½R, curl of solid body 1.4, `is_irrotational`
true for pure strain, rectangle-loop Γ = ω × area, disc area π, ∯x·n = 4π, …). C.6 scripts:
`test_scripts_V1_drawing_helpers_run` calls `ring_arrows`, `paddle`, `flow_lines_figure`, `shear_elements_frames`
(ABCD keeps its height), `leibniz_strips_figure`, `rtt_blob_figure`; all 10 scripts exit 0.

Contract details the notebook and explainers rely on, each pinned by a test: λ ascending and the +0.5 axis (1, 1)/√2;
semi-axes descending; `material_line_rotation_rate(G, 30°) = −0.25`; `vorticity_from_gradient` always 3 components;
`streakline` column order = release order; `thermal_front_terms` ∂T/∂y = −G; `frame_acceleration_terms(0, 1.5, 1, 1,
0)["local"][1] = −0.592593`; `rtt_field("warming")` = 1 + 0.5x₁ + 0.3t; `rtt_field("carried")` exists (the design's
reported gap is closed); `leibniz_terms.lower` is the subtracted term.

## Convergence studies (scheme | grids | observed order | design order)
| scheme | steps | observed (pairwise) | design |
|---|---|---|---|
| material derivative, central stencils | h = 0.1 … 0.0125 | 2.000 (1.999, 2.000, 2.000) | 2 |
| velocity gradient, central stencils | h = 0.1 … 0.0125 | 2.000 (1.999, 2.000, 2.000) | 2 |
| Lagrangian u, a stencils | h = 0.2 … 0.025 | 2.00 / 2.00 | 2 |
| Taylor remainder of (3.10) | \|dx\| = 0.1 … 0.0125 | 2.0 | 2 |
| polar vorticity stencils (3.23) | h = 0.1 … 0.0125 | 2.0 | 2 |
| sector Γ/area → ω_z | size 0.2 … 0.025 | 2.000 (1.999, 2.000, 2.000) | 2 |
| RTT: FD of (3.31) vs (3.35) | Δt = 0.1 … 0.0125 | 1.998 (1.996, 1.999, 2.000) | 2 |
| (3.32) T4 = ∫_ΔV Δt ∂F/∂t | Δt = 0.08 … 0.01 | 2.005 (2.010, 2.005, 2.002) | 2 |
| (3.34) sliver error | Δt = 0.08 … 0.01 | 2.002 (2.004, 2.002, 2.001) | 2 |
| (3.31) one-sided quotient − total | Δt = 0.08 … 0.01 | 1.008 (1.015, 1.007, 1.004) | 1 |
| Leibniz measured rate | Δt = 0.2 … 0.025 | 2.000 (1.999, 2.000, 2.000) | 2 |
| tracked segment stretch / corner closing / spin / area | dt = 1e-2 … 1.25e-3 | 1.002 (stretch); closing, spin, area within ±0.15 of 1 | 1 |
| strain ellipse, first order vs exact | t = 0.08 … 0.01 | 2.009 | 2 |
| (1/δV)∮u·n dA → ∇·u (material sphere) | R = 0.2 … 0.025 | 2.0 | 2 |
| path line error vs rtol | 1e-4 … 1e-10 | monotone decrease, final 1e-8 | — |

## Conservation / invariant residuals
Closed control surfaces ∮n dA < 1e-12 and ∮x·n dA = d·V (rel 1e-12) for box, sphere, cylinder, cone, ellipse;
Gauss for material volumes: surface flux = ∫∇·u dV (linear 1e-10, nonlinear 1e-8); incompressible rigid motion: flux
< 1e-12; stream tube: equal flux through two sections (rel 1e-12); ψ drift along computed streamlines < 1e-8;
circulation of the line vortex round any enclosing circle = 2πB (1e-10), non-enclosing < 1e-12; Γ_ABCD = 0 (1e-12);
path dependence of φ around the line vortex = 2πB (rel 3e-6); Jacobi det e^{Gt} = e^{t tr G} (1e-12); "carried" field
RTT total 0 (1e-13).

## Benchmarks used (value | our value | source + DOI/URL | date verified)
| value | ours | source | verified |
|---|---|---|---|
| α = r_max²/σ² = 1.256 (3 digits), Lamb–Oseen profile Eq. (10) | x* = 1.2564312086 (inside [1.2555, 1.2565]; their normalisation (1 + 1/2α)(1 − e^{−α}) = 1 ⇔ e^α = 1 + 2α; profiles identical to 1e-12 when α = x*) | Canivete Cuissa & Steiner 2022, SWIRL I, arXiv:2210.05223, Sect. 2.4 Eq. (10) (PDF text read) | 2026-09-23 |
| Rankine v_θ = (Γ/2π) r/a², (Γ/2π)/r; ω = 2Ω, 0 | identical (rel 1e-14) | Wikipedia "Rankine vortex" | 2026-09-23 |
| Lamb–Oseen v_θ = Γ/2πr (1 − e^{−r²/4νt}), ω = Γ/4πνt e^{−r²/4νt} | `gaussian_vortex` with σ² = 4νt identical (1e-12) | Wikipedia "Lamb–Oseen vortex" | 2026-09-23 |
| RTT d/dt∫f dV = ∫∂f/∂t dV + ∮(v_b·n) f dA, n outward; fixed region → derivative inside | sign and fixed-region reduction reproduced | Wikipedia "Reynolds transport theorem" | 2026-09-23 |
| Leibniz: + f(b)b′ − f(a)a′ + ∫∂f/∂x | signs reproduced; "+ lower" slip disagrees with the measured rate by > 0.1 | Wikipedia "Leibniz integral rule" | 2026-09-23 |
Library cross-checks (not literature): `scipy.special.lambertw` (x* to 1e-13), `scipy.linalg.expm`/`numpy.linalg.svd`,
`scipy.integrate.solve_ivp` as the independent route for D/Dt.

## Numbers from the text (book vs ours)  [private values redacted to relative errors]
| book item | our route | relative error | comment |
|---|---|---|---|
| Ex. 3.1: path/streak circle centres and radius at 4 drawing angles | `example_3_1` | < 1e-14 | closed forms identical |
| Ex. 3.1: path line x(t), y(t) and streak line in the release time | numerical `pathline` / `streakline` | < 1e-9 / < 1e-7 | ODE tolerance |
| §3.5 parallel shear flow: ω₃, average spin, S, S̄, principal angle | `parallel_shear_kinematics` (3 values of γ) | 0 (exact) | |
| (3.22)–(3.24) solid body: ω_z, Γ(r) | `polar_vorticity_z`, `circulation_circle` | < 1e-9 / < 1e-14 | |
| (3.25)–(3.27) line vortex: ω_z = 0, Γ = 2πB, 2B/r², Γ_ABCD = 0 | vortices + sector circulation | < 1e-12 | |
| Gaussian vortex radius of maximum speed (6 figures) | `gaussian_vortex_max_radius` | 3.2e-6 (0.0003 %) | book rounding |
| Rankine maximum at σ | grid argmax | < 1e-4 (grid spacing) | |
| Ex. 3.2: dV/dt direct and the RTT intermediate form (sample h, r₀, ṙ) | `example_3_2` | < 1e-12 | book's "b = 0 on the base" replaced by b·n = 0 (tested) |

## Figures reproduced with our code (outputs/ch03/verify/, local)
| figure | file | one-sentence visual verdict |
|---|---|---|
| Fig. 3.7 (Ex. 3.1) | `fig3_7_flow_lines_verify.png` | At ωt′ = 30° the straight streamline and the two unit circles (path centre (−0.5, 0.866), streak centre (0.5, −0.866)) touch at the origin with a common tangent, and the ODE markers sit on the closed forms (radius error 5e-11). |
| Fig. 3.14 | `fig3_14_shear_elements_verify.png` | The aligned square ABCD becomes a parallelogram of unchanged height (slanted sides 1.044 m at t = 0.3 s), while the 45° square PQRS stretches along 45°/compresses along 135° (sides 0.820 / 0.610 m) and both turn clockwise. |
| (3.28)–(3.29) profiles | `vortex_profiles_verify.png` | Rankine u_θ rises linearly to a kink of 1 m/s at r = σ and joins the Γ/2πr tail; the Gaussian peaks at 0.638 m/s at r* = 1.1209σ and merges with the tail by r ≈ 2.5σ; ω_z is a step (2 → 0) vs a smooth Gaussian from 2 s⁻¹. |
| Figs. 3.15–3.16 | `fig3_15_16_elements_verify.png` | The solid-body element keeps its shape and turns 45° per 45° of revolution (pair-average rate +1.000000 rad/s); the line-vortex element is sheared into a sliver, but two perpendicular threads through its centre have pair-average rate −6.4e-3 → −4.0e-4 rad/s as their length shrinks 4× (∝ size², i.e. zero spin, as the book states). |
| convergence summary | `convergence_verify.png` | Ten straight log–log lines with slopes 2.00–2.01 (stencils, sector, RTT, Leibniz, ellipse) and 1.00–1.01 (tracked segments, one-sided quotient), none flattening at round-off. |

## Deviations & justifications
No `# DEVIATION` comment exists in the ch03 code. Book typos implemented as corrected physics and tested: (3.6) carries
|u| ∂F/∂s (`test_streamwise_derivative_V1_…`, `…_V2_dimensions_expose_book_typo`); Ex. 3.2 b·n = 0 on the base while
b ≠ 0 there (`test_example_3_2_V1_…`); the "[?]" in Ex. 3.2's integrand is the factor z (`test_example_3_2_V2_…` — the
integral without z differs); (3.19) follows from (3.15), not (3.14) (docstring of `relative_velocity_split`; D15 test);
Fig. 3.16's C is the B of (3.25) (code uses B). Test-side choices recorded honestly: the solid-body streamline circle is
checked at rtol = 1e-12 (the default rtol = 1e-10 gives 1.6e-9 radius drift over 3 m of arc, consistent with the
tolerance); the line-vortex curl is taken with h = 1e-5 because the h = 1e-4 truncation near r = 0.2 is ∝ h²/r⁴.

## Open items
- **O1 (designer, Part F D06 check_src):** the sympy cell as written ends with `sp.simplify((a − b).doit())` and, under
  sympy 1.14, prints `[-Subs(Derivative(f(xp, yp, t), t), …) + Subs(Derivative(f(-U_1*t + x, -U_2*t + y, _xi_3), _xi_3), _xi_3, t), …]`
  instead of `[0, 0]` — two representations of the same ∂f/∂t′. The mathematics is right; the notebook cell would show a
  confusing non-zero. Fix used in the test: substitute generic cubic polynomials for f and g
  (`repl = {f: Lambda((a, b, c), Σ p_n a^i b^j c^k), g: …}`) and `sp.expand((a − b).subs(repl).doit()) == 0`.
- **O2 (designer, Part F D19 check):** "u_θ(5) = 0.19999998 ≈ 1/5" is wrong in the 8th digit: (1 − e^{−25})/5 =
  0.19999999999722. Write "0.2 (to 11 digits)" or the exact value.
- **O3 (designer, Part F D23 check / E7):** the live number "∮u·n ds = 1.257 = ∫∇·u dA = 0.2 × 6.283" needs an ellipse of
  area 2π (e.g. a = 2, b = 1); the E7 / `rtt_blob_figure` defaults a = 1, b = 0.6 give area 1.885 m² and 0.377 m²/s.
  Confirm which ellipse the explainer uses. (Not a code issue.)
- **O4 (implementer, cosmetic):** docstrings of `linear_strain_rate`, `shear_strain_rate` and `volumetric_strain_rate`
  label the tracked-segment check "V4"; it is a first-order convergence in dt (V3, as `measured_strain_rates` says).
- **O5 (orchestrator):** `tests/test_machinery.py::test_viz_library_inlined_and_template_lints` fails with 14 STALE
  files because `assets/viz_lib.js`, `assets/viz_base.css` (and `tools/shot.py`) have uncommitted edits in the working
  tree; run `tools/viz_inline.py --all` after that library work is finished. Not a ch03 physics item.
- Nothing is labelled `qualitative` or `unverified`.

## Verdict: PASS
All 10 scripts ran; `tests/test_ch03.py` 120/120 pass; every computable CORE row has ≥ 2 independent levels (one of
V1/V2/V3/V5), every coded NOTE row ≥ 1, every Part C function is at least smoke-tested with a physical assertion, and
every ★★/★★★ derivation is re-derived symbolically. The single full-suite failure (O5) is outside chapter 3.
