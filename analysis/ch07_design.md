# Chapter 7 — Gravity Waves: lesson design
(from `analysis/ch07_curation.md` (A 16 · B 193 · C 27 = 236 rows; CORE 16 · NOTE 199 · RECAP 20 · SKIP 1; 37
derivations D01–D37 written out (★ 9 · ★★ 24 · ★★★ 4); E1–E9 + backup B1) and `analysis/ch07.md` (§2b derivations
a-D01…a-D68, §4 implementation rows I1–I51, §5 core module `core/waves.py`, §9 conventions R1–R11 and slips T1–T9);
equations re-read on the rendered pages `chapters/pages/ch07/` (printed = pdf − 27): (7.1)–(7.2) p254, (7.3)–(7.9)
p255, Fig. 7.1 p256, (7.10)–(7.21) p257–p258, (7.22)–(7.27) p259, the unnumbered dynamic line
p259, (7.28)–(7.32) p260, (7.33)–(7.36) p261, (7.37) p262, (7.38)–(7.39) p263, (7.40)–(7.42) p264, (7.43)–(7.45)
p265, (7.46)–(7.48) p266, (7.49)–(7.52) p267, (7.53)–(7.55) p269, (7.56)–(7.59) p270, (7.60)–(7.62) p271, (7.63)–(7.64)
p272, (7.65)–(7.67) p273, (7.68) p274, (7.69)–(7.71) p275, (7.72)–(7.75) p277, (7.76)–(7.79) p278, (7.80)–(7.81) p282,
the jump energy line and (7.82)–(7.83) p283, (7.84a, b) p284, (7.85)–(7.87) p285, (7.88) p286, (7.89)–(7.92) p287,
(7.93)–(7.95) p288, (7.96) p289, (7.97)–(7.102) p290, (7.103)–(7.113) p291, (7.114)–(7.119) p292, (4.9, 4.10) p293,
(7.120)–(7.127) p294, (7.128)–(7.134) p295, (7.135)–(7.138) p296, (7.139)–(7.141) p297, (7.142) p298, (7.143)–(7.146)
p299, (7.147)–(7.152) p302, (7.153)–(7.158) p303, (7.159) p304; 2026-09-24, lesson-designer. The implementer works in parallel
from analysis §4 + curation §8–§9; **Part C is written first and is the contract both sides keep.**)

**Binding conventions for every builder (analysis §9 R1–R11, curation decisions 4–6, §8).**
1. **Imports and aliases.** `from fluidpy import ch07_gravity_waves as ch07`; `from fluidpy.core import waves as W`.
   **`ch07` re-exports every public name of `core/waves.py`** (the ch05/ch06 pattern), so explainer parity rows write
   `ch07.<name>` only (`tools/shot.py` puts only `np`, `math`, `fluidpy`, `core` and `ch07` in the namespace). Units SI:
   x, z, a, H, λ, L [m]; k, l, m, K [rad/m]; ω, N [rad/s]; T [s]; c, c_g, u, w [m/s]; φ, ψ [m²/s]; p, p′ [Pa]; ρ
   [kg/m³]; σ [N/m]; E [J/m²] (surface, interface) or [J/m³] (internal); F [W/m] (surface, per metre of crest) or
   [W/m²] (internal); Q [m²/s] (per unit width); ΔE [J/kg]; angles in radians inside functions, `_deg` at the interface.
2. **Defaults.** `g = G_BOOK = 9.81` m/s² in **both** `core.W` and `ch07` (said in every docstring; CUMULATIVE flags the
   G0/G_BOOK split); ρ = 1000 kg/m³; **H = `np.inf` means deep water** and must work everywhere (closed deep branch, no
   inf·0); surface tension default σ = 0 in the gravity-wave functions and σ = 0.0727 N/m in the capillary ones (the
   notebook passes σ = `ch01.surface_tension_water(293.15)` = 0.07274 N/m and ρ = 998.2 kg/m³ explicitly for C07's
   worked numbers). Overflow-safe hyperbolic ratios everywhere (analysis §9 numerical risks).
3. **Scalar-callable and parity-friendly.** Every public function accepts Python floats and returns a float or a
   `dict` of floats (arrays only when the input is an array). Parity `py:` expressions may use only `ch07.…`, `np.pi`,
   `np.inf`, numbers, strings, lists and dicts (no builtins, no lambdas), so every callable that an explainer mirrors
   takes its dispersion relation as **data** (`H`, `sigma`, `rho`, `g`, a `kind` string), never as a callable; the
   callable-taking general tools (`linear_evolve`, `group_velocity_numeric`, `ray_trace`) also accept `omega_fn=None`
   with the gravity/capillary relation built from those keywords.
4. **Coordinates and signs** (⚠️ callouts where first used): z **up**, still surface z = 0, bottom z = −H; §7.7 two-layer
   problem: origin at the mean free surface, interface at z = −H. Plane ψ with **u = ∂ψ/∂z, w = −∂ψ/∂x** (ch04's with
   y → z). **η is the surface elevation** here but the level-set function in ch04 (4.90) (R1): the surface function is
   f = z − η, and `core.interfaces.kinematic_bc_residual` is always called with f = z − η. p is **gauge**; p′ = p + ρgz
   *(7.30)* in §7.2 but p′ = p − p̄(z) *(7.124)* in §7.8. **ω is a frequency** (a vorticity in ch02–ch06). Code names
   `zeta_particle` (ξ, ζ excursions, §7.2/§7.6/(7.149)) vs `zeta_interface` (§7.7); `theta_phase` (7.72) vs `theta_K`
   (angle of K **above the horizontal**, (7.139)) vs `beam_from_vertical` (= θ_K; R10). `eps2_density` for ε² of
   (7.95). **Internal waves use ∣k∣ and ∇_K ω** (R9): the printed $\omega=kN/K$ (7.138) and
   $\mathbf c_g=\frac{Nm}{K^3}(m\mathbf e_x-k\mathbf e_z)$ (7.145) hold for k > 0 only; Fig. 7.29 draws k < 0.
5. **Reduced gravity has two conventions** (R7): the book's $g'=g(\rho_2-\rho_1)/\rho_2$ (7.117) (`W.reduced_gravity_book(
   rho1, rho2, g, ref="lower")`) vs ch04's `core.similarity.reduced_gravity(rho1, rho2, g)` with ρ₁ below the line
   (`ref="upper"`). Numbers for the callout: ocean ρ₁ = 1025, ρ₂ = 1027 kg/m³ → 0.01910 vs 0.01914 m/s² (0.2 %); oil
   800 over water 1000 → 1.962 vs 2.453 m/s² (25 %).
6. **Energy normalisations switch** (R5): surface and interface E, E_k, E_p are per unit **horizontal area**
   (depth-integrated) and F per metre of crest; internal-wave E is per unit **volume** and F per unit area. Surface E
   uses $\overline{\eta^2}$ (= a²/2 only for a sinusoid).
7. **Book slips taught in corrected form** (each gets "the book prints X; the correct form is Y" where it is used and a
   discriminating wrong-variant option in code, never asserted as physics): **T1** (7.66) prints ½Δω **x** → ½Δω **t**
   (`W.beat_wave(..., printed=True)` gives the frozen envelope; N65, D19 step 4, E5); **T2** "u from (7.28)" before
   (7.44) → (7.27) (N36, D14 step 5); **T3** (7.105) prints $e^{i(kz-\omega t)}$ → $e^{i(kx-\omega t)}$ (`ch07.
   two_layer_sympy()["printed_7105_residual"]` ≠ 0; N118, D29 step 5); **T4** (7.98) ∂φ₁/dz → ∂φ₁/∂z (N111); **T5** §7.1
   "y = 0" → z = 0 (N12); **T6** Ursell remark "(7.88)" → (7.87) (N95); **T8** Exercise 7.6 lake period stays in the
   private JSON (N63); **T9** T = 10 s ↔ 156 m (not 150 m) in deep water (N39).
8. **Colours (text, figures, derivation terms, explainers — one meaning each; curation §5):** free surface `blue` ·
   particles/orbits `teal` · phase (crests, c) `orange` · group/energy (envelope, c_g, F) `accent` purple · pressure
   `amber` · surface tension `rose` · gravity `blue` · lower/denser layer navy fill, upper layer light · deep-water
   limit muted dashed · shallow-water limit muted dotted · residuals and errors black · printed-slip ghosts `rose`
   dashed.
9. **Every book equation is shown in full next to its number** (CLAUDE.md rule 3) — in markdown, derivation steps
   ("substitute (7.28), $\omega^2=gk\tanh kH$"), traps, recaps, notes and every explainer text. Builders reuse
   `show_eqs(text, EQ)` from `notebooks/build_ch06.py` with an `EQ` dict for all 162 labels of ch07 plus the earlier
   labels cited here ((1.5), (1.29), (3.8), (4.9), (4.10), (4.17), (4.56), (4.75), (4.83), (4.86), (4.91), (4.104));
   JS explainers use a local `showEqs`. `tools/eq_refs.py` must list 0 offenders. Exercises are cited as "Exercise 7.19",
   never as bare equations.
10. **nbkit behaviour** (ch04–ch06 lesson): `nb.recap(...)` and `nb.section(...)` end the current CORE block, so every
    RECAP sits **before** the `nb.core` call of the block that uses it (R01–R06 before C02; R07 before C04; R08 before
    C05; R09 before C07; R10 before C11; R11 before C13; R12–R20 before C15). Every `nb.derivation` sits inside its
    curation CORE block with `ref=` a bare equation label ("7.28"). Items taught in another block get a one-line
    `nb.pointer` in their own section (N47, N48, N174, N175 refraction: §7.2 → C10).
11. **Parsing.** Part E is the only part with table rows after its heading; Part F has no line starting with a table
    bar (absolute values are written \lvert … \rvert or in words). Explainer headings are exactly `### E1 ·
    dispersion_relation` … `### E9 · internal_wave_beams` and `### B1 · linearised_free_surface`. Primer terms in Part A
    are the exact Concept text of the Part E rows marked "primer" (coverage_check matches the first 18 characters).
    "Explained by" never names an earlier chapter's C-number (it cites "Ch. 4 §4.10" or a primers.md entry).
12. **Public repo.** No exercise text; the book's ocean numbers (150 m), c_min/λ_m of (7.59), the Exercise 7.6 answer,
    Fig. 7.33's frequency ratio and every exercise answer stay in `tests/book_values_ch07.json`. The notebook shows
    **our** numbers from our inputs; figures are drawn by our code (the St Andrew's cross is our superposition field,
    labelled as an illustration, not the photograph).

Order of parts: C (contract) · A (notebook storyboard) · B (explainers) · D (runtime) · E (prerequisite ledger) · F
(derivations).

---

## Part C — functions the builders will call (the implementer's contract)

**Status column.** **I#** = planned in `analysis/ch07.md` §4 row I# (signature kept or refined here). **§9** = added by
the curation's §9 table. **NEW** = added by this design (in neither) — flagged as the phase asks. **reuse** = exists
(ch01–ch06) and is only called. Docstrings cite § and Eq. (from the rendered page), list symbols with units and
assumptions, and carry the validation label. Return types: `float`/`ndarray` for single quantities, `dict` (keys listed)
for several named results.

### C.0 Reused (existing; called, not changed)
| # | Callable (signature) | Returns | Used by | Status |
|---|---|---|---|---|
| 0.1 | `style.setup_notebook() -> bool` (FAST) · `style.COLORS` · `style.savefig(fig, "ch07", name)` | FAST · palette · path | setup, every figure | reuse |
| 0.2 | `anim.animate(update, frames, fig, interval)` · `anim.show_animation(anim, player="video"\|"frames")` | HTML | C05 (A1), C08 (A3), C09 (A2), C11 (A5), C12 (A4), C16 (A6) | reuse |
| 0.3 | `interact.slider_figure(fn, name, values, *, unit, xlabel, ylabel, title, xrange, yrange, height)` · `interact.animate_figure(frame_fn, times, …)` · `interact.live(fn, **widgets)` | plotly Figure / widget | C04 (IF1, IF3, live), C07 (IF2), C08 (IF7), C09 (IF4), C10 (IF5), C14 (IF6), C15 (IF8) | reuse |
| 0.4 | `embed.show_viz("ch07", slug)` | display | 9 explainer cells | reuse |
| 0.5 | `tools.convergence.observed_order(h, err) -> float` | log–log slope | C02 (residual slope 2), C11 (N94 slope 4), C16 (N157 ε-limit) | reuse |
| 0.6 | `core.interfaces.kinematic_bc_residual(eta, u=None, x=None, t=0.0, h=1e-4, ht=None, **p)` (called with the **level-set** f = z − η, R1) · `surface_normal_speed(eta, x, t, h, ht, **p)` · `laplace_jump_from_balance(sigma, R1, R2, zeta=1e-6, exact=False)` · `capillary_length(sigma, rho, g, rho_other=0.0)` | residual [m/s] · Δp [Pa] · ℓ [m] | R04 (C02), R09 (C07) | reuse |
| 0.7 | `core.bernoulli.unsteady_bernoulli_pressure(dphi_dt, speed, z=0.0, rho=1000.0, g=G0, C=0.0)` | p [Pa] | R06 (C02, speed = 0, g = 9.81 passed) | reuse |
| 0.8 | `core.potential.laplacian_residual(fn, x, y, h=1e-3)` (y ↦ z) | [1/s] | R02 (C02) | reuse |
| 0.9 | `core.streamfunction.velocity_from_streamfunction_2d(psi, x, y, rho=1.0, h=1e-5, **p)` (y ↦ z; u = ∂ψ/∂z, w = −∂ψ/∂x) | (u, w) | N30 (C05), N59 (C08) | reuse |
| 0.10 | `core.kinematics.pathline(u, r0, t0, t_eval, rtol=1e-10, atol=1e-12, …)` · `ch03.example_3_1(t_prime, xi0=1.0, omega=1.0, n=200)` | dict | R08 (C05) | reuse |
| 0.11 | `core.similarity.froude_number(U, l, g=G0)` · `core.similarity.reduced_gravity(rho1, rho2, g=G0)` (ρ₁ denominator) · `internal_froude_number(U, g_prime=None, l=1.0, N=None)` | [–], [m/s²] | R10 (C11), N128 (C14 callout) | reuse |
| 0.12 | `core.stratification.brunt_vaisala_sq(rho0, drho_dz, drho_a_dz, g=G0)` · `parcel_displacement(t, zeta0, N2, w0=0.0)` | N² [1/s²] · ζ(t) [m] | R18, N145 (C15) | reuse |
| 0.13 | `core.navier_stokes.continuity_residual(rho, u, x, t=0.0, h=1e-4, ht=None)` · `buoyancy(rho_pert, rho0, g=G0)` | residual · b [m/s²] | R13, R14 (C15) | reuse |
| 0.14 | `ch01.surface_tension_water(T)` (IAPWS R1-76(2014), V5) · `core.thermo.G_BOOK` | σ [N/m] | C07 | reuse |
| 0.15 | `ch04.linear_wave_fields(a=0.05, k=1.0, g=G, H=np.inf)` · `ch04.linear_wave_surface(x, z, t, a, k, g, H)` · `ch04.bore_speed(h_in, h_out, g=G)` · `ch04.bore_outlet_velocity(h_in, h_out, g)` | callables / dict / [m/s] | N24 parity (C03), N43 limit (C04), N82 moving jump (C11) | reuse |
| 0.16 | `ch05.vortex_sheet_strength(u_above, u_below, convention="ccw")` | γ [m/s] | N107 (C13) | reuse |

### C.1 `fluidpy/core/waves.py` (NEW module, `core.W`; re-exported by `ch07`) — the wave toolkit for Ch. 7, 11, 13, 15
| # | Callable (signature) | Implements (Eq.) | Returns / units | Used by | Status |
|---|---|---|---|---|---|
| 1.1 | `sinusoid(x, t, a, k, omega, direction=+1)` | (7.1)–(7.2) $\eta=a\cos(kx-\omega t)$ ($\eta=a\cos(kx+\omega t)$ (7.61) for direction = −1) | η [m] | C01, C08 | I1 |
| 1.2 | `wave_parameters(*, k=None, lam=None, omega=None, T=None, nu=None, c=None)` — any two consistent inputs (one of k/λ and one of ω/T/ν/c, or λ and c) → all; raises `ValueError` on under/over-determined or inconsistent input | (7.1), (7.4) $c=\omega/k=\lambda\nu$, (7.7) | dict(k, lam, omega, T, nu, c) | C01 | I2 |
| 1.3 | `crest_positions(t, k, omega, n=0)` | (7.3) $kx_{crest}-\omega t=2n\pi$ | x_crest [m] | C01 (D01, from scratch) | I3 |
| 1.4 | `plane_wave(X, K, omega, a, t=0.0)` — X shape (…, 2) or (…, 3), K (2,) or (3,) | (7.5) $\eta=a\cos(\mathbf K\cdot\mathbf x-\omega t)$, (7.6) | η [m] | C01 (Fig. 7.1 remake) | I4 (t added) |
| 1.5 | `phase_velocity_vector(K, omega)` · `trace_velocities(K, omega)` (inf where a component is 0) | (7.8) $\mathbf c=(\omega/K)\mathbf e_K$; trace speeds ω/k, ω/l, ω/m | (2,)/(3,) [m/s] · tuple [m/s] | C01 (N09, N10) | I5 |
| 1.6 | `doppler_frequency(omega, U, K)` | (7.9) $\omega_0=\omega+\mathbf U\cdot\mathbf K$ | ω₀ [rad/s] | C01 (N11) | I6 |
| 1.7 | `omega_gravity(k, H=np.inf, g=G_BOOK)` (k > 0; deep branch √(gk) exactly when H = inf) | (7.28) $\omega=\sqrt{gk\tanh kH}$ | ω [rad/s] | C03, C04, C08, C09, C10, E1, E4, E5, E6 | I7 |
| 1.8 | `omega_capillary_gravity(k, H=np.inf, sigma=0.0727, rho=1000.0, g=G_BOOK)` | (7.56) $\omega=\sqrt{k(g+\sigma k^2/\rho)\tanh kH}$ | ω [rad/s] | C07, N72, E3, E5 | I8 |
| 1.9 | `phase_speed(k, H=np.inf, g=G_BOOK, sigma=0.0, rho=1000.0)` | (7.29) $c=\sqrt{\frac gk\tanh kH}$, (7.57) with σ, (7.45), (7.49), (7.60) with g = 0 | c [m/s] | C04, C07, E1, E3 | I9 |
| 1.10 | `period_from_wavelength(lam, H=np.inf, g=G_BOOK)` · `wavelength_from_period(T, H=np.inf, g=G_BOOK)` · `wavenumber_from_omega(omega, H=np.inf, g=G_BOOK, sigma=0.0, rho=1000.0)` (`brentq` with the physics bracket; closed forms for deep σ = 0; residual asserted < 1e-12 relative) | (7.28) T–λ form $T=\sqrt{\frac{2\pi\lambda}{g}\coth\frac{2\pi H}{\lambda}}$ and its inverse | T [s] · λ [m] · k [rad/m] | C03, C04 (N39), C08, C10, E1, E6 | I10 |
| 1.11 | `fenton_mckee_kh(omega, H, g=G_BOOK)` · `guo_kh(omega, H, g=G_BOOK)` | explicit approximations of (7.28) (V5 cross-check) | kH [–] | C03 (V5 table) | I11 |
| 1.12 | `group_velocity(k, H=np.inf, g=G_BOOK, sigma=0.0, rho=1000.0)` (analytic dω/dk of (7.56)) | (7.69) $c_g=\frac c2\big[1+\frac{2kH}{\sinh 2kH}\big]$, (7.70), c_g = 3c/2 (g = 0) | c_g [m/s] | C06, C09, C10, E1, E3, E5, E6 | I12 |
| 1.13 | `group_velocity_numeric(omega_fn=None, k=1.0, h=None, **disp)` (complex step if ω accepts complex, else 4th-order central; `omega_fn=None` → `omega_capillary_gravity(**disp)`) · `group_velocity_vector(omega_fn, K, h=None)` | (7.67) $c_g=d\omega/dk$; N71 $c_{gi}=\partial\omega/\partial K_i$ | c_g [m/s] · (n,) [m/s] | C09 (from scratch), C16, N71 | I13 |
| 1.14 | `depth_regime(k, H, deep_kH=2.0, shallow_H_over_lambda=0.07)` | book thresholds; errors $1-\sqrt{\tanh kH}$ and $1-\sqrt{\tanh kH/kH}$ | dict(regime ∈ {"deep", "intermediate", "shallow"}, kH, H_over_lambda, deep_error, shallow_error) | C04, E1 | I14 |
| 1.15 | `capillary_minimum(sigma=0.0727, rho=1000.0, g=G_BOOK)` · `min_group_velocity(sigma=0.0727, rho=1000.0, g=G_BOOK)` | (7.58) $c_{min}=(4g\sigma/\rho)^{1/4}$ at $\lambda_m=2\pi\sqrt{\sigma/\rho g}$; c_g,min at $\sigma k^2/\rho g=2/\sqrt3-1$ | dict(c_min, lam_m, k_m) · dict(cg_min, lam, k) | C07, C09 (N72), E3, E5 | I15 (dicts) |
| 1.16 | `beat_wave(x, t, k1, k2, a=1.0, H=np.inf, g=G_BOOK, sigma=0.0, rho=1000.0, omega_fn=None, printed=False)` (`printed=True` = the T1 slip ½Δω x) | (7.66) $\eta=2a\cos(\tfrac12\Delta k\,x-\tfrac12\Delta\omega\,t)\cos(kx-\omega t)$ | dict(eta, envelope, carrier, c, cg_finite, k, omega, dk, domega) | C09 (D19), E5 | I16 (dispersion as data; `printed`) |
| 1.17 | `gaussian_packet(x, t, a, k0, sigma_x, omega_fn=None, order=2, H=np.inf, g=G_BOOK, sigma=0.0, rho=1000.0)` — closed-form packet with ω(k) Taylor-expanded to `order` about k₀ (order 2 = exact chirped Gaussian; order 1 = D20's $a(x-c_gt)$) | (7.68) $\eta=a(x-c_gt)\cos(kx-\omega t)$, D20 | dict(eta, envelope, carrier, c, cg, omega2 = ω″(k₀)) | C09 (D20 check), E5, IF4 | §9 (dispersion as data) |
| 1.18 | `linear_evolve(eta0, x, t, omega_fn=None, eta_t0=None, **disp)` (periodic FFT; each mode with its own ω(k); both directions from η₀ and ∂η/∂t; ω(0) mode handled) · `envelope(eta)` (`scipy.signal.hilbert` modulus) · `linear_evolve_2d(field0, x, z, t, omega_fn_K)` (complex analytic packet; each mode × e^{−iω(k, m)t}; returns the real part) | N02 Fourier superposition; (7.68) numerically; Fig. 7.32 | η (M, N) [m] · envelope · field (M, Nz, Nx) | C01 (N02), C09, C10 (N183), C16 (A6) | I17 |
| 1.19 | `ray_trace(omega_fn=None, x0=(0.0, 0.0), k0=(0.01, 0.0), t_span=(0.0, 100.0), H_fn=None, dim=2, g=G_BOOK, n_out=200, rtol=1e-9)` (`solve_ivp` RK45 on $d\mathbf x/dt=\nabla_{\mathbf k}\omega$, $d\mathbf k/dt=-\nabla_{\mathbf x}\omega$; `omega_fn=None` → gravity ω(∣k∣, H_fn(x))) | (7.79) extended to rays (ours, labelled); D24 | dict(t, x (dim, n), k (dim, n), omega (n,), omega_drift) | C10, E6 | I18 |
| 1.20 | `wave_energy_density(a, rho=1000.0, g=G_BOOK, drho=None)` · `viscous_decay(a0, k, nu, t)` | (7.42) $E=\tfrac12\rho ga^2$ (or ½Δρga², (7.96)); $a_0e^{-2\nu k^2t}$ | E [J/m²] · a [m] | C06, C09, C13, E5 | I19 |
| 1.21 | `interface_omega(k, rho1, rho2, g=G_BOOK)` (ρ₁ > ρ₂ → NaN + `RuntimeWarning`) | (7.95) $\omega=\sqrt{gk\frac{\rho_2-\rho_1}{\rho_2+\rho_1}}$ | ω [rad/s] | C13, E8 | I20 |
| 1.22 | `two_layer_free_surface_omega(k, H, rho1, rho2, g=G_BOOK)` → (ω_bt, ω_bc) · `two_layer_long_wave_speed(H, rho1, rho2, g=G_BOOK)` · `reduced_gravity_book(rho1, rho2, g=G_BOOK, ref="lower")` | (7.111) $\omega^2=gk$, (7.113), (7.116)–(7.117) $c=\sqrt{g'H}$, $g'=g(\rho_2-\rho_1)/\rho_2$ | [rad/s] · [m/s] · [m/s²] | C14, E8 | I21 |
| 1.23 | `internal_wave_omega(k, m, N, l=0.0)` (∣k∣, sign-safe) · `beam_angle(omega, N)` (from the vertical) · `internal_wave_velocities(k, m, N)` (c and $\mathbf c_g=\nabla_{\mathbf K}\omega$ of N∣k∣/K) | (7.137)–(7.139) $\omega=N\cos\theta$; (7.144)–(7.145) | ω [rad/s] · θ [rad] · dict(c (2,), cg (2,), dot) [m/s] | C15, C16, E9 | I22 |
| 1.24 | `real_field(amp, phase)` | N99 $\zeta=\mathrm{Re}\{a\,e^{i(kx-\omega t)}\}$ | real array | C13, C14, C16 | I23 |

### C.2 `fluidpy/ch07_gravity_waves.py` (chapter module; re-exports C.1)
| # | Callable (signature) | Implements (Eq.) | Returns / units | Used by | Status |
|---|---|---|---|---|---|
| 2.1 | `surface_normal(eta_x)` · `surface_velocity(eta_t)` | (7.14) $\mathbf n=(-\eta_x\mathbf e_x+\mathbf e_z)/\sqrt{\eta_x^2+1}$, (7.15) | (2,) [–] · (2,) [m/s] | C02 (D02) | I24 |
| 2.2 | `linear_bernoulli_pressure(dphi_dt, z, rho=1000.0, g=G_BOOK)` | (7.20) $\partial\phi/\partial t+p/\rho+gz\cong0$, (7.30) | dict(p, p_prime) [Pa] | C02, C04 | I25 |
| 2.3 | `surface_wave_sympy()` (cached) | D05–D06: (7.22)–(7.28) | dict(trial, ode, general, B_over_A, A, B, phi, regroup_identity (0), laplace_residual (0), bottom_residual (0), kinematic_residual (0), dynamic_line, dispersion) — sympy | C03 (D05, D06 checks) | I26 |
| 2.4 | `wave_fields(x, z, t, a, k, H=np.inf, g=G_BOOK, rho=1000.0, direction=+1)` | (7.26), (7.27), (7.31), (7.37), (7.46)–(7.48) | dict(eta, phi, psi, u, w, p_prime, omega) (SI) | C03, C04, C05, C08, C12 | I27 |
| 2.5 | `free_surface_residuals(x, t, a, k, H=np.inf, g=G_BOOK, rho=1000.0, sigma=0.0, h=1e-6)` | (7.11), (7.12), (7.16), (7.17), (7.18), (7.19)–(7.21), (7.55) | dict(laplace, bottom, kinematic_exact, kinematic_17, kinematic_linear, dynamic_exact, dynamic_linear) [m/s or m²/s²] | C02, B1 | I28 (`linear` flag dropped: all residuals returned) |
| 2.6 | `free_surface_residual_scan(ka_values, kH=1.0, g=G_BOOK, n_x=64)` | size of the dropped terms vs ka at fixed k: the exact conditions evaluated with the linear solution | dict(ka, kinematic_abs (m/s, ÷ nothing), dynamic_abs (m²/s²), kinematic_rel (÷ aω), dynamic_rel (÷ ga), slope_abs (≈ 2 vs ka), slope_rel (≈ 1 vs ka)) | C02 figure, B1 | §9 (keys fixed here) |
| 2.7 | `pressure_response(k, z, H=np.inf)` | (7.31) ratio $\cosh k(z+H)/\cosh kH$; (7.48) $e^{kz}$; (7.52) → 1 | [–] | C04, E1 | I29 |
| 2.8 | `dispersion_state(lam, H=np.inf, g=G_BOOK, rho=1000.0, a=1.0)` | (7.28), (7.29), (7.45), (7.49), (7.31), (7.69) | dict(k, kH, H_over_lambda, omega, T, c, cg, c_deep, c_shallow, regime, deep_error, shallow_error, p_bottom_fraction, p_surface_amp (ρga, Pa), t_cross_1000km_c_h, t_cross_1000km_cg_h) | E1, C04 live | §9 (keys fixed here) |
| 2.9 | `particle_path(x0, z0, t_eval, a, k, H=np.inf, g=G_BOOK, model="exact")` (`model` ∈ {"exact" (7.33), "linear" (7.34), "taylor1" (7.84)}; DOP853, rtol 1e-10, atol 1e-12) | (7.32)–(7.34), (7.84a, b) | dict(t, x, z) | C05, C12, E2 | I30 |
| 2.10 | `orbit_linear(x0, z0, t, a, k, H=np.inf, g=G_BOOK)` · `orbit_semi_axes(z0, a, k, H=np.inf)` | (7.35a, b), (7.36), (7.46), (7.50) | (ξ, ζ) [m] · dict(A, B, focal_half, sense ("cw")) [m] | C05, E2 | I31 |
| 2.11 | `orbit_state(z0, a, k, H=np.inf, g=G_BOOK, periods=1, model="linear")` | (7.36), (7.85)/(7.86), N92 | dict(A, B, focal_half, B_over_A, omega, T, drift_speed, drift_per_period, eulerian_mean, drift_numeric (NaN unless model = "exact")) | E2 | §9 (keys fixed here) |
| 2.12 | `dyed_line(z0s, t, a, k, H=np.inf, g=G_BOOK, x0=0.0)` (exact path lines, vectorised over z0s) | Fig. 7.22 | dict(x, z) [m] | C12 (A4), E2 | §9 |
| 2.13 | `wave_energy(a, k, H=np.inf, g=G_BOOK, rho=1000.0, method="closed")` (`"quad"` = `dblquad` of (7.38), (7.40); deep truncated at z = −40/k) | (7.38)–(7.42) | dict(Ek, Ep, E) [J/m²] | C06 | I32 |
| 2.14 | `energy_flux(a, k, H=np.inf, g=G_BOOK, rho=1000.0, method="closed")` (`"quad"`: depth–time average of p′u from `wave_fields`) | (7.43)–(7.44), (7.71) | F [W/m] | C06, C09 | I33 |
| 2.15 | `curvature(eta_x, eta_xx, linear=False)` · `capillary_surface_pressure(eta_xx, sigma, eta_x=None, p_a=0.0)` | (7.53), (7.54) $(p)_{z=\eta}=-\sigma\,\partial^2\eta/\partial x^2$ | 1/R [1/m] · p [Pa] | C07 (D15) | I34 |
| 2.16 | `capillary_state(lam, sigma=0.0727, rho=1000.0, H=np.inf, g=G_BOOK)` | (7.56)–(7.58), (7.60) | dict(k, omega, c, cg, gravity_term (g/k, m²/s²), tension_term (σk/ρ, m²/s²), tension_ratio, regime ∈ {"capillary" (ratio > 2), "crossover", "gravity" (ratio < 0.5)}, c_min, lam_m, k_m, cg_min, p_crest_per_a (σk², Pa/m)) | E3 | §9 (keys fixed here) |
| 2.17 | `standing_wave_fields(x, z, t, a, k, H=np.inf, g=G_BOOK)` | $\eta=2a\cos kx\cos\omega t$, (7.62), (7.63) | dict(eta, psi, u, w) | C08, E4 | I35 |
| 2.18 | `seiche_modes(L, H, n=0, g=G_BOOK)` · `basin_modes(L, b, H, m, n, g=G_BOOK)` | (7.64) $\lambda=2L/(n+1)$, (7.65) | dict(k, lam, omega, T) | C08, E4 | I36 |
| 2.19 | `seiche_state(L, H, n=0, g=G_BOOK)` | (7.64)–(7.65), shallow limit $T=2L/((n+1)\sqrt{gH})$ | dict(k, lam, kH, omega, T, T_shallow, shallow_error, T_min) | E4 | §9 (keys fixed here) |
| 2.20 | `packet_state(k0, H=np.inf, g=G_BOOK, sigma=0.0, rho=1000.0, a=1.0, distance=1.0e6, group_length=None)` (group_length default 10λ) | (7.67)–(7.71) | dict(omega, lam, c, cg, ratio, E, F, t_crest_cross (group_length/∣c − c_g∣, s), arrival_h (distance/c_g, h), regime ∈ {"cg<c", "cg=c", "cg>c"}) | E5, C09 worked number | **NEW** |
| 2.21 | `pond_ripples(x, t, width=0.01, n_modes=256, k_max=None, H=np.inf, sigma=0.0727, rho=1000.0, g=G_BOOK)` — direct cosine sum η = Σ A(k) cos kx cos ω(k)t for a released Gaussian hump (A(k) its cosine transform); JS mirrors it exactly | N72, Fig. 7.16 | η [m] (x, t scalar or array) | C09 (N72), E5 stone mode | **NEW** |
| 2.22 | `local_wavenumber_frequency(theta_fn, x, t, h=1e-4)` · `crest_conservation_residual(theta_fn, x, t, H_fn=None)` | (7.73), (7.74), (7.75), (7.79) | (k, ω) · dict(r774, r775, r779) | C10 | I37 |
| 2.23 | `snell_ray_plane_beach(x, alpha0, x0, T, slope, g=G_BOOK)` — plane beach H = slope·x (x offshore distance, shore at x = 0); ray from (x0, 0) at angle α₀ to the shore normal; k(x) by `brentq` from ω fixed; y(x) = ∫ tan α dx | D24 $k\sin\alpha=\text{const}$ with (7.28) | dict(x, y, H, k, alpha, snell) | C10, E6, IF5 | §9 |
| 2.24 | `refraction_state(T, alpha0, H0, H, g=G_BOOK)` — one point of the plane-beach ray at depth H | D24, (7.28), (7.69) | dict(omega, k0, k, alpha, alpha_deg, c, cg, snell, lam) | E6 | **NEW** |
| 2.25 | `nonlinear_wavelet_speed(eta, H, g=G_BOOK, model="simple")` (`"book"` = $\sqrt{gH}+u$ with linear u; `"simple"` = $3\sqrt{g(H+\eta)}-2\sqrt{gH}$) · `simple_wave_evolve(eta0, x, t, H, g=G_BOOK)` | N81 (our Riemann simple wave, labelled) | c [m/s] · dict(x_points (M, N), eta (N,), t_break) | C11 (A5) | I38 |
| 2.26 | `hydraulic_jump(H1, u1=None, Fr1=None, g=G_BOOK)` · `jump_momentum_residual(H1, H2, Q, g=G_BOOK)` | (7.80), (7.81), the jump energy change | dict(H2, ratio, u1, u2, Q, Fr1, Fr2, dE [J/kg], head_loss [m]) · residual [m³/s²] | C11 | I39 |
| 2.27 | `jump_state(H1, Fr1, g=G_BOOK, rho=1000.0, frame="stationary")` | (7.80)–(7.81); moving bore speed $\sqrt{gH_2(H_1+H_2)/2H_1}$ | dict(H2, ratio, u1, u2, Q, Fr2, dE, head_loss, power_loss (ρgQ·head_loss, W/m), mom_in (ρQu₁), mom_out (ρQu₂), p_in (½ρgH₁²), p_out (½ρgH₂²), residual, allowed (bool), bore_speed) [SI, per metre of width] | E7 | §9 (keys fixed here) |
| 2.28 | `stokes_wave_profile(x, t, a, k, g=G_BOOK, order=3)` · `stokes_wave_speed(k, a, g=G_BOOK)` · `STOKES_LIMIT_STEEPNESS` = 0.1410633 (H/λ) | (7.82), (7.83) | η [m] · c [m/s] | C12 | I40 |
| 2.29 | `stokes_expansion_sympy()` (deep water: second-harmonic coefficient from the kinematic condition at O((ka)²), frequency correction from the dynamic condition at O((ka)³)) | N86, N87 (a-D40, stated) | dict(alpha (= 1/2), gamma (= 1)) — sympy | C12 (sympy cell) | **NEW** (curation §4c asks for the cell; a function keeps it tested) |
| 2.30 | `stokes_drift(z0, a, k, H=np.inf, g=G_BOOK)` · `stokes_drift_numeric(z0, a, k, H=np.inf, g=G_BOOK, periods=20)` · `eulerian_mean_u(z, a, k, H=np.inf, g=G_BOOK)` | (7.85), (7.86), N92 | ū_L [m/s] | C12, E2 | I41 |
| 2.31 | `kdv_rhs(eta, x, H, g=G_BOOK)` · `kdv_solve(eta0, x, t_out, H, g=G_BOOK, dt=None)` (IF-RK4 pseudo-spectral, 2/3 dealiasing; cache to `outputs/ch07/kdv_run.npz`) · `kdv_invariants(eta, x)` · `kdv_linear_phase_speed(k, H, g=G_BOOK)` · `ursell_number(a, lam, H)` | (7.87), N94 $c=c_0(1-\tfrac16k^2H^2)$, N95 | arrays · dict(mass, momentum, energy) · c [m/s] · [–] | C11 | I42 |
| 2.32 | `solitary_wave(x, t, a, H, g=G_BOOK)` · `cnoidal_wave(x, t, H, height, m, g=G_BOOK)` (optional) · `kdv_residual_sympy()` | (7.88) $\eta=a\,\mathrm{sech}^2[(\frac{3a}{4H^3})^{1/2}(x-ct)]$, $c=c_0(1+\frac a{2H})$ | η [m] · sympy 0 | C11, E7 | I43 |
| 2.33 | `interface_fields(x, z, t, a, k, rho1, rho2, g=G_BOOK)` · `interface_residuals(x, t, a, k, rho1, rho2, g=G_BOOK)` · `interface_energy(a, k, rho1, rho2, g=G_BOOK, method="closed")` | (7.89)–(7.96), N105, N107 | dict(zeta_interface, phi1, phi2, u1, u2, w1, w2, A, B, gamma_sheet, omega) · dict of residuals of (7.90)–(7.94) · dict(Ek, Ep, E) [J/m²] | C13, E8 | I44 |
| 2.34 | `two_layer_modes(k, H, rho1, rho2, g=G_BOOK, a=1.0, mode="baroclinic", long_wave=False)` · `two_layer_residuals(…)` · `two_layer_sympy()` (cached) · `two_layer_rigid_lid_omega(k, h1, h2, rho1, rho2, g=G_BOOK)` | (7.97)–(7.119), Exercise 7.20 form | dict(omega, A, B, C, b, eta_over_zeta, p_prime_check) · residuals · dict(A, B, C, b, bc_residuals (4 zeros), pressure_residual, common_factor, dispersion_7110, factored, printed_7105_residual) · ω [rad/s] | C14, E8 | I45 |
| 2.35 | `two_layer_state(k, H, rho1, rho2, g=G_BOOK, mode="baroclinic")` | (7.95), (7.110)–(7.117) | dict(omega_bt, omega_bc, omega, c_bt, c_bc, T, eta_over_zeta, g_prime_lower, g_prime_upper, c_long, eps2_density, omega_two_deep) | E8 | §9 (keys fixed here) |
| 2.36 | `boussinesq_linear_sympy()` (cached) | D32–D33: (7.123)–(7.134) | dict(r7126, r7128, r7129, r7130, r7131, r7132, r7133, r7134, w_equation) — each residual 0 | C15 (D32, D33 checks) | I46 |
| 2.37 | `internal_wave_fields(x, z, t, k, m, N, w0, rho0=1000.0, g=G_BOOK, residuals=False)` · `w_equation_residual(w_fn, x, z, t, N, h=1e-3)` · `layered_flow_check(u_fn, v_fn, x, y)` | (7.128)–(7.141), (7.153), (7.142) | dict(u, w, p_prime, rho_prime, zeta_particle, K_dot_u) (+ residuals) · residual · dict | C15, C16 | I47 |
| 2.38 | `internal_wave_energy(k, m, N, w0, rho0=1000.0, g=G_BOOK)` · `internal_energy_budget_residual(x, z, t, k, m, N, w0, rho0=1000.0, g=G_BOOK)` | (7.147), (7.154)–(7.159) | dict(Ek, Ep, E [J/m³], F (2,) [W/m²], cgE (2,)) · residual [W/m³] | C16 | I48 |
| 2.39 | `internal_pe_interface_limit(a, rho1, rho2, g=G_BOOK, rho0=1000.0, eps=1.0)` | (7.150)–(7.152) with a tanh profile of width ε | ⟨∫E_p dz⟩ [J/m²] | C16 (N157) | I49 |
| 2.40 | `internal_wave_state(omega_over_N, N, K=1.0, k_sign=1.0, m_sign=1.0)` (k_sign = −1 gives the book's Fig. 7.29 geometry) | (7.139), (7.144)–(7.146), N150 | dict(theta_K, theta_K_deg, beam_from_vertical_deg, beam_from_horizontal_deg, k, m, omega, cx, cz, cgx, cgz, c, cg, dot, T, T_N) | E9, IF8 | §9 (`direction` refined into k_sign, m_sign) |
| 2.41 | `st_andrews_cross(x, z, t, omega, N, width, amp=1.0)` (four Gaussian-envelope beams with c ⟂ c_g — our illustration, labelled) | N150, Fig. 7.33 | ρ′-like field [–] | C16, E9 | I50 |
| 2.42 | `G_BOOK` (re-exported), `wave_regime_label(kH)` → "deep"/"intermediate"/"shallow" (thin wrapper of `depth_regime`) | book thresholds | str | E1 status parity | **NEW** (tiny) |

### C.3 `scripts/ch07_*.py` (runnable demos, each `--no-show` for headless runs; drawing helpers carry no physics)
`ch07_wave_basics.py` (Fig. 7.1) · `ch07_dispersion.py` (Fig. 7.7, 7.10, c(λ), pressure) · `ch07_orbits.py` (Figs. 7.3–7.5)
· `ch07_energy.py` · `ch07_standing.py` (Figs. 7.11–7.12) · `ch07_groups.py` (Figs. 7.13–7.17) · `ch07_refraction.py`
(Figs. 7.8, 7.9, 7.18) · `ch07_nonlinear.py` (Figs. 7.19, 7.21–7.23) · `ch07_hydraulic_jump.py` (Fig. 7.20) ·
`ch07_interface.py` (Figs. 7.24, 7.27–7.28) · `ch07_internal.py` (Figs. 7.29, 7.31–7.33) — I51. One drawing helper module
**`scripts/ch07_drawings.py`** (NEW): `tank(ax, H, L)` (bottom, still level, z-axis), `wave_surface(ax, x, eta)`,
`orbit_ghosts(ax, x0s, z0s, A, B)`, `cv_box(ax, H1, H2)`, `beam_labels(ax, theta)`, `two_layer_sketch(ax, H)` — the
notebook imports them with `sys.path.insert(0, "scripts")`.

**Flag (functions not in analysis §4):** 2.6, 2.8, 2.11, 2.12, 2.16, 2.19, 2.23, 2.27, 2.35, 2.40 and 1.17 come from
curation §9 (keys fixed above); **2.20 `packet_state`, 2.21 `pond_ripples`, 2.24 `refraction_state`, 2.29
`stokes_expansion_sympy`, 2.42 `wave_regime_label` and `scripts/ch07_drawings.py` are new in this design.** Refinements
of §4 signatures: 1.4 `plane_wave` gains `t`; 1.15 returns dicts; 1.16 takes the dispersion as data and `printed=`;
1.19 keyword defaults; 2.5 drops `linear=` (returns all residuals); 2.40 `direction` → `k_sign`, `m_sign`.
---

## Part A — notebook storyboard (`notebooks/build_ch07.py` → `notebooks/ch07_gravity_waves.ipynb`)

**One line per book section** (cell numbers are estimates for the builder's budget; ≈ 640 cells in all):
- §7.1 → C01 (N01–N11, N167 · D01 · primer P165 · figure, Fig. 7.1 remake) — cells ≈ 8–45
- §7.2 → R01 R02 R03 R04 R05 R06, C02 (N12–N17, N168 · D02 D03 D04 · P166), C03 (N18–N25 · D05 D06 · P167 P168), R07,
  C04 (N26, N37, N38, N39, N42, N43, N46, N173 · D07 D08 D09 · IF1 IF3 · live · **E1**), R08, C05 (N27–N30, N40, N41,
  N44, N45, N169–N171 · D10 D11 · **A1**), C06 (N31–N36, N172 · D12 D13 D14), pointer N47 N48 N174 N175 → C10 — cells ≈ 46–265
- §7.3 → R09, C07 (N49–N57, N176 · D15 D16 · P169 P170 · IF2 · **E3**) — cells ≈ 266–305
- §7.4 → C08 (N58–N63, N177, N178 · D17 D18 · P171 · IF7 · **A3** · **E4**) — cells ≈ 306–345
- §7.5 → C09 (N64–N73, N179–N182 · D19 D20 D21 · P172 P173 · IF4 · **A2** · **E5**), C10 (N47, N48, N74–N80, N174,
  N175, N183, N184 · D22 D23 D24 · P174 P175 · IF5 · **E6**) — cells ≈ 346–445
- §7.6 → R10, C11 (N81–N85, N93–N97, N185, N186, N189 · D25 D26 · **A5** · **E7**), C12 (N86–N92, N187, N188 · D27 ·
  **A4** · **E2**) — cells ≈ 446–520
- §7.7 → R11, C13 (N98–N108, N190–N192 · D28 · P176), C14 (N109–N131, N193, N194 · D29 D30 D31 · IF6 · **E8**) — cells ≈ 521–580
- §7.8 → R12 R13 R14 R15 R16 R17 R18 R19 R20, C15 (N132–N142, N145–N147, N196 · D32 D33 D34 · P177 · IF8), C16 (N143,
  N144, N148–N166, N195, N197–N199 · D35 D36 D37 · P178 · **A6** · **E9**), S01 pointer, summary — cells ≈ 581–645

Every CORE block below follows the order: problem in plain words → idea → primers → maths (notes and derivations, Part F)
→ tiny example → code (fluidpy) + "What does the code above do?" → from-scratch check → visual(s) → notes and "What would
change if…". Code drafts are intent + exact calls; the builder writes every line commented (novice grade, units, the
equation written out next to its number). *expect* gives the numbers the executed cell must print (sanity values computed
for this design with g = 9.81 m/s², ρ = 1000 kg/m³ unless stated). *see / read / change* are the three figure notes.

### A.0 Front matter
1. `nb.title(big_idea=…, roadmap=[…16…], prerequisites=[…])`. **Big idea (draft):** "Drop a stone in a pond, watch swell
   roll onto a beach, or feel a lake slosh after a storm: the water surface wants to be flat, gravity pulls every bump
   back, and the water's inertia overshoots — so the surface oscillates and the oscillation travels. This chapter turns
   that picture into one formula, the dispersion relation $\omega=\sqrt{gk\tanh kH}$ (7.28), and then reads everything
   off it: how fast crests move and why long waves outrun short ones, how the water under a wave goes round in orbits,
   how much energy a wave holds and at which speed that energy travels (the group velocity, not the crest speed), why
   waves turn toward a beach, what happens when they steepen into a hydraulic jump, and how a small density step or a
   continuous stratification carries slow, huge internal waves whose energy leaves at right angles to their crests.
   Ocean swell, tides, tsunamis, thermocline waves and the internal waves of the atmosphere are all here." **Roadmap (one
   line per CORE):** C01 the sinusoidal wave and its vocabulary · C02 the linear free-surface problem · C03 the dispersion
   relation · C04 phase speed, deep and shallow water · C05 particle orbits · C06 wave energy and its flux · C07
   capillary–gravity waves and the slowest ripple · C08 standing waves and seiches · C09 group velocity · C10 rays and
   refraction · C11 the hydraulic jump (and solitons) · C12 Stokes waves and Stokes drift · C13 waves on a density
   interface · C14 barotropic and baroclinic modes, reduced gravity · C15 internal waves: ω = N cos θ · C16 c ⟂ c_g,
   beams and F = c_g E. **Prerequisites:** velocity potential and Laplace's equation (Ch. 6 §6.2) · unsteady Bernoulli
   (Ch. 4 §4.9) · kinematic and dynamic conditions at a moving surface, surface tension (Ch. 4 §4.10, Ch. 1 §1.6) ·
   control-volume momentum (Ch. 4 §4.4) · path lines (Ch. 3 §3.2) · Boussinesq equations and N² (Ch. 4 §4.9, Ch. 1
   §1.10) · complex exponentials (Ch. 1, Ch. 6).
2. `nb.explainer_index([...])` — 9 rows: ("dispersion_relation", "Why do long waves outrun short ones?", "ω(k), c(λ),
   the deep and shallow limits and the pressure under a wave on one clock") · ("particle_orbits", "What does the water
   under a wave actually do?", "orbits shrink with depth, circles become flat ellipses, and at finite amplitude they fail
   to close: Stokes drift") · ("capillary_gravity_waves", "Why is there a slowest ripple?", "two restoring forces,
   gravity for long waves and surface tension for short ones, leave a minimum speed of 23 cm/s") ·
   ("seiche_standing_waves", "Which waves can live in a lake?", "two opposite waves make fixed nodes; the walls pick
   the periods") · ("group_velocity_packets", "Where does a wave's energy go?", "crests move at c, the envelope and the
   energy at c_g = dω/dk") · ("wave_rays_refraction", "Why do waves arrive parallel to the beach?", "frequency is fixed
   along a ray, the wavelength shrinks with depth, the crest swings round") · ("hydraulic_jump", "How high does the water
   jump?", "momentum sets the height, energy sets the direction") · ("two_layer_modes", "What is a baroclinic mode?", "a
   layer over deep water rings in two ways: surface and interface in phase, or in antiphase and slow") ·
   ("internal_wave_beams", "Why does energy leave at right angles to the crests?", "ω = N cos θ, c ⟂ c_g, the St Andrew's
   cross").
3. `nb.setup()`.
4. `nb.code` — **chapter imports** (outside any block): `import numpy as np` · `import matplotlib.pyplot as plt` ·
   `import sympy as sp` · `from fluidpy import ch07_gravity_waves as ch07` · `from fluidpy.core import waves as W` ·
   `from fluidpy import ch01_introduction as ch01, ch03_kinematics as ch03, ch04_conservation_laws as ch04,
   ch05_vorticity_dynamics as ch05` · `from fluidpy.core import interfaces, bernoulli, potential, streamfunction,
   kinematics, similarity, stratification` · `from fluidpy.core.interact import slider_figure, animate_figure, live` ·
   `from fluidpy.core.anim import animate` · `from fluidpy.core.style import savefig` · `import sys; sys.path.insert(0,
   "scripts"); from ch07_drawings import tank, wave_surface, orbit_ghosts, cv_box, beam_labels, two_layer_sketch` ·
   `G = ch07.G_BOOK` (9.81). *explain:* one line per import ("`ch07` re-exports the new wave toolkit `core.waves`, so
   `ch07.omega_gravity` and `W.omega_gravity` are the same function").
5. `nb.md` — **⚠️ Conventions in this chapter** (a two-column table "symbol | meaning here (and before)", then numbers):
   ω = angular **frequency** [rad/s] (vorticity in Ch. 2–6); η(x, t) = surface **elevation** (the level-set function in
   Ch. 4 (4.90): here the surface function is f = z − η); z up, still surface z = 0, bottom z = −H; p gauge, and the
   perturbation $p'\equiv p+\rho gz$ (7.30) in §7.2 but $p=\bar p(z)+p'$ (7.124) in §7.8; ζ = a particle's vertical
   excursion (§7.2, §7.6, §7.8) but the interface displacement (§7.7); θ = the local phase θ(x, t) (7.72) but also K's
   angle above the horizontal (7.139); reduced gravity $g'=g(\rho_2-\rho_1)/\rho_2$ (7.117) (lower density below the
   line) vs Ch. 4's ρ₁ (0.2 % apart in the ocean, 25 % for oil over water); energies per unit **area** for surface and
   interface waves, per unit **volume** for internal waves; g = 9.81 m/s² in every ch07 function (Ch. 1's G0 = 9.80665
   elsewhere). "We compute with the book's conventions and name the other one wherever it differs."
6. `nb.md` — **🔁 Tools from earlier chapters used in this one** (one line each, primer number): partial derivative (P25),
   first-order Taylor (P26), multivariable Taylor (P98), definite integral (P27), solve_ivp (P31/P94), RK4 by hand (P95),
   sympy (P40), sympy expand/series (P117), linear second-order ODE and the e^{λz} trial (P44), Euler's formula (P45),
   chain rule (P49/P91), limits and orders of smallness (P68), level sets and the normal (P75), meshgrid (P76),
   broadcasting (P77), contour/quiver/streamplot (P78), eigenvalues (P80, the idea only), complex conjugate (P81),
   quad/dblquad (P87), parametric curves (P92), expm1 and cancellation (P107), brentq (P108), momentum flux (P114),
   Schwarz's theorem (P121), power of a force (P127), reduced gravity (P131), moving level set (P132), Fourier modes and
   the FFT (P142), 2-D delta (P149), integrals of sines and cosines over a period (P151), the complex plane in numpy
   (P153), sympy for complex algebra (P158), Vieta (P71), animate (P16), slider_figure (P17), show_viz (P18), live
   widgets (P47), assert np.allclose (P15), power laws and log–log plots (P13). Each block repeats the ones it uses in a
   one-line reminder at first use.

### A.1 §7.1 Introduction — C01 (+N01–N11, N167 · D01)
1. `nb.section("7.1", "Introduction", intro="**What is this section about?** Before any fluid mechanics, the words: a
   wave's amplitude, wavelength, period, frequency, wavenumber and phase, the speed at which a crest moves, the
   wavenumber vector of a wave travelling in any direction, and what an observer on a current measures. Everything later
   in the chapter is written in these words.")`
2. `nb.core("C01", "The sinusoidal travelling wave $\\eta=a\\cos[\\frac{2\\pi}{\\lambda}(x-ct)]$ (7.1) and its vocabulary",
   question="What exactly moves at 'the wave speed' — and how do λ, T, k, ω and c fit together?")`
3. `nb.md` — **The problem in plain words:** "Stand on a pier and watch swell arrive: a buoy tied below you rises and
   falls every eight seconds, but it does not travel toward the beach. What travels is the *shape* — the crest. We want
   exact words for that shape and one clean definition of how fast it moves, because every later result (how fast a
   tsunami crosses an ocean, how fast its energy arrives) is a statement about those speeds."
4. `nb.note` — **N01 [C]** "Three families of waves appear in this book: waves on an interface (restored by gravity and
   surface tension — this chapter's first half), internal waves inside a stratified fluid (restored by buoyancy — its
   end) and compression waves (restored by compressibility — Ch. 15). Water waves are neither purely longitudinal nor
   transverse: the water goes round in loops (C05). The whole chapter assumes small amplitude (so everything is linear)
   and waves fast compared with the Earth's rotation (rotation joins in Ch. 13: Poincaré, Kelvin and Rossby waves)."
5. `nb.md` — **The idea** (ASCII + table):
   ```
   snapshot at t:       η(x) = a cos(kx − ωt)      crest where the phase kx − ωt = 0, 2π, 4π, …
   a moment Δt later:   the same cosine, slid right by c·Δt
   phase fixed:  kx − ωt = const  ⇒  x = (ω/k)t + const  ⇒  c = ω/k
   ```
   | symbol | name | unit | from the others |
   |---|---|---|---|
   | a | amplitude (crest height above the still level) | m | — |
   | λ | wavelength (crest to crest) | m | λ = 2π/k |
   | k | wavenumber (radians of phase per metre) | rad/m | k = 2π/λ |
   | T | period (crest to crest at a fixed point) | s | T = λ/c |
   | ν | cyclic frequency | Hz | ν = 1/T |
   | ω | angular frequency (radians of phase per second) | rad/s | ω = 2πν |
   | c | phase speed (speed of a crest) | m/s | c = ω/k = λν |
   "**A wave is a moving phase, not moving water.**"
6. `nb.md` — "> ⚠️ **Common confusion:** in Ch. 3–6 ω was the vorticity (a spin rate, 1/s). From here on ω is the wave's
   angular frequency (rad/s): how many radians of phase pass a fixed point per second."
7. `nb.primer("phase of a wave", "The phase is the argument of the cosine, kx − ωt, measured in radians: 0 at a crest,
   π at a trough, 2π at the next crest. Radians per metre (k) and radians per second (ω) are the natural units; cycles
   (1/λ, ν) differ by 2π. A point that keeps its phase fixed moves to the right when ω/k > 0.", code="import numpy as np
   # numbers\nk, w = 2*np.pi/100, 2*np.pi/8   # λ = 100 m, T = 8 s\nx = np.array([0.0, 25.0, 50.0])   # three
   positions [m]\nprint(k*x - w*0.0, np.cos(k*x))   # phases [0, 1.571, 3.142] rad → cos [1, 0, −1]")` (**P165**).
8. `nb.note` — **N03 [B]** "The same wave in k and ω, $\eta(x,t)=a\cos[kx-\omega t]$ (7.2): with $k=2\pi/\lambda$ and
   $\omega=2\pi/T$ this is (7.1) with $c=\omega/k$. Number: λ = 100 m, T = 8 s → k = 0.0628 rad/m, ω = 0.785 rad/s,
   c = 12.5 m/s." equation `\eta(x,t)=a\cos[kx-\omega t]`, ref "7.2".
9. `nb.derivation("D01", …)` — Part F D01 (5 steps), ref "7.4".
10. `nb.note` — **N04 [B]** "The crest condition $\frac{2\pi}{\lambda}(x_{crest}-ct)=2n\pi=kx_{crest}-\omega t$ (7.3) is D01's
    step 2; the integer n only labels which crest." equation (7.3), ref "7.3". **N05 [B]** "$c=\omega/k=\lambda\nu$ (7.4)
    is D01's result: a crest's Δx/Δt." equation (7.4), ref "7.4".
11. `nb.worked_example("a swell with λ = 100 m and T = 8 s", "1. $k=2\\pi/\\lambda=6.2832/100=0.0628$ rad/m. 2.
    $\\omega=2\\pi/T=6.2832/8=0.785$ rad/s. 3. $\\nu=1/T=0.125$ Hz. 4. $c=\\omega/k=0.785/0.0628=12.5$ m/s — the same as
    $\\lambda\\nu=100\\times0.125$. 5. In 2 s a crest moves $c\\,\\Delta t=25$ m: a quarter wavelength.")`
12. `nb.code` — `p = ch07.wave_parameters(lam=100.0, T=8.0)` · `print(p)` · `t = np.array([0.0, 1.0, 2.0])` ·
    `print(ch07.crest_positions(t, p["k"], p["omega"], n=0))`. *expect:* `{'k': 0.06283, 'lam': 100.0, 'omega': 0.7854,
    'T': 8.0, 'nu': 0.125, 'c': 12.5}`; crest at `[0.0, 12.5, 25.0]` m. *explain:* 1. `wave_parameters` accepts any
    consistent pair and returns all six quantities; 2. `crest_positions` solves (7.3) $kx_{crest}-\omega t=2n\pi$ for
    $x_{crest}$ at each time; 3. the crest advances 12.5 m per second.
13. `nb.check_agree` — **from scratch (curation §7):** `x = np.linspace(-10, 110, 120001)` (1 mm spacing) · for each t in
    `[0, 0.5, 1.0, 1.5, 2.0]`: `eta = W.sinusoid(x, t, 1.0, p["k"], p["omega"])`; keep `x[np.argmax(eta * (x < 50))]`
    (the crest that started at 0) · `slope = np.polyfit(ts, xs, 1)[0]` · `assert np.allclose(slope, p["c"], rtol=1e-3)`.
    Markdown: "Tracking the highest point of the cosine on a fine grid gives 12.500 m/s — the crest really moves at ω/k."
14. `nb.figure` — (7 × 3 in) η(x) at t = 0 (blue) and t = 2 s (blue dashed) for λ = 100 m, a = 1 m, x ∈ [0, 200] m; the
    tracked crest as orange dots with an orange arrow of length cΔt = 25 m; a buoy (teal square) at x = 60 m drawn at its
    two heights (it moved only vertically). Title "The shape moves; the water at the buoy only goes up and down". *see:*
    "two snapshots of the same cosine, shifted right by 25 m; the buoy changed height, not position." *read:* "measure
    the shift of any crest and divide by the time: 25 m / 2 s = 12.5 m/s = ω/k." *change:* "…T were 16 s with the same
    λ: ω halves, the crest moves 12.5 m in 2 s — half as fast (c = λ/T)."
15. `nb.note` — **N02 [B]** "**Why sinusoids.** For linear waves any surface shape is a sum of sinusoids,
    $\eta=\sum_k\hat\eta_ke^{i(kx-\omega(k)t)}+\text{c.c.}$, and each one moves with its own ω(k) (the book leaves the
    general solution to Exercise 7.3; our `W.linear_evolve` does it with the FFT). So one sinusoid is the building block."
    + `nb.code`: two deep-water modes k₁ = 0.1, k₂ = 0.4 rad/m (`x = np.linspace(0, 2000, 4096, endpoint=False)` hmm —
    choose a periodic domain L = 2π/0.02 ≈ 314.16 m with k = 5·0.02 and 20·0.02); `eta = W.linear_evolve(eta0, x,
    [0.0, 10.0], H=np.inf)`; measure each mode's shift by its FFT phase. *expect:* shifts after 10 s: long mode 99.0 m
    (c = 9.90 m/s), short mode 49.5 m (4.95 m/s) — the long one twice as fast (deep water c ∝ 1/√k, C04). *explain:* the
    FFT splits η₀ into modes, `linear_evolve` multiplies each by $e^{-i\omega(k)t}$, the inverse FFT adds them back.
16. `nb.note` — **N06 [B], N07 [B], N08 [B], N09 [B]** "**Waves in any direction.** $\eta=a\cos(kx+ly+mz-\omega t)=a\cos(
    \mathbf K\cdot\mathbf x-\omega t)$ (7.5) with the wavenumber vector $\mathbf K=(k,l,m)$ and $K^2=k^2+l^2+m^2$ (7.6).
    Crests are the planes $\mathbf K\cdot\mathbf x-\omega t=2n\pi$, perpendicular to K and $\lambda=2\pi/K$ (7.7) apart
    along K; they advance along K at $\mathbf c=(\omega/K)\mathbf e_K$, $\mathbf e_K=\mathbf K/K$ (7.8) (met again in C16,
    where c and the energy velocity are perpendicular)." equations (7.5)–(7.8).
17. `nb.code` — `K = np.array([1.0, 1.0])`; `print(np.hypot(*K), 2*np.pi/np.hypot(*K))`; `print(ch07.
    phase_velocity_vector(K, 1.0), ch07.trace_velocities(K, 1.0))`. *expect:* K = 1.414 rad/m, λ = 4.443 m; c = (0.5,
    0.5) m/s (|c| = 0.707 m/s); trace speeds (1.0, 1.0) m/s.
18. `nb.figure` — **Fig. 7.1 remake (N167)** (5 × 5 in): contours η = a (crests, orange) of `ch07.plane_wave` on a 2-D grid
    for K = (1, 1) rad/m, the K arrow (black), the spacing λ = 4.44 m along K and 6.28 m along each axis marked, inset: the
    trace speeds c_x, c_y (grey) vs the phase velocity c (orange). *see:* "straight crest lines at 45°, 4.44 m apart
    along K but 6.28 m apart along x and y." *read:* "along an axis the crests are further apart and pass a fixed point
    at the same rate ω, so they seem to move faster: c_x = ω/k = 1.0 > c = ω/K = 0.707 m/s." *change:* "…K turned to
    (1, 0): crests perpendicular to x, c_x = c, and c_y = ∞ (the crest lines never cross the y-axis)."
19. `nb.md` — **N10 [B]** "> ⚠️ **Common confusion:** the trace speeds $c_x=\omega/k$, $c_y=\omega/l$, $c_z=\omega/m$ are
    each ≥ c = ω/K, so they are **not** the components of c (components would be smaller). Reciprocals do not add like
    components: $c_x\mathbf e_x+c_y\mathbf e_y\ne\mathbf c$ — 1.41 m/s vs 0.71 m/s here."
20. `nb.note` — **N11 [B]** "**A current shifts the frequency you see.** In water moving at uniform U the crests are
    carried along, so a fixed probe sees $\omega_0=\omega+\mathbf U\cdot\mathbf K$ (7.9) (the observed frequency; ω is
    the intrinsic one, measured drifting with the water — one substitution x′ = x − Ut in the phase). A frozen pattern
    (ω = 0) swept past at U gives ω₀ = Uk. From here on every ω is intrinsic." equation (7.9) + `nb.code`: `w = 2*np.pi/8;
    k = w**2/G` (deep swell, C04) · `for U in (-1.0, 0.0, 1.0): w0 = ch07.doppler_frequency(w, np.array([U, 0.0]),
    np.array([k, 0.0])); print(U, w0, 2*np.pi/w0)`. *expect:* U = −1 m/s (against the swell): ω₀ = 0.7225 rad/s, T₀ = 8.70 s;
    U = 0: 8.00 s; U = +1: 0.8483 rad/s, 7.41 s. "Ch. 13 uses exactly this for mountain waves and Rossby waves in a wind."
21. `nb.md` — **What would change if…** "…the medium itself set a rule between ω and k? For water it does — (7.28)
    $\omega=\sqrt{gk\tanh kH}$ — and that single rule makes long waves faster than short ones (C03–C04). The vocabulary
    stays; only the relation changes from wave to wave (capillary C07, interfacial C13, internal C15)."

### A.2 §7.2 Linear Liquid-Surface Gravity Waves — R01–R06, C02 (+N12–N17, N168 · D02 D03 D04), C03 (+N18–N25 · D05 D06), R07, C04 (+N26, N37–N39, N42, N43, N46, N173 · D07 D08 D09 · E1), R08, C05 (+N27–N30, N40, N41, N44, N45, N169–N171 · D10 D11), C06 (+N31–N36, N172 · D12 D13 D14), pointer N47 N48 N174 N175
1. `nb.section("7.2", "Linear Liquid-Surface Gravity Waves", intro="**What is this section about?** The heart of the
   chapter. We set up the flow under small waves on a layer of depth H — Laplace's equation with conditions at the
   bottom and at the moving surface — make it linear, solve it, and read off the dispersion relation, the speed of
   crests in deep and shallow water, the pressure under a wave, the orbits of the water and the energy a wave carries.")`
2. `nb.recap("R01", "The velocity potential", "Where the flow has no vorticity the velocity is the gradient of one
   scalar, $u=\\partial\\phi/\\partial x$, $w=\\partial\\phi/\\partial z$ (7.10). Waves started from rest by gravity stay
   irrotational (Kelvin), so this is exact here.", where="Ch. 6 §6.2")`
3. `nb.recap("R02", "Laplace's equation", "Incompressibility $\\partial u/\\partial x+\\partial w/\\partial z=0$ with (7.10)
   gives $\\partial^2\\phi/\\partial x^2+\\partial^2\\phi/\\partial z^2=0$ (7.11) — the same Laplacian as Ch. 6 with y
   renamed z.", where="Ch. 6 §6.2")` + `nb.code`: `phi = lambda x, z: np.exp(z) * np.sin(x)` (a deep-water-shaped potential,
   k = 1) · `print(potential.laplacian_residual(phi, 0.3, -0.5))` (the ch06 stencil with y ↦ z). *expect:* ≈ 1e-9 (stencil
   round-off). "(Ch. 6 P-reminder: a harmonic function has zero Laplacian.)"
4. `nb.recap("R03", "No flow through the bottom", "A solid flat bottom at z = −H lets no water through: $w=\\partial
   \\phi/\\partial z=0$ on $z=-H$ (7.12) — the wall-is-a-streamline condition of Ch. 6 (6.16).", where="Ch. 6 §6.2")`
5. `nb.recap("R04", "The kinematic condition at a moving surface", "Water on the surface stays on it: the fluid velocity
   normal to the surface equals the surface's own normal speed, $(\\mathbf n\\cdot\\mathbf u)_{z=\\eta}=\\mathbf n\\cdot
   \\mathbf U_s$ (7.13) — Ch. 4 wrote it as $D(\\text{surface function})/Dt=0$ (4.91). ⚠️ In Ch. 4 η was that surface
   *function*; here η(x, t) is the *elevation* and the function is f = z − η.", where="Ch. 4 §4.10")` + `nb.code`:
   `f = ch04.linear_wave_fields(a=0.05, k=1.0, H=2.0)` (the Ch. 4 test field — this chapter derives it in C03) · the
   residual of `interfaces.kinematic_bc_residual` called with the level-set function z − η (as `ch04.
   wave_kinematic_residual` does) at x = 0.3, t = 0. *expect:* nonzero, of order a²kω ≈ 8×10⁻³ m/s: the *exact* condition is not met by
   the linear solution — the gap C02 is about.
6. `nb.recap("R05", "The free surface feels only the air", "With no surface tension and no shear, stress continuity
   (Ch. 4 §4.10) says the water pressure just below the surface equals the air's: $(p)_{z=\\eta}=0$ (7.19) in **gauge**
   pressure (atmospheric = 0).", where="Ch. 4 §4.10")`
7. `nb.recap("R06", "Unsteady Bernoulli", "For irrotational flow $\\partial\\phi/\\partial t+\\tfrac12\\lvert\\nabla\\phi
   \\rvert^2+p/\\rho+gz=C(t)$ (4.83); dropping the small square gives $\\partial\\phi/\\partial t+p/\\rho+gz\\cong0$ (7.20)
   (the constant is zero on the still surface far away).", where="Ch. 4 §4.9")` + `nb.code`: `print(bernoulli.
   unsteady_bernoulli_pressure(dphi_dt=-0.5, speed=0.0, z=-1.0, rho=1000.0, g=9.81))` and `print(ch07.
   linear_bernoulli_pressure(-0.5, -1.0))`. *expect:* both p = 10 310 Pa (= −ρ(φ_t + gz) = −1000(−0.5 − 9.81)); `ch07`
   also returns p′ = p + ρgz = 500 Pa.
8. `nb.core("C02", "The linear free-surface problem: move the conditions to z = 0 — $(\\partial\\phi/\\partial z)_{z=0}
   \\cong\\partial\\eta/\\partial t$ (7.18)", question="The surface is itself an unknown that moves. How can we impose
   conditions on it — and what do we give up?")`
9. `nb.md` — **The problem in plain words:** "A leaf floating on swell rides up and down with the surface and never
   leaves it; the air above pushes on the water with the same pressure everywhere. Those are our two surface conditions.
   The trouble: they hold on the surface z = η(x, t), which is part of the answer. For small, gentle waves we can
   apply them on the flat level z = 0 instead — and we want to know exactly what that costs."
10. `nb.note` — **N12 [C]** "Set-up (the book says 'y = 0' once for the still surface; it means z = 0): a liquid of
    constant density and uniform depth H, small slope a/λ ≪ 1 and small amplitude a/H ≪ 1, no surface tension (C07 adds
    it), the air ignored, irrotational motion caused only by the waves; x along the propagation, z up."
11. `nb.figure` — **Fig. 7.2 remake (N168)** with `tank(ax, H=1.0, L=4.0)` + `wave_surface`: the still level z = 0
    (dashed), the surface η = a cos kx (blue), the bottom z = −H (hatched), a, H and η(x, t) labelled, g pointing down.
    *see:* "the geometry of every problem in §7.2–§7.6." *read:* "z = 0 is where the surface would be at rest; the
    conditions of this block are first written on the blue curve, then moved to the dashed line." *change:* "…H → ∞:
    the bottom disappears and (7.12) becomes 'φ → 0 as z → −∞' (deep water, D08)."
12. `nb.md` — **The idea**: "Any smooth quantity evaluated on the surface can be Taylor-expanded about the flat level:
    ```
    condition lives on z = η(x, t)            (unknown, moving)
    F(η) = F(0) + η F′(0) + …                 (Taylor about z = 0)
    keep F(0):  error ≈ η F′ ≈ ka × (kept term)  ⇒  negligible when ka ≪ 1
    ```
    **Linear = drop products of small quantities, then apply the conditions on z = 0.**"
13. `nb.primer("Taylor transfer of a boundary condition", "A condition on a moving boundary z = η can be written on the
    fixed level z = 0 by expanding each quantity F about z = 0: F(η) ≈ F(0) + ηF′(0). If F varies on the scale 1/k and
    η ~ a, the correction is about ka times F — so for gentle waves we keep only F(0), and both of the book's
    approximations (7.17) and (7.18) are errors of the same order ka.", code="import numpy as np   # numbers\nk, eta
    = 1.0, 0.05   # decay rate [1/m] and a small surface height [m]\nF = lambda z: np.exp(k*z)   # a quantity that varies
    on the scale 1/k\nprint(F(eta), F(0) + eta*k*F(0), F(0))   # 1.05127, 1.05, 1.0: dropping ηF′ costs 5 % = ka")`
    (**P166**).
14. `nb.derivation("D02", …)` — Part F D02 (8 steps), ref "7.16".
15. `nb.note` — **N13 [B]** "The upward unit normal of the surface f = z − η(x, t) = 0 is
    $\mathbf n=\nabla f/\lvert\nabla f\rvert=(-(\partial\eta/\partial x)\mathbf e_x+\mathbf e_z)/\sqrt{(\partial\eta/
    \partial x)^2+1}$ (7.14) (D02 step 2). ⚠️ In Ch. 4 the level-set function was called η; here η is the height and the
    function is z − η." equation (7.14) + `nb.code`: `print(ch07.surface_normal(0.1))` → (−0.0995, 0.9950). **N14 [B]**
    "$\mathbf U_s=(\partial\eta/\partial t)\mathbf e_z$ (7.15): only the normal part of the surface's velocity matters, so
    we may follow the surface point straight above x." **N15 [B]** "The exact kinematic condition
    $\big(\frac{\partial\phi}{\partial z}\big)_{z=\eta}=\frac{\partial\eta}{\partial t}+\frac{\partial\eta}{\partial x}
    \big(\frac{\partial\phi}{\partial x}\big)_{z=\eta}$ (7.16) is D02's result — it is $D(z-\eta)/Dt=0$."
16. `nb.derivation("D03", …)` — Part F D03 (8 steps), ref "7.18".
17. `nb.note` — **N16 [B]** "Dropping the slope × velocity product gives $(\partial\phi/\partial z)_{z=\eta}\cong\partial
    \eta/\partial t$ (7.17) (D03 step 3): the dropped term is about ka·aω, the kept ones aω."
18. `nb.derivation("D04", …)` — Part F D04 (8 steps), ref "7.21".
19. `nb.note` — **N17 [B]** "The linearised dynamic condition $(\partial\phi/\partial t)_{z=\eta}\cong(\partial\phi/
    \partial t)_{z=0}\cong-g\eta$ (7.21) is D04's result." equation (7.21).
20. `nb.md` — **The linear problem in one box** (the whole of C03 starts from it):
    $$\nabla^2\phi=0\ \text{(7.11)},\qquad \frac{\partial\phi}{\partial z}=0\ \text{on } z=-H\ \text{(7.12)},\qquad
    \Big(\frac{\partial\phi}{\partial z}\Big)_{z=0}=\frac{\partial\eta}{\partial t}\ \text{(7.18)},\qquad
    \Big(\frac{\partial\phi}{\partial t}\Big)_{z=0}=-g\eta\ \text{(7.21)}.$$
    "Four linear equations on fixed boundaries: sums of solutions are solutions (N02)."
21. `nb.worked_example("how much do we throw away?", "A swell with a = 0.5 m, λ = 50 m. 1. $k=2\\pi/50=0.126$ rad/m. 2.
    $ka=0.063$. 3. The kept terms of the kinematic condition are about $a\\omega$; the dropped slope × velocity term is
    about $ka\\cdot a\\omega$ — 6 % of the kept one. 4. A steep wave with ka = 0.3 (close to breaking, C12) loses 30 %:
    linear theory is then only a first guess.")`
22. `nb.code` — `for a in (0.05, 0.3): r = ch07.free_surface_residuals(0.3, 0.0, a, 1.0, H=2.0); print(a, {name:
    f"{v:.2e}" for name, v in r.items()})`. *expect:* for both a: `laplace`, `bottom`, `kinematic_linear`,
    `dynamic_linear` ≲ 1e-8 (central-difference round-off: the linear solution solves the linear problem exactly);
    `kinematic_exact`, `kinematic_17` and `dynamic_exact` of order 1e-3 at a = 0.05 and ≈ 36× larger at a = 0.3 (six
    times the amplitude, squared). *explain:* 1. the residual of each condition of the box, evaluated for the linear wave
    of C03 (a, k = 1 rad/m, H = 2 m); 2. the linear conditions are met to round-off, the exact ones (on z = η, with the
    dropped products) are not; 3. their misfit grows like a².
23. `nb.check_agree` — **from scratch:** the exact kinematic residual (7.16) at one point by hand: `h = 1e-6`; `F = lambda
    x, z, t: ch07.wave_fields(x, z, t, 0.05, 1.0, H=2.0)`; `eta = F(x, 0, t)["eta"]`; `eta_t = (F(x, 0, t+h)["eta"] -
    F(x, 0, t-h)["eta"])/(2*h)`; `eta_x = (F(x+h, 0, t)["eta"] - F(x-h, 0, t)["eta"])/(2*h)`; `res = F(x, eta, t)["w"] -
    eta_t - eta_x*F(x, eta, t)["u"]`; `assert np.allclose(res, ch07.free_surface_residuals(x, t, 0.05, 1.0, H=2.0)[
    "kinematic_exact"], rtol=1e-5)`. Markdown: "Our three-line evaluation of (7.16) on the true surface matches the
    library's residual."
24. `nb.figure` — (6.5 × 3.4 in) log–log: `s = ch07.free_surface_residual_scan(np.logspace(-3, -0.5, 12), kH=1.0)`;
    `kinematic_abs` and `dynamic_abs` (black, solid and dashed) with slope-2 guide, `kinematic_rel` (grey) with slope-1
    guide; vertical markers "swell ka = 0.063" and "steep ka = 0.3"; printed slopes from `observed_order` (expect 2.00 and
    1.00). Title "What linearisation drops grows like (ka)²". *see:* "two straight lines of slope 2 (absolute misfit) and
    one of slope 1 (misfit per unit of kept term)." *read:* "at the swell marker the relative error is ~6 %, at the
    steep-wave marker ~30 %: the approximation is as good as ka is small." *change:* "…the waves were in shallow water
    (kH = 0.2): the same slopes, but a/H ≪ 1 becomes the stricter condition (the bottom term of D03 step 7)."
25. `nb.md` — "> ⚠️ **Common confusion:** 'linear' does not mean 'the flow is slow'. It means every product of two wave
    quantities (slope × velocity, velocity²) is dropped **and** the conditions are moved to z = 0; both errors are of the
    same order ka, so doing one without the other is inconsistent."
26. `nb.md` — **What would change if…** "…the surface had surface tension? The dynamic condition gains a curvature term,
    $(\partial\phi/\partial t)_{z=0}=\frac\sigma\rho\frac{\partial^2\eta}{\partial x^2}-g\eta$ (7.55) (C07). …there were a
    second fluid above? Both conditions become two-sided, (7.93)–(7.94) (C13) — and with a shear flow they give the
    Kelvin–Helmholtz instability of Ch. 11."
27. `nb.core("C03", "The dispersion relation $\\omega=\\sqrt{gk\\tanh kH}$ (7.28)", question="Given a wavelength and a
    depth, how fast must the surface oscillate — and why is there only one answer?")`
28. `nb.md` — **The problem in plain words:** "In a wave tank, push the paddle slowly and long waves come out; push it
    fast and short ones do. The water decides the wavelength that goes with each frequency. We want that rule — the
    dispersion relation — because it is the chapter's hub: speeds, orbits, energy flux, seiches, refraction and group
    velocity are all read off it." + recap sentence: "Ch. 4 already *used* $\omega^2=gk\tanh kH$ as a test field
    (`ch04.linear_wave_surface`); here we derive it for the first time."
29. `nb.md` — **The idea** (ASCII):
    ```
    try   φ = f(z) · sin(kx − ωt)        (the surface is a cosine, so φ must be a sine)
    Laplace (7.11)          ⇒  f″ = k² f      ⇒  f = A e^{kz} + B e^{−kz}
    bottom (7.12)           ⇒  B/A fixed      ⇒  f ∝ cosh k(z + H)
    kinematic (7.18)        ⇒  size fixed by a and ω
    dynamic (7.21)          ⇒  only one ω for each k   ←  the dispersion relation
    ```
30. `nb.primer("separation of variables for a PDE", "When a linear PDE and its boundaries are straight in x and z, try a
    product: one function of z times a wave in x. Substituting turns ∂²/∂x² into −k², so the PDE becomes an ordinary
    differential equation for the z-part (Ch. 1 P44 solved those with an e^{λz} trial). The x-part is chosen by the
    boundary conditions — here the surface cosine.", code="import sympy as sp   # symbolic algebra\nx, z, k = sp.symbols('x
    z k', positive=True)\nf = sp.Function('f')   # the unknown depth profile\nphi = f(z)*sp.sin(k*x)   # product trial\n
    print(sp.simplify((sp.diff(phi, x, 2) + sp.diff(phi, z, 2))/sp.sin(k*x)))   # -k**2*f(z) + f''(z): an ODE in z")`
    (**P167**).
31. `nb.primer("hyperbolic functions cosh, sinh, tanh", "cosh x = (eˣ + e⁻ˣ)/2, sinh x = (eˣ − e⁻ˣ)/2, tanh = sinh/cosh,
    coth = 1/tanh, sech = 1/cosh. Rules used in this chapter: cosh² − sinh² = 1, cosh² + sinh² = cosh 2x, 2 sinh x cosh x
    = sinh 2x, (cosh)′ = sinh, (sinh)′ = cosh, (tanh)′ = sech². Small x: cosh ≈ 1, sinh ≈ tanh ≈ x. Large x: cosh ≈ sinh ≈
    eˣ/2, tanh → 1 (tanh 2 = 0.964). Here x = kH: tanh kH switches between shallow and deep water.", code="import numpy as
    np   # numbers\nx = np.array([0.1, 1.0, 2.0, 5.0])\nprint(np.cosh(x), np.sinh(x), np.tanh(x))   # tanh: 0.0997 0.762
    0.964 1.000\nprint(np.cosh(x)**2 - np.sinh(x)**2)   # [1. 1. 1. 1.] — the identity")` (**P168**).
32. `nb.note` — **N18 [C]** "The linear problem also needs an initial surface shape; the book picks η(x, 0) = a cos kx
    (the general case is N02's Fourier sum, Exercise 7.3)."
33. `nb.note` — **N19 [B]** "The separable trial $\phi(x,z,t)=f(z)\sin(kx-\omega(k)t)$ (7.22) is D05's step 1." equation
    (7.22), ref "7.22".
34. `nb.derivation("D05", …)` — Part F D05 (11 steps), ref "7.26".
35. `nb.note` — **N20 [B], N21 [B], N22 [B], N23 [B]** "D05 passes through $\frac{d^2f}{dz^2}-k^2f=0$, $f=Ae^{kz}+Be^{-kz}$
    and $\phi=(Ae^{kz}+Be^{-kz})\sin(kx-\omega t)$ (7.23) (steps 2–4), the bottom condition $k(Ae^{-kH}-Be^{+kH})\sin(kx-
    \omega t)=0$, i.e. $B=Ae^{-2kH}$ (7.24) (step 5), the kinematic condition $k(A-B)=\omega a$ (7.25) (step 7), and ends at
    $\phi=\frac{a\omega}{k}\frac{\cosh k(z+H)}{\sinh kH}\sin(kx-\omega t)$ (7.26)." equations (7.23)–(7.26).
36. `nb.note` — **N24 [B]** "One derivative each gives the velocities $u=a\omega\frac{\cosh k(z+H)}{\sinh kH}\cos(kx-
    \omega t)$, $w=a\omega\frac{\sinh k(z+H)}{\sinh kH}\sin(kx-\omega t)$ (7.27) (∂/∂x of the sine gives k cos; d/dz of
    cosh k(z + H) gives k sinh k(z + H)). They came from the kinematic conditions alone: the pressure follows afterwards
    from Bernoulli." + `nb.code`: `x, z, t = 0.7, -0.4, 1.3` · `f7 = ch07.wave_fields(x, z, t, 0.05, 1.0, H=2.0)` ·
    `f4 = ch04.linear_wave_surface(x, z, t, a=0.05, k=1.0, H=2.0)` · `assert np.allclose([f7["u"], f7["w"]], [f4["u"],
    f4["w"]], rtol=1e-12)` (key names as in the ch04 function). "Ch. 4's test field is exactly this derived solution."
37. `nb.derivation("D06", …)` — Part F D06 (7 steps), ref "7.28".
38. `nb.note` — **N25 [B]** "The line that produced (7.28): $\big(\frac{\partial\phi}{\partial t}\big)_{z=0}=-\frac{a
    \omega^2}{k}\frac{\cosh kH}{\sinh kH}\cos(kx-\omega t)\cong-g\eta=-ag\cos(kx-\omega t)$ (D06 steps 2–4)." equation.
39. `nb.worked_example("H = 1 m, k = 1 rad/m", "1. $kH=1$, $\\tanh 1=0.762$. 2. $\\omega^2=gk\\tanh kH=9.81\\times1
    \\times0.762=7.47$ s⁻². 3. $\\omega=2.73$ rad/s. 4. $T=2\\pi/\\omega=2.30$ s. 5. $\\lambda=2\\pi/k=6.28$ m and
    $c=\\omega/k=2.73$ m/s. Check with deep water: $\\sqrt{gk}=3.13$ rad/s would be 15 % faster — the 1 m bottom slows this
    6 m wave.")`
40. `nb.code` — `print(ch07.omega_gravity(1.0, 1.0))` · `print(ch07.period_from_wavelength(50.0, 10.0),
    ch07.period_from_wavelength(156.0))` · `k = ch07.wavenumber_from_omega(2*np.pi/10, 30.0); print(k, 2*np.pi/k)` ·
    `s = ch07.surface_wave_sympy(); print(s["dispersion"], s["laplace_residual"], s["bottom_residual"],
    s["kinematic_residual"], s["regroup_identity"])`. *expect:* 2.733 rad/s; 6.138 s (H = 10 m, λ = 50 m: kH = 1.257)
    and 9.996 s (deep, λ = 156 m); k = 0.045764 rad/m, λ = 137.3 m; `Eq(omega**2, g*k*tanh(H*k))` and four zeros.
    *explain:* 1. the forward relation ω(k, H); 2. its T–λ form $T=\sqrt{\frac{2\pi\lambda}{g}\coth\frac{2\pi H}{\lambda}}$
    (7.28); 3. the inverse by `brentq` between the deep and shallow roots (P108 reminder: a bracket where the function
    changes sign); 4. sympy re-does D05–D06 and confirms every condition. Gloss in the explain: "cosh and sinh overflow
    near kH ≈ 710; the library writes the ratio as $e^{kz}(1+e^{-2k(z+H)})/(1-e^{-2kH})$ and uses `np.tanh`, which
    saturates safely at 1."
41. `nb.check_agree` — **from scratch (Newton's method, T = 10 s, H = 30 m):** `w0 = 2*np.pi/10; k = w0**2/G` (deep
    guess) · `for i in range(5): f = w0**2 - G*k*np.tanh(k*30); fp = -(G*np.tanh(k*30) + G*k*30/np.cosh(k*30)**2); k -=
    f/fp; print(i, k)` · `assert np.allclose(k, ch07.wavenumber_from_omega(w0, 30.0), rtol=1e-10)`. *expect:* 0.045749,
    0.0457642, 0.045764159 (three iterations to 10 digits). Gloss: "Newton: k ← k − f/f′ (Ch. 6 P152 did it for complex
    zeros)."
42. `nb.code` — **V5 cross-check** (published explicit approximations, Fenton 2006): `for H in (2.0, 30.0, 300.0): w =
    2*np.pi/10; exact = ch07.wavenumber_from_omega(w, H)*H; print(H, exact, ch07.fenton_mckee_kh(w, H), ch07.guo_kh(w,
    H))`. *expect:* relative differences ≤ 1.5 % (Fenton–McKee) and ≤ 0.7 % (Guo), both exact in the two limits.
43. `nb.figure` — (6.5 × 3.8 in) ω(k) for H = 1, 10, 100 m (three blues), k ∈ [10⁻³, 1] rad/m on log–log axes, the deep
    asymptote √(gk) (muted dashed) and the shallow lines k√(gH) (muted dotted, one per H), a dot at kH = 1 on each curve.
    Title "Each depth follows k√(gH) for long waves, then joins √(gk)". *see:* "three curves that start on straight lines
    of slope 1 and bend onto the common slope-½ line." *read:* "left of kH ≈ 0.3 the wave feels the bottom (ω ∝ k: all
    long waves move together); right of kH ≈ 3 it does not (ω ∝ √k)." *change:* "…the depth were halved: each curve's
    bend moves to twice the k; the deep line does not move at all."
44. `nb.md` — **What would change if…** "…the surface had tension σ? Every g in D06 becomes g + σk²/ρ, (7.56) (C07).
    …there were a lighter fluid above instead of air? g is multiplied by (ρ₂ − ρ₁)/(ρ₂ + ρ₁), (7.95) (C13)."
45. `nb.recap("R07", "Perturbation pressure", "Split the pressure into its still-water part and the wave's part:
    $p'\\equiv p+\\rho gz$ (7.30) (in gauge pressure the still water has p = −ρgz). Ch. 4 made the same split, p = p_s + p′,
    for the Boussinesq equations.", where="Ch. 4 §4.9")`
46. `nb.core("C04", "Phase speed $c=\\sqrt{\\frac gk\\tanh kH}$ (7.29): dispersion, deep water and shallow water",
    question="When does the depth matter — and why do long waves all travel together?")`
47. `nb.md` — **The problem in plain words:** "A tsunami crosses the Pacific in less than a day; the swell from a storm in
    the Southern Ocean needs about a week to reach California. Both are gravity waves on the same water. The difference is
    the wavelength compared with the depth — we want the rule, its two limits and the numbers."
48. `nb.md` — **The idea** (table):
    | regime | condition (book) | tanh kH ≈ | phase speed | pressure below | dispersive? |
    |---|---|---|---|---|---|
    | deep | kH > 2 (H > 0.32λ) | 1 | $\sqrt{g/k}$ (7.45) | decays like $e^{kz}$ (7.48) | yes: c ∝ √λ |
    | intermediate | 0.44 < kH < 2 | tanh kH | (7.29) | cosh ratio (7.31) | yes |
    | shallow | H < 0.07λ (kH < 0.44) | kH | $\sqrt{gH}$ (7.49) | hydrostatic ρgη (7.52) | no |
    "**Depth is invisible to a wave shorter than about twice the depth; for longer waves it sets the speed.**"
49. `nb.note` — **N37 [B]** "The two limits of the hyperbolic functions (primer P168): tanh x → 1 for large x (0.964 at
    x = 2 — the deep-water threshold), tanh x ≈ sinh x ≈ x and cosh x ≈ 1 for small x." + `nb.figure` **Fig. 7.7 remake
    (N173)** cosh, sinh, tanh on 0 ≤ x ≤ 2.3 with the dashed line y = 1 and a marker at x = 2 (tanh 2 = 0.964). *see:*
    "cosh starts at 1, sinh and tanh at 0 with slope 1; tanh flattens under 1." *read:* "left: sinh ≈ tanh ≈ x (shallow
    water); right: tanh ≈ 1 (deep water)." *change:* "…plot to x = 5: cosh and sinh become indistinguishable (both ≈ eˣ/2)."
50. `nb.derivation("D07", …)` — Part F D07 (5 steps), ref "7.29".
51. `nb.derivation("D08", …)` — Part F D08 (7 steps), ref "7.45".
52. `nb.note` — **N38 [B]** "Deep water: $c=\sqrt{g/k}=\sqrt{g\lambda/2\pi}$ (7.45), within 2 % once kH > 2 (H > 0.32λ,
    'deeper than a third of a wavelength'); still dispersive — c grows like √λ." **N42 [B]** "Deep-water pressure
    $p'=\rho gae^{kz}\cos(kx-\omega t)$ (7.48): at depth λ/2 only $e^{-\pi}=4.3\,\%$ is left, so a bottom pressure
    gauge in deep water is a low-pass filter — it sees swell only when the water is shallow compared with λ."
53. `nb.derivation("D09", …)` — Part F D09 (7 steps), ref "7.49".
54. `nb.note` — **N43 [B]** "Shallow water: $c=\sqrt{gH}$ (7.49) — the same speed for every long wave (nondispersive),
    within 3 % when H < 0.07λ; it matches the control-volume bore speed of Ch. 4 Example 4.3 in the small-step limit." +
    `nb.code`: `print(ch07.phase_speed(2*np.pi/2e5, 4000.0), np.sqrt(G*4000), ch04.bore_speed(4000.0, 4000.001))`.
    *expect:* 197.57, 198.09 and 198.09 m/s (the Ch. 4 bore speed for a vanishing step, same g). **N46 [B]**
    "$p'=\rho ga\cos(kx-\omega t)=\rho g\eta$ (7.52): in shallow water the pressure at every depth is just the weight of
    the extra water above — **shallow waves are hydrostatic**. This is the hydrostatic shallow-water model of Ch. 13
    (tides, tsunamis, Kelvin waves)."
55. `nb.note` — **N26 [B]** "The full pressure under a wave, from (7.20) with (7.30) and (7.26), then (7.28) to remove ω²:
    $p'=-\rho\frac{\partial\phi}{\partial t}=\rho\frac{a\omega^2}{k}\frac{\cosh k(z+H)}{\sinh kH}\cos(kx-\omega t)=\rho ga
    \frac{\cosh k(z+H)}{\cosh kH}\cos(kx-\omega t)$ (7.31) (two substitutions; the book writes them in one line). It
    decays with depth at a rate set by k: short waves are felt only near the surface." equation (7.31).
56. `nb.worked_example("swell, a tsunami and a bottom gauge", "1. T = 10 s in the open ocean (deep): $\\lambda=gT^2/2\\pi=
    9.81\\times100/6.283=156$ m (the book rounds to 150 m) and c = λ/T = 15.6 m/s. 2. The same swell on a 30 m shelf:
    solve (7.28) → λ = 137 m, c = 13.7 m/s (kH = 1.37: it has started to feel the bottom). 3. A tsunami, H = 4 km:
    $c=\\sqrt{gH}=\\sqrt{39\\,240}=198$ m/s ≈ 713 km/h — the Pacific (≈ 15 000 km) in ≈ 21 h. 4. λ = 50 m, a = 1 m on
    H = 10 m: ρga = 9.81 kPa at z = 0; at the bottom divide by cosh kH = cosh 1.257 = 1.90 → 5.17 kPa.")`
57. `nb.note` — **N39 [B]** "Ocean scales: wind waves of T ≈ 10 s are ≈ 156 m long, so over a 100 m shelf and a 4 km
    ocean they are deep-water waves until close to shore; tides and tsunamis (hundreds of km) are always shallow-water
    waves." + `nb.code`: `for H in (np.inf, 4000.0, 100.0, 30.0, 10.0): lam = ch07.wavelength_from_period(10.0, H);
    print(H, round(lam, 1), ch07.depth_regime(2*np.pi/lam, H))`. *expect:* λ = 156.1, 156.1, 156.0, 137.3, 92.4 m; regimes
    deep, deep, deep (kH = 4.03), intermediate (1.37), intermediate (0.68), with `deep_error` and `shallow_error` printed.
58. `nb.code` — `for lam in (5.0, 50.0, 500.0): k = 2*np.pi/lam; print(lam, ch07.pressure_response(k, np.array([0.0,
    -5.0, -10.0]), 10.0))`. *expect:* λ = 5 m: [1, 0.0019, 0.0000] (a bottom gauge sees nothing), λ = 50 m: [1, 0.63,
    0.53], λ = 500 m: [1, 0.99, 0.99] (hydrostatic). *explain:* the cosh ratio of (7.31) at three depths.
59. `nb.figure` — two panels (9 × 3.6 in). Left: $c/\sqrt{gH}$ vs kH on log x (0.01…10) with the deep limit
    $\sqrt{\tanh kH/kH}\to1/\sqrt{kH}$ (muted dashed) and the shallow limit 1 (muted dotted), shaded bands where each limit
    is within 2 % / 3 %, the book's thresholds kH = 0.44 and kH = 2 as vertical lines. Right: $p'/\rho ga$ vs z/H for kH =
    0.3, 1, 3 (amber shades) with the hydrostatic line 1 (dotted) and $e^{kz}$ (dashed). Title "Where each limit holds".
    *see:* "the speed curve leaves the shallow value near kH ≈ 0.4 and joins the deep curve near kH ≈ 2; the pressure
    profiles go from vertical lines to fast decays." *read:* "read the error of a limit as the gap between the solid curve
    and the ghost: 3 % at kH = 0.44, 1.8 % at kH = 2." *change:* "…the depth doubled at fixed λ: every point slides right
    by log 2 — the wave becomes 'deeper'."
60. `nb.plotly` — **IF1:** `slider_figure(lambda H: {"finite depth": (lam, ch07.phase_speed(2*np.pi/lam, H)), "deep
    limit (7.45)": (lam, np.sqrt(G*lam/(2*np.pi))), "shallow limit (7.49)": (lam, np.sqrt(G*H) + 0*lam)}, "H",
    np.logspace(np.log10(0.5), np.log10(5000), 30 if not FAST else 12), unit="m", xlabel="wavelength λ [m]", ylabel="c
    [m/s]", title="Depth only matters for waves longer than about 2H")` with `lam = np.logspace(-1, 6, 300)` (log axes set
    on the returned figure). *explain:* every slider position is computed up front (P17), so it works on the page.
61. `nb.plotly` — **IF3:** slider over λ (1 … 1000 m, 25 steps) at H = 20 m: $p'/\rho ga$ vs z from `pressure_response`
    with the hydrostatic and $e^{kz}$ ghosts. Title "When a bottom gauge sees the wave".
62. `nb.live` — `live(lambda lam=100.0, H=20.0: print(ch07.dispersion_state(lam, H)), lam=(1.0, 2000.0, 1.0), H=(0.5,
    4000.0, 0.5))` — prints k, kH, ω, T, c, c_g, the regime and the bottom-pressure fraction (paired with IF1, P47).
63. `nb.explainer("dispersion_relation", heading="Why long waves outrun short ones", why="Three linked pictures of one
    number: the tank with a crest marker moving at c and the pressure fading with depth, the c(λ) curve with its two
    limits, and the pressure profile. Dragging λ across kH ≈ 0.3…3 moves all three at once — a static figure can show
    only one.", tries=["Start at the *wind swell* preset, then drag H down to 10 m: watch λ and c shrink while T stays
    10 s.", "Choose *tsunami*: the status says shallow — c = √(gH) whatever λ you pick.", "Click the tank at the bottom
    for the pressure arithmetic; compare λ = 50 m and λ = 500 m.", "Open the Derivation tab (D06) and step to 'cancel the
    cosine'."])`
64. `nb.md` — "> ⚠️ **Common confusion:** 'deep water' is not a depth in metres. A 10 m pond is deep for 1 m ripples
    and shallow for a 1 km seiche; the ocean is deep for swell and shallow for tides. Always compare H with λ."
65. `nb.md` — **What would change if…** "…the Earth rotated fast enough to matter (periods of hours)? Long waves in
    Ch. 13 add the Coriolis frequency f: $\omega^2=f^2+gHk^2$ (Poincaré waves) — this section's √(gH) is the limit
    f → 0."
66. `nb.recap("R08", "Path lines", "A fluid particle's path x_p(t) obeys $\\frac{dx_p}{dt}=u(x_p,z_p,t)$,
    $\\frac{dz_p}{dt}=w(x_p,z_p,t)$ (7.32) — the path-line equations (3.8); Ch. 3 Example 3.1 linearised an oscillating
    flow exactly as we are about to.", where="Ch. 3 §3.2")` + `nb.code`: `r = ch03.example_3_1(0.5)` (one line; the
    oscillating path line of Ch. 3 is a closed loop too).
67. `nb.core("C05", "Particle orbits: closed ellipses $\\xi^2/A^2+\\zeta^2/B^2=1$ (7.36)", question="If the wave moves but
    the water does not travel, what path does a water parcel follow?")`
68. `nb.md` — **The problem in plain words:** "Throw a cork on swell and release a neutrally buoyant float a few metres
    down: both go round and come back. Divers feel the surge forward under a crest and backward under a trough. We want
    the shape of those loops at every depth — the picture behind Stokes drift (C12), ocean floats (Ch. 13) and sediment
    moved to and fro on the sea bed."
69. `nb.md` — **The idea** (ASCII):
    ```
    deep water:     circles of radius a e^{kz₀}, shrinking fast with depth
    intermediate:   ellipses, flatter as you go down
    shallow water:  flat ellipses of the same width a/kH at every depth, height → 0 at the bottom
    all parcels of one vertical column are at the same point of their loops (in phase); clockwise for a wave to the right
    ```
70. `nb.note` — **N27 [B]** "With (7.27) inserted, the path lines read $\frac{dx_p}{dt}=a\omega\frac{\cosh k(z_p+H)}
    {\sinh kH}\cos(kx_p-\omega t)$, $\frac{dz_p}{dt}=a\omega\frac{\sinh k(z_p+H)}{\sinh kH}\sin(kx_p-\omega t)$ (7.33) —
    nonlinear in the unknown position." equation (7.33).
71. `nb.derivation("D10", …)` — Part F D10 (9 steps), ref "7.35a".
72. `nb.note` — **N28 [B], N29 [B]** "Freezing the right sides at the mean position gives $\frac{d\xi}{dt}\cong a\omega
    \frac{\cosh k(z_0+H)}{\sinh kH}\cos(kx_0-\omega t)$, $\frac{d\zeta}{dt}\cong a\omega\frac{\sinh k(z_0+H)}{\sinh kH}
    \sin(kx_0-\omega t)$ (7.34a, 7.34b) (D10 step 4), and one integration gives the excursions $\xi\cong-a\frac{\cosh
    k(z_0+H)}{\sinh kH}\sin(kx_0-\omega t)$, $\zeta\cong a\frac{\sinh k(z_0+H)}{\sinh kH}\cos(kx_0-\omega t)$ (7.35a,
    7.35b) — purely oscillatory, so the assumption 'the parcel stays near (x₀, z₀)' is self-consistent."
73. `nb.derivation("D11", …)` — Part F D11 (6 steps), ref "7.36".
74. `nb.note` — **N40 [B], N41 [B]** "Deep water (kH > 2): $\xi\cong-ae^{kz_0}\sin(kx_0-\omega t)$, $\zeta\cong
    ae^{kz_0}\cos(kx_0-\omega t)$ (7.46) — circles whose radius halves every $\ln2/k=0.11\lambda$ of depth; the velocity
    $u=a\omega e^{kz}\cos(kx-\omega t)$, $w=a\omega e^{kz}\sin(kx-\omega t)$ (7.47) turns clockwise at ω with constant size
    aωe^{kz}." **N44 [B], N45 [B]** "Shallow water (kH ≪ 1): $\xi\cong-\frac{a}{kH}\sin(kx_0-\omega t)$, $\zeta\cong
    a(1+\frac zH)\cos(kx_0-\omega t)$ (7.50) — flat ellipses of width a/kH at every depth; $u=\frac{a\omega}{kH}\cos(kx-
    \omega t)$, $w=a\omega(1+\frac zH)\sin(kx-\omega t)$ (7.51), so w/u ≲ kH (0.2 for kH = 0.2)."
75. `nb.worked_example("orbits under a 20 m wave in 50 m of water", "1. $k=2\\pi/20=0.314$ rad/m, $kH=15.7$: deep. 2.
    Radius $ae^{kz_0}$ with a = 1 m: 1.00 m at the surface. 3. At z₀ = −2 m: $e^{-0.628}=0.533$ → 0.53 m. 4. At z₀ = −10 m
    (half a wavelength): $e^{-\\pi}=0.043$ → 4 cm. 5. The bottom (−50 m) does not feel the wave: $e^{-15.7}\\approx10^{-7}$.")`
76. `nb.code` — `for kH in (3.0, 1.0, 0.3): print(kH, ch07.orbit_semi_axes(np.array([0.0, -kH/2]), 1.0, 1.0, kH))`
    (k = 1, H = kH) · `p_lin = ch07.particle_path(0.0, -0.5, np.linspace(0, 20, 401), 0.05, 1.0, 1.0, model="linear")` ·
    `p_ex = ch07.particle_path(0.0, -0.5, np.linspace(0, 20, 401), 0.05, 1.0, 1.0, model="exact")` · `print(p_lin["x"][-1],
    p_ex["x"][-1])`. *expect:* kH = 3: A = 1.005, B = 1.000 at the surface; A = 0.235, B = 0.213 at mid-depth; kH = 1:
    1.313/1.000 and 0.960/0.443; kH = 0.3: 3.433/1.000 and 3.321/0.494; focal half-distance a/sinh kH = 0.0998, 0.851,
    3.284 — the same at every depth. The linear path returns to its start; the exact one ends a few mm ahead (the Stokes
    drift of C12). *explain:* 1. semi-axes of (7.36) at two depths; 2. `particle_path` integrates (7.34) or (7.33) with
    `solve_ivp` (DOP853, tight tolerances — gloss: a high-order method so the tiny drift is not numerical error).
77. `nb.check_agree` — **from scratch (RK4 by hand, P95):** step (7.34a, b) with dt = T/400 over one period from ξ = ζ = 0
    at the mean position and compare with `ch07.orbit_linear` shifted by its starting value: `assert np.allclose(mine,
    lib, atol=1e-6)`. Markdown: "RK4 on the frozen-argument equations draws the same ellipse as the closed form (7.35)."
78. `nb.figure` — **Fig. 7.3 and 7.4 remake (N169, N170)** four panels (10 × 3.4 in): (a) one orbit with the mean position
    (x₀, z₀), the excursion (ξ, ζ) as two arrows and a clockwise arrowhead; (b)–(d) kH = 3, 1, 0.3: surface (blue) and
    orbits (teal) at five depths drawn at their positions, bottom hatched. Title "Circles in deep water, flat ellipses in
    shallow water". *see:* "circles that shrink fast with depth; ellipses; flat ellipses reaching the bottom." *read:* "the
    horizontal axis in (d) is the same at every depth (a/kH) and the vertical axis falls linearly to zero at the bottom."
    *change:* "…a doubled: every orbit doubles, the shapes do not change (linear theory) — until ka ~ 0.3, when the loops
    open (C12)."
79. `nb.note` — **N30 [B]** "Integrating u = ∂ψ/∂z in z (and checking −∂ψ/∂x = w) gives the stream function
    $\psi=\frac{a\omega}{k}\frac{\sinh k(z+H)}{\sinh kH}\cos(kx-\omega t)$ (7.37) (Exercise 7.4): ψ = 0 on the bottom and
    on the vertical lines below the zeros of η. Under a crest the water moves forward at every depth, under a trough
    backward." + `nb.code`: `u, w = streamfunction.velocity_from_streamfunction_2d(lambda x, z: ch07.wave_fields(x, z, 0.0,
    0.05, 1.0, H=2.0)["psi"], 0.4, -0.7)` (y ↦ z) vs `ch07.wave_fields(0.4, -0.7, 0.0, 0.05, 1.0, H=2.0)` → equal to 1e-8.
80. `nb.figure` — **Fig. 7.5 remake (N171)** (7 × 3.2 in): contours of (7.37) at t = 0 over one wavelength (kH = 1), the
    ψ = 0 lines (bold) along the bottom and below the zeros of η, arrows of (u, w), the surface on top. *see:* "loops that
    start and end on the surface, forward under the crest and backward under the trough." *read:* "streamlines are a
    snapshot; the particles themselves go round the orbits of the previous figure." *change:* "…t advanced by T/4: the
    whole pattern slides a quarter wavelength right (it moves with the wave)."
81. `nb.animation` — **A1** (curation §6): three panels kH = 3, 1, 0.3; `animate(update, frames=60 if not FAST else 30,
    fig=fig, interval=60)` — the surface (blue) and 5 × 7 tracer particles on their linear orbits (teal dots) with faint
    orbit ghosts; one column's particles joined by a line to show they move in phase. `show_animation(anim)`. *explain:*
    then *see / read / change* in a markdown cell: "the column stays straight and swings as a whole; in shallow water the
    swing is almost purely horizontal."
82. `nb.md` — **What would change if…** "…the amplitude were not small? The parcel is slightly further forward at the top
    of its loop, where the forward velocity is larger, than it is backward at the bottom: the loop does not close. That
    second-order creep is the Stokes drift of C12."
83. `nb.core("C06", "Wave energy $E=\\tfrac12\\rho ga^2$ (7.42), equipartition and the energy flux (7.44)", question="How
    much energy does a wave hold per square metre of sea — and how fast is it delivered to a beach?")`
84. `nb.md` — **The problem in plain words:** "A storm thousands of kilometres away sends swell that breaks on a beach and
    keeps a surfer moving. The energy crossed the ocean without the water crossing it. We want how much energy a wave
    stores (moving water and lifted water) and the power it carries per metre of crest — the numbers wave-energy
    engineers and swell forecasters use."
85. `nb.md` — **The idea:** "Kinetic energy lives in the orbiting water under the surface; potential energy in the humps
    and hollows of the surface itself. Average both over a wavelength and a period:
    ```
    E_k = ½ρg·mean(η²)   (moving water)        E_p = ½ρg·mean(η²)   (lifted water)
    E   = ρg·mean(η²) = ½ρga²                  F = E × (a speed)  —  which speed?
    ```
    Gloss: two kinds of average: the overbar $\overline{\eta^2}=\frac1\lambda\int_0^\lambda\eta^2dx$ over a wavelength
    and ⟨·⟩ over a period; for a sinusoid both give ⟨cos²⟩ = ½ (Ch. 6 P151)."
86. `nb.note` — **N31 [B]** "The kinetic energy per unit horizontal area, depth-integrated and wavelength-averaged:
    $E_k=\frac{\rho}{2\lambda}\int_0^\lambda\int_{-H}^0(u^2+w^2)\,dz\,dx$, which with (7.27) becomes (7.38) (D12 step 3). The
    z-integral stops at 0, not at η: the slab in between changes E only at third order in a." equation (7.38).
87. `nb.derivation("D12", …)` — Part F D12 (10 steps), ref "7.39".
88. `nb.note` — **N32 [B]** "$E_k=\tfrac12\rho g\overline{\eta^2}$ (7.39): gravity appears through (7.28)."
89. `nb.note` — **N33 [B]** "The potential energy is the work to deform the flat surface:
    $E_p=\frac{\rho g}{\lambda}\int_0^\lambda\int_{-H}^\eta z\,dz\,dx-\frac{\rho g}{\lambda}\int_0^\lambda\int_{-H}^0 z\,dz\,dx
    =\frac{\rho g}{\lambda}\int_0^\lambda\int_0^\eta z\,dz\,dx=\frac{\rho g}{2\lambda}\int_0^\lambda\eta^2dx$ (7.40)."
    equation (7.40) + `nb.figure` **Fig. 7.6 sketch (N172)**: a column A in a trough and B on a crest (hatched), dz and dx
    marked. *see:* "lifting column A (mass ρη dx per unit width) by η onto B builds the surface." *read:* "over half a
    wavelength the swaps build a whole wavelength of surface: work ρgη² dx." *change:* "…the fluid above were not air but
    a lighter liquid ρ₁: each swap also moves ρ₁ down — the work uses ρ₂ − ρ₁ (N105)."
90. `nb.derivation("D13", …)` — Part F D13 (7 steps), ref "7.42".
91. `nb.note` — **N34 [B]** "$E_p=\tfrac12\rho g\overline{\eta^2}$ (7.41) $=E_k$: **equipartition** — small oscillations of
    a conservative system share their energy equally between motion and height. It fails once Coriolis forces act (Ch. 13:
    geostrophic adjustment leaves potential energy locked in)."
92. `nb.note` — **N35 [B]** "The energy flux across a vertical line is the pressure work: $F=\big\langle\int_{-H}^0
    pu\,dz\big\rangle=\big\langle\int_{-H}^0p'u\,dz\big\rangle-\rho g\langle u\rangle\int_{-H}^0z\,dz=\big\langle\int_{-H}^0
    p'u\,dz\big\rangle$ (7.43) (the background part drops because ⟨u⟩ = 0 at every depth)."
93. `nb.derivation("D14", …)` — Part F D14 (10 steps), ref "7.44".
94. `nb.note` — **N36 [B]** "$F=\big[\tfrac12\rho ga^2\big]\big[\frac c2\big(1+\frac{2kH}{\sinh2kH}\big)\big]$ (7.44):
    energy × a speed that is **not** c. The book calls it the group speed and returns to it in §7.5 (C09: $c_g=\frac c2
    [1+\frac{2kH}{\sinh 2kH}]$ (7.69)). (Where the book says 'u from (7.28)' before (7.44) it means (7.27): (7.28) is the
    dispersion relation.)"
95. `nb.worked_example("energy and power of a metre of swell", "1. a = 1 m: $E=\\tfrac12\\rho ga^2=\\tfrac12\\times1000
    \\times9.81\\times1=4905$ J/m², half of it kinetic. 2. T = 10 s deep water: c = 15.6 m/s, and the second factor of
    (7.44) is c/2 = 7.81 m/s. 3. $F=4905\\times7.81=38\\,300$ W per metre of crest ≈ 38 kW/m. 4. If the energy moved at c
    we would have claimed 76.6 kW/m — twice too much.")`
96. `nb.code` — `for kH in (0.5, 1.0, 3.0): k = 1.0; H = kH; c = ch07.wave_energy(1.0, k, H, method="closed"); q =
    ch07.wave_energy(1.0, k, H, method="quad"); print(kH, c, q)` · `k = (2*np.pi/10)**2/G; print(ch07.energy_flux(1.0, k),
    ch07.energy_flux(1.0, k, method="quad"))`. *expect:* for each kH: Ek = Ep = 2452.5 J/m², E = 4905 J/m² in both methods
    (quad to ~1e-8 relative); F = 38 291 W/m both ways. *explain:* closed forms (7.39), (7.41), (7.42) vs `dblquad` of the
    integrals (7.38), (7.40); flux (7.44) vs the depth–period average of p′u.
97. `nb.check_agree` — **from scratch:** midpoint double sum of ½ρ(u² + w²) from `ch07.wave_fields` on a 400 × 200 grid
    over one wavelength and the depth (kH = 1), divided by λ: `assert np.allclose(Ek_mine, ch07.wave_energy(1.0, 1.0, 1.0)[
    "Ek"], rtol=1e-4)`.
98. `nb.figure` — two panels (9 × 3.4 in). Left: the depth-integrated kinetic energy density and ½ρgη² along one
    wavelength at t = 0 (teal and blue) with their equal means (dashed). Right: F/(Ec) = c_g/c vs kH (purple, log x) from 1
    (shallow) to ½ (deep), dots at the swell (kH = 4.0) and at the 30 m shelf (kH = 1.37). Title "Energy is shared equally
    and travels slower than the crests". *see:* "two curves with the same average; a ratio falling from 1 to ½." *read:*
    "in deep water energy arrives at half the crest speed (C09 explains why)." *change:* "…the wave were in 2 m of water:
    kH ≈ 0.2, the ratio ≈ 0.99 — energy and crests travel together."
99. `nb.md` — "> ⚠️ **Common confusion:** surface-wave E is energy per square metre of sea surface (the whole water column
    below it), and F is power per metre of crest. For internal waves in C16 the same letters mean per cubic metre and per
    square metre."
100. `nb.md` — **What would change if…** "…the Earth's rotation joined in (long waves, Ch. 13)? Kinetic and potential
     energy stop being equal: geostrophic adjustment keeps some potential energy in a balanced current."
101. `nb.pointer("**Refraction (Figs. 7.8–7.9: waves turning toward a beach and wrapping round an island) belongs to this
     section of the book**; we teach it in C10 (§7.5), where the rays that explain it are derived (N47, N48, N174,
     N175).")`

### A.3 §7.3 Influence of Surface Tension — R09, C07 (+N49–N57, N176 · D15 D16 · E3)
1. `nb.section("7.3", "Influence of Surface Tension", intro="**What is this section about?** A curved surface pulls
   itself flat, like a stretched membrane — a second restoring force besides gravity. It only matters for short waves,
   but there it wins, and together the two forces leave a slowest possible wave on water.")`
2. `nb.recap("R09", "The Laplace pressure jump", "Across a curved interface the pressure jumps by $\\Delta p=\\sigma(1/R_1+
   1/R_2)$ (1.5), higher on the side of the centres of curvature; Ch. 4 re-derived it as a force balance on a small
   patch.", where="Ch. 1 §1.6, Ch. 4 §4.10")` + `nb.code`: `print(interfaces.laplace_jump_from_balance(0.0727, 0.01,
   np.inf))` → 7.27 Pa (a 1 cm radius in one direction only, as under a straight crest).
3. `nb.core("C07", "Capillary–gravity waves: $c=\\sqrt{(\\frac gk+\\frac{\\sigma k}{\\rho})\\tanh kH}$ (7.57) and the
   minimum speed (7.58)", question="Why do the smallest waves on a pond run faster than slightly longer ones — and is
   there a slowest wave?")`
4. `nb.note` — **N49 [C]** "Surface tension acts like a stretched membrane: a second restoring force. Pure capillary waves
   (no gravity) are rare on Earth, so we extend §7.2 rather than replace it."
5. `nb.md` — **The problem in plain words:** "Rain on a lake makes tiny ripples that race ahead of the wind waves; a
   stone makes a ring of ripples with a calm patch in the middle. Short waves feel the surface's tension, long ones
   gravity. We want the speed with both forces and the wavelength where the speed is smallest."
6. `nb.md` — **The idea** (table):
   | restoring force | stronger for | term in c² (deep water) | branch |
   |---|---|---|---|
   | gravity (weight of the hump) | long waves | g/k | $c=\sqrt{g\lambda/2\pi}$ grows with λ |
   | surface tension (curvature) | short waves | σk/ρ | $c=\sqrt{2\pi\sigma/\rho\lambda}$ grows as λ shrinks |
   "**One force favours long waves, the other short ones — so the speed has a minimum in between.**"
7. `nb.primer("curvature of a plane curve", "For a curve z = η(x) the curvature (1/radius of the osculating circle) is
   $1/R=\\eta_{xx}/(1+\\eta_x^2)^{3/2}$; for gentle slopes 1/R ≈ η_xx. Sign: under a crest η_xx < 0 and the centre of
   curvature lies below, inside the water.", code="import numpy as np   # numbers\na, k = 0.001, 2*np.pi/0.01   # a
   1 mm ripple, λ = 1 cm\netaxx_crest = -a*k**2   # η = a cos kx at x = 0: η'' = −ak²\nprint(etaxx_crest,
   1/abs(etaxx_crest))   # −394.8 1/m: radius 2.5 mm at the crest")` (**P169**).
8. `nb.derivation("D15", …)` — Part F D15 (10 steps), ref "7.56".
9. `nb.note` — **N50 [B], N51 [B], N52 [B], N53 [B]** "D15 passes through the curvature condition
   $p_a-(p)_{z=\eta}=\sigma\frac1R=\sigma\frac{\partial^2\eta/\partial x^2}{[1+(\partial\eta/\partial x)^2]^{3/2}}\cong\sigma
   \frac{\partial^2\eta}{\partial x^2}$ (7.53), its gauge form $(p)_{z=\eta}=-\sigma\frac{\partial^2\eta}{\partial x^2}$
   (7.54), the new dynamic condition $\big(\frac{\partial\phi}{\partial t}\big)_{z=0}=\frac\sigma\rho\frac{\partial^2\eta}
   {\partial x^2}-g\eta$ (7.55) and ends at $\omega=\sqrt{k\big(g+\frac{\sigma k^2}{\rho}\big)\tanh kH}$ (7.56): only the
   dispersion relation changes, φ keeps its form (7.26)." equations (7.53)–(7.56).
10. `nb.derivation("D16", …)` — Part F D16 (7 steps), ref "7.58".
11. `nb.note` — **N54 [B]** "$c_{min}=\big[\frac{4g\sigma}{\rho}\big]^{1/4}$ at $\lambda_m=2\pi\sqrt{\frac{\sigma}{\rho g}}$
    (7.58); for λ < λ_m surface tension dominates."
12. `nb.worked_example("clean water at 20 °C", "1. σ = 72.74 mN/m (IAPWS, Ch. 1), ρ = 998.2 kg/m³, g = 9.81 m/s². 2.
    $\\lambda_m=2\\pi\\sqrt{0.07274/(998.2\\times9.81)}=2\\pi\\times2.73\\times10^{-3}=0.0171$ m = 1.71 cm. 3.
    $c_{min}=(4\\times9.81\\times0.07274/998.2)^{1/4}=(2.86\\times10^{-3})^{1/4}=0.231$ m/s = 23.1 cm/s. 4. Under the crest
    of a 1 mm ripple with λ = 1 cm: $k=2\\pi/0.01=628.3$ rad/m, $k^2=3.95\\times10^5$ m⁻², so
    $p=-\\sigma\\eta_{xx}=\\sigma ak^2=0.07274\\times0.001\\times3.95\\times10^5=+28.7$ Pa, pushing the crest down.")`
13. `nb.note` — **N55 [B]** "Air–water at 20 °C (7.59): our numbers above, c_min ≈ 23.1 cm/s at λ_m ≈ 1.71 cm (the book's
    printed values are in the private JSON)." **N56 [C]** "Only ripples shorter than about 7 cm feel σ, and below about
    4 mm gravity hardly matters; soap and oil films lower σ and change the ripples." **N57 [B]** "Dropping g in deep water
    gives pure capillary waves $c=\sqrt{\frac{2\pi\sigma}{\rho\lambda}}$ (7.60): shorter is faster." equation (7.60).
14. `nb.primer("scipy.optimize.minimize_scalar", "Finds the minimum of a function of one variable inside a bracket
    (method='bounded'); we use it to check a formula for a minimum by brute force.", code="from scipy.optimize import
    minimize_scalar   # 1-D minimiser\nr = minimize_scalar(lambda x: (x - 2)**2 + 1, bounds=(0, 5), method='bounded')\n
    print(r.x, r.fun)   # 2.0 1.0")` (**P170**).
15. `nb.code` — `sig, rho = ch01.surface_tension_water(293.15), 998.2` · `print(sig, ch07.capillary_minimum(sig, rho))` ·
    `for lam in (0.004, 0.0171, 0.07, 1.0): k = 2*np.pi/lam; print(lam, ch07.phase_speed(k, np.inf, sigma=sig, rho=rho),
    ch07.phase_speed(k, np.inf))` · `print(ch07.capillary_surface_pressure(-0.001*(2*np.pi/0.01)**2, sig),
    ch07.curvature(0.0, -394.8, linear=True))`. *expect:* σ = 0.07274 N/m; c_min = 0.2312 m/s, λ_m = 0.01712 m; c with σ
    for λ = 4 mm, 1.71 cm, 7 cm, 1 m ≈ 0.35, 0.231, 0.34, 1.249 m/s vs without σ 0.079, 0.163, 0.33, 1.249 m/s; pressure
    +28.7 Pa. *explain:* the branch that wins at each λ.
16. `nb.check_agree` — **from scratch:** scan c(λ) on `np.logspace(-3.5, 0, 4001)`, take the smallest, then refine with
    `minimize_scalar` on c² in k; `assert np.allclose([cmin_mine, lam_mine], [d["c_min"], d["lam_m"]], rtol=1e-6)`.
17. `nb.figure` — **Fig. 7.10 remake (N176)** log–log (6.5 × 4 in): c(λ) for σ = 0.0727 N/m, H = 1 m (black), the capillary
    branch √(2πσ/ρλ) (rose dashed), the deep gravity branch √(gλ/2π) (blue dashed), the shallow plateau √(gH) (muted
    dotted), the minimum marked (dot + labels c_min, λ_m), bands "capillary < 4 mm", "affected by σ < 7 cm" (N56). Title
    "Two restoring forces leave a slowest wave". *see:* "a valley at 1.7 cm between the rising capillary branch and the
    gravity branch, then the √(gH) plateau for very long waves." *read:* "at any λ the speed is the square root of the sum
    of the two branches' squares (times tanh kH)." *change:* "…soap halves σ: the valley moves to 1.2 cm and 19 cm/s
    (both ∝ σ^{1/2} and σ^{1/4})."
18. `nb.plotly` — **IF2:** `slider_figure` over σ ∈ [0, 0.5] N/m (26 steps) of c(λ) with the minimum marker (from
    `capillary_minimum`), log axes. Title "The minimum slides with σ".
19. `nb.explainer("capillary_gravity_waves", heading="Why is there a slowest ripple?", why="Sweep λ over four decades and
    change the liquid: the gravity and surface-tension term bars trade places, the minimum moves, and the status names the
    branch. Only interaction shows the two forces competing.", tries=["Press *λ = λ_m*: the two bars are equal.", "Switch
    the liquid to mercury: the minimum moves to 1.2 cm and 19 cm/s.", "Drag λ to 5 mm and read the crest pressure in
    the Explain tab.", "Step the Derivation (D16) to 'set the slope to zero'."])`
20. `nb.md` — **What would change if…** "…the group speed were asked for instead of the phase speed? It also has a minimum,
    ≈ 17.8 cm/s at λ ≈ 4.4 cm — which is why the ring of ripples from a stone leaves a calm centre (C09, N72)."

### A.4 §7.4 Standing Waves — C08 (+N58–N63, N177, N178 · D17 D18 · A3 · E4)
1. `nb.section("7.4", "Standing Waves", intro="**What is this section about?** Two equal waves travelling in opposite
   directions add into a pattern that does not travel: fixed nodes, sloshing in between. Walls allow only the patterns
   whose nodes fit — so a lake, a harbour or a bathtub rings at its own set of periods (seiches).")`
2. `nb.core("C08", "Standing waves $\\eta=2a\\cos kx\\cos\\omega t$ and seiches (7.64)–(7.65)", question="Which waves can
   live in a closed basin, and with which periods?")`
3. `nb.md` — **The problem in plain words:** "Push the water in a bathtub once and it sloshes back and forth at a fixed
   rhythm. Lake Geneva does the same after a storm, with a period of about an hour; harbours do it and ships break their
   moorings. We want to know which sloshing patterns a basin allows and their periods."
4. `nb.md` — **The idea** (ASCII):
   ```
   right-going  a cos(kx − ωt)  +  left-going  a cos(kx + ωt)  =  2a cos kx · cos ωt
   the shape cos kx never moves; it only breathes in time (cos ωt)
   nodes of η at kx = π/2, 3π/2, …   — walls must sit where u = 0 (antinodes of η)
   ```
5. `nb.primer("sum-to-product identities", "cos A + cos B = 2 cos((A − B)/2) cos((A + B)/2) and cos A − cos B = −2 sin((A +
   B)/2) sin((A − B)/2): a sum of two waves becomes a product of a slow and a fast factor. They turn two travelling waves
   into a standing wave here and into beats in C09.", code="import numpy as np   # numbers\nA, B = 1.1, 0.3   # any two
   angles [rad]\nprint(np.cos(A) + np.cos(B), 2*np.cos((A - B)/2)*np.cos((A + B)/2))   # both 1.4089")` (**P171**).
6. `nb.note` — **N58 [B]** "The left-going wave $\eta(x,t)=a\cos[kx+\omega t]$ (7.61) is (7.2) with ω → −ω: its crests move
   to −x." equation (7.61) + `nb.code`: `print(ch07.sinusoid(0.0, [0, 1], 1.0, 1.0, 1.0, direction=-1))`.
7. `nb.derivation("D17", …)` — Part F D17 (8 steps), ref "7.63".
8. `nb.note` — **N59 [B], N61 [B]** "The standing wave's stream function $\psi=\frac{a\omega}{k}\frac{\sinh k(z+H)}{\sinh kH}
   [\cos(kx-\omega t)-\cos(kx+\omega t)]=\frac{2a\omega}{k}\frac{\sinh k(z+H)}{\sinh kH}\sin kx\sin\omega t$ (7.62) and its
   velocity $u=2a\omega\frac{\cosh k(z+H)}{\sinh kH}\sin kx\sin\omega t$ (7.63): where η has a node, u is largest."
   equations (7.62), (7.63).
9. `nb.figure` — **Fig. 7.11 remake (N177)** contours of (7.62) at ωt = π/2 (largest flow) over one wavelength, kH = 1, the
   surface at ωt = π (drawn with exaggerated amplitude), ψ = 0 lines bold. *see:* "flow from under the crests into the
   troughs' side, vertical lines at crests and troughs." *read:* "the vertical ψ = 0 lines at kx = 0, π are where a wall
   can stand (u = 0 there)." *change:* "…the snapshot were at ωt = 0: η is largest and the flow is zero everywhere — the
   energy is all potential at that instant."
10. `nb.note` — **N60 [B]** "A seiche is a standing oscillation of a closed basin (length L, depth H, vertical walls): only
    patterns with u = 0 at both walls survive — an eigenvalue problem (the idea of Ch. 2 P80: the boundary picks discrete
    solutions)."
11. `nb.derivation("D18", …)` — Part F D18 (5 steps), ref "7.65".
12. `nb.note` — **N62 [B], N63 [B]** "Walls at x = 0 and L allow $kL=(n+1)\pi$, $\lambda=\frac{2L}{n+1}$ (7.64) — the longest
    is 2L, then L, 2L/3 … — and the lake's natural frequencies $\omega=\sqrt{\frac{\pi g(n+1)}{L}\tanh\big[\frac{(n+1)\pi
    H}{L}\big]}$ (7.65). A rectangular basin L × b rings with $k^2=(m\pi/L)^2+(n\pi/b)^2$ (Exercises 7.5–7.6)."
13. `nb.worked_example("a bathtub and a lake", "1. Bathtub L = 1.5 m, H = 0.2 m, n = 0: $k=\\pi/L=2.094$ rad/m, kH = 0.419,
    tanh = 0.396. 2. $\\omega^2=9.81\\times2.094\\times0.396=8.14$ → ω = 2.85 rad/s, T = 2.20 s. 3. Shallow estimate
    $T=2L/\\sqrt{gH}=3/1.40=2.14$ s (3 % short). 4. A lake L = 50 km, H = 100 m: kH = 0.0063 — shallow;
    $T_0=2L/\\sqrt{gH}=100\\,000/31.3=3193$ s ≈ 53 min.")`
14. `nb.code` — `print(ch07.seiche_modes(1.5, 0.2, 0), ch07.seiche_modes(50e3, 100.0, 0))` · `for n in range(4):
    print(n, ch07.seiche_modes(50e3, 100.0, n)["T"]/60)` · `print(ch07.basin_modes(1.5, 0.7, 0.2, 1, 1))` ·
    `x = np.linspace(0, 2, 5); s = ch07.standing_wave_fields(x, -0.1, 0.4, 0.01, np.pi, 0.2)` ·
    `r = ch07.wave_fields(x, -0.1, 0.4, 0.01, np.pi, 0.2); l = ch07.wave_fields(x, -0.1, 0.4, 0.01, np.pi, 0.2,
    direction=-1)` · `assert np.allclose(s["eta"], r["eta"] + l["eta"]) and np.allclose(s["u"], r["u"] + l["u"])`.
    *expect:* T = 2.203 s and 3192.8 s; lake modes 53.2, 26.6, 17.7, 13.3 min (shallow: T ∝ 1/(n + 1)); the basin mode;
    the assert passes. *explain:* the standing wave is exactly the sum of the two travelling waves of C03.
15. `nb.check_agree` — **from scratch:** `Ts = [2*np.pi/np.sqrt(G*(n+1)*np.pi/L*np.tanh((n+1)*np.pi*H/L)) for n in
    range(5)]` for the bathtub · `assert np.allclose(Ts, [ch07.seiche_modes(L, H, n)["T"] for n in range(5)], rtol=1e-12)`.
16. `nb.figure` — **Fig. 7.12 remake (N178)** u(x) at the surface for n = 0 and n = 1 in a basin of length L (walls drawn),
    the η(x) shapes above as thin lines, nodes marked. *see:* "u is zero at both walls; mode 0 has one bulge of u, mode 1
    two opposite ones." *read:* "where u is largest η has its node: water rushes through the node from one side to the
    other." *change:* "…n = 2: three bulges of u and λ = 2L/3; the period shortens to T₀/3 in a shallow lake."
17. `nb.plotly` — **IF7:** `slider_figure` over n = 0…5: η(x) and u(x) of seiche mode n in a unit basin, nodes marked.
18. `nb.animation` — **A3:** standing wave in a basin (kH = 1, n = 1) with streamlines of (7.62) redrawn each frame and the
    two travelling components as faint orange ghosts (`frames=48 if not FAST else 24`). Then *see / read / change*: "the
    ghosts slide through each other; their sum never moves sideways; the streamlines reverse every half period."
19. `nb.explainer("seiche_standing_waves", heading="Which waves can live in a lake?", why="Toggle the two travelling waves
    and watch their sum pin its nodes; then change the basin and the mode and compare the period with real basins in the
    table.", tries=["Switch off the left-going wave: the pattern starts to travel.", "Choose the 50 km lake and step n
    from 0 to 3: the period drops as 1/(n + 1).", "Drag H up in the bathtub until the shallow estimate fails by 10 %."])`
20. `nb.md` — **What would change if…** "…the basin were a bay open to the sea? Its open end becomes a node of η, so the
    longest wave is 4L (a quarter wave) — the harbour-resonance rule; tides in long gulfs resonate this way (Ch. 13)."

### A.5 §7.5 Group Velocity, Energy Flux, and Dispersion — C09 (+N64–N73, N179–N182 · D19 D20 D21 · A2 · E5), C10 (+N47, N48, N74–N80, N174, N175, N183, N184 · D22 D23 D24 · E6)
1. `nb.section("7.5", "Group Velocity, Energy Flux, and Dispersion", intro="**What is this section about?** When the speed
   depends on the wavelength, a group of waves moves at a different speed from its crests. That group speed, dω/dk, is
   the speed of the energy — the answer to C06's cliff-hanger. Then we let the depth change slowly: frequency stays
   constant along rays while wavelength and direction change, which is why waves turn toward beaches.")`
2. `nb.core("C09", "Group velocity $c_g=d\\omega/dk$ (7.67) and the energy flux $F=Ec_g$ (7.71)", question="Why do crests
   appear at the back of a group of swell and vanish at its front — and which speed carries the energy?")`
3. `nb.note` — **N64 [C]** "Speeds that depend on wavelength are common for waves on interfaces between materials
   (Rayleigh and Stoneley waves in solids, liquid–liquid interfaces); we stay with water."
4. `nb.md` — **The problem in plain words:** "Surfers wait for 'sets': swell arrives in groups of a few big waves, then a
   lull. Watch one crest in a group and it runs forward through the group and fades at the front, while new crests grow
   at the back. The group — and the energy — is slower than the crests. We want that speed and the proof that energy
   rides with it."
5. `nb.md` — **The idea** (ASCII):
   ```
   ω(k) curve:   slope of the chord from the origin to (k, ω)   = ω/k   = c    (crests)
                 slope of the tangent at (k, ω)                 = dω/dk = c_g  (groups, energy)
   deep water:   ω = √(gk) bends down  ⇒  tangent flatter than chord  ⇒  c_g = c/2
   ```
6. `nb.derivation("D19", …)` — Part F D19 (7 steps), ref "7.67". (Gloss in D19's *why*: a derivative is the limit of the
   difference quotient Δω/Δk, the chord slope → tangent slope (Ch. 1 P19).)
7. `nb.note` — **N65 [B]** "Two nearby waves make beats: $\eta=2a\cos\big(\tfrac12\Delta k\,x-\tfrac12\Delta\omega\,t\big)
   \cos(kx-\omega t)$ (7.66). The book prints ½Δω **x** in (7.66); it must be ½Δω **t** — the next line of the text and
   the figure use t, and with x the envelope could not move." equation (7.66) + `nb.code`: `b = ch07.beat_wave(x, t, 0.9,
   1.1)` and `ch07.beat_wave(x, t, 0.9, 1.1, printed=True)` at t = 0 and 20 s: the node of the correct envelope moves by
   c_g,finite·20 s, the printed one does not move. *expect:* `b["cg_finite"]` = 1.568 m/s vs dω/dk at k = 1: 1.566 m/s
   (deep, g = 9.81: c_g = ½√(g/k)), c = 3.132 m/s.
8. `nb.figure` — **Fig. 7.13 remake (N179)** (7 × 3 in): the beat pattern (blue), its envelope ±2a|cos(…)| (purple), a node
   (black dot) with a purple arrow c_g and a crest (orange dot) with an orange arrow c; the printed-slip envelope as a
   rose dashed ghost that stays put when time advances (shown at two times). *see:* "fast crests inside slow bulges." *read:*
   "nodes travel at Δω/Δk: no energy crosses a node, so the energy travels at the node speed." *change:* "…Δk halved: the
   groups become twice as long, the node speed tends to dω/dk."
9. `nb.figure` — **Fig. 7.14 remake (N180)** ω(k) deep water with, at k = 1 rad/m, the chord from the origin (orange, slope
   c = 3.13 m/s) and the tangent (purple, slope c_g = 1.57 m/s). *see:* "the tangent is flatter than the chord." *read:*
   "c_g < c wherever ω(k) bends down (gravity waves), c_g > c where it bends up (capillary waves)." *change:* "…shallow
   water: ω = k√(gH) is a straight line through the origin — chord and tangent coincide, c_g = c."
10. `nb.primer("Fourier integral and a packet's spectrum", "A real wave group is a continuous sum of sinusoids,
    $\\eta(x,t)=\\int A(k)e^{i(kx-\\omega(k)t)}dk$ (real part understood), where A(k) is its spectrum. A narrow spectrum of
    width δk around k₀ means a long group of length ~1/δk, and near k₀ we may Taylor-expand ω(k) (Ch. 1 P26). This is Ch. 5
    P142's discrete Fourier modes with the sum turned into an integral.", code="import numpy as np   # numbers\nx =
    np.linspace(-200, 200, 4001)   # [m]\nk = np.linspace(0.9, 1.1, 201); A = np.exp(-(k - 1)**2/(2*0.02**2))   # narrow
    spectrum\neta = (A[:, None]*np.cos(k[:, None]*x)).sum(0)*(k[1] - k[0])   # the integral as a sum\nprint(x[np.argmax(eta)],
    np.abs(eta[np.abs(x) > 150]).max()/eta.max())   # peak at 0; tails < 1 %: a group about 1/0.02 = 50 m long")`
    (**P172**).
11. `nb.note` — **N66 [B]** "A wave packet: all wavenumbers in a narrow band δk around k; in space a nearly sinusoidal
    wave whose amplitude dies away over a length ∝ 1/δk." + `nb.figure` **Fig. 7.15 remake (N181)** two panels: a Gaussian
    packet from `ch07.gaussian_packet` at t = 0 with its envelope, and |FFT| vs k (width δk marked). *see:* "a long packet
    has a narrow spectrum." *read:* "length × width ≈ constant (for a Gaussian, σ_x σ_k = 1)." *change:* "…the packet
    halved in length: its spectrum doubles in width — and it spreads faster (D20's dropped term)."
12. `nb.derivation("D20", …)` — Part F D20 (12 steps, ★★★, with `check_src`), ref "7.68".
13. `nb.note` — **N67 [B]** "For short times the packet keeps its shape and its envelope moves at c_g:
    $\eta=a(x-c_gt)\cos(kx-\omega t)$ (7.68) (the book cites Phillips; D20 writes it out). The nodes of the envelope move
    at c_g and no energy crosses them."
14. `nb.derivation("D21", …)` — Part F D21 (11 steps), ref "7.71".
15. `nb.note` — **N68 [B], N69 [B], N70 [B]** "$c_g=\frac c2\big[1+\frac{2kH}{\sinh(2kH)}\big]$ (7.69); $c_g=c/2$ (deep
    water) and $c_g=c$ (shallow water) (7.70); pure capillary waves have $c_g=3c/2$ (Exercise 7.9, checked by
    `group_velocity(g=0)`); and $F=E\frac c2\big[1+\frac{2kH}{\sinh(2kH)}\big]=Ec_g$ (7.71), E = ρga²/2 — **the rate of
    energy transmission is energy times group velocity**." **N71 [B]** "In 3-D, $c_{gi}=\partial\omega/\partial K_i$: the
    group velocity is the gradient of ω in wavenumber space (used for internal waves in C16, (7.143))." + `nb.code`:
    `print(ch07.group_velocity_vector(lambda K: np.sqrt(G*np.hypot(K[0], K[1])), np.array([0.6, 0.8])))` → c_g along K
    with size ½√(g/K) = 1.566 m/s: (0.940, 1.253) m/s.
16. `nb.worked_example("swell from a distant storm", "1. Deep water, λ = 156 m (T = 10 s): $c=\\sqrt{g\\lambda/2\\pi}=15.6$
    m/s. 2. $c_g=c/2=7.8$ m/s. 3. A storm 1000 km away: the first energy arrives after $10^6/7.8=128\\,000$ s ≈ 35.6 h —
    not 17.8 h, as the crest speed would suggest. 4. Longer swell (T = 20 s) arrives first: c_g = 15.6 m/s, 17.8 h — a
    forecaster reads the storm's distance from how fast the period falls.")`
17. `nb.code` — `k = 2*np.pi/156.0` · `print(ch07.phase_speed(k), ch07.group_velocity(k), ch07.group_velocity(k,
    30.0))` · `print(ch07.group_velocity_numeric(None, k), ch07.group_velocity(2*np.pi/0.005, sigma=0.0727, g=0.0)/
    ch07.phase_speed(2*np.pi/0.005, sigma=0.0727, g=0.0))` · `ps = ch07.packet_state(k, distance=1e6); print(ps)` ·
    `assert np.isclose(ch07.energy_flux(1.0, k), ch07.wave_energy_density(1.0)*ch07.group_velocity(k))`. *expect:* c =
    15.61, c_g = 7.80 m/s; the same k on 30 m (kH = 1.21): c = 14.27, c_g = 10.24 m/s (> the deep value: in intermediate depth c_g/c rises toward 1); numeric c_g
    equal to 1e-10; capillary ratio 1.500; `packet_state` with arrival_h = 35.6 h; the flux assertion passes (F = E c_g).
    *explain:* gloss "complex-step derivative: Im ω(k + ih)/h has no subtraction, so no cancellation error (Ch. 3 P107)".
18. `nb.primer("envelope with scipy.signal.hilbert", "The modulus of the 'analytic signal' `np.abs(scipy.signal.hilbert
    (eta))` follows the slowly varying amplitude of a fast oscillation — the envelope — without fitting. We use it to
    track a packet's peak.", code="import numpy as np, scipy.signal as ss   # numbers and signal tools\nx =
    np.linspace(0, 100, 2001); eta = np.exp(-((x - 50)/10)**2)*np.cos(3*x)   # a packet\nenv = np.abs(ss.hilbert(eta))   #
    its envelope\nprint(x[np.argmax(env)], env.max())   # 50.0, ≈ 1.0")` (**P173**).
19. `nb.check_agree` — **from scratch:** (1) `h = 1e-4; cg_mine = (ch07.omega_gravity(k + h) - ch07.omega_gravity(k -
    h))/(2*h)`; `assert np.allclose(cg_mine, ch07.group_velocity(k), rtol=1e-8)`. (2) Envelope-peak speed: evolve a
    Gaussian packet (k₀ = 1 rad/m, σ_x = 30 m, deep) with `ch07.linear_evolve` on a periodic 4096-point (FAST 1024) domain
    of 2000 m, find `x[np.argmax(ch07.envelope(eta))]` at 11 times up to 200 s, fit the slope; `assert np.allclose(slope,
    ch07.group_velocity(1.0), rtol=1e-2)`. *expect:* 1.566 m/s both ways (within 1 %). Gloss: "`np.fft.rfft` works on a
    periodic domain: make it long enough that nothing wraps round (Ch. 5 P142)."
20. `nb.animation` — **A2** (two parts in one cell, `frames=100 if not FAST else 50`): (a) the deep-water Gaussian packet:
    carrier (blue), envelope (purple), one crest tracked by an orange dot that is born at the rear, runs through the group
    and dies at the front (N73), with c and c_g speed bars; (b) caption card, then the stone in a pond (N72). Then *see /
    read / change*.
21. `nb.plotly` — **IF4:** `animate_figure` over t (0…200 s, 30 frames): the packet from `ch07.gaussian_packet(order=2)`
    with markers moving at c (orange) and at c_g (purple). Title "Crests at c, the group at c_g" (page-surviving A2).
22. `nb.note` — **N72 [B]** "**A stone in a pond** (Fig. 7.16): the impact contains all wavelengths; each travels at its
    own c_g, so the train sorts itself — longest gravity waves in front, and because capillary–gravity waves have a
    **minimum group velocity** of ≈ 17.8 cm/s (at λ ≈ 4.4 cm, Exercise 7.10), nothing is left inside a circle of radius
    ≈ 0.178 m/s × t: a calm centre. Heights fall as the train lengthens; viscosity finally damps each wave like
    $a_0e^{-2\nu k^2t}$ (Exercise 7.11, Ch. 8's tools)." + `nb.code`: `print(ch07.min_group_velocity(0.07274, 998.2))`;
    `eta = ch07.pond_ripples(x, t, width=0.01)` for t = 0.5, 1, 2 s on x ∈ [0, 1] m; `print(ch07.viscous_decay(1.0,
    2*np.pi/0.01, 1e-6, 10.0))`. *expect:* cg_min = 0.1776 m/s at λ = 4.35 cm; viscous decay of a 1 cm ripple over 10 s:
    a = e^{−2·1e-6·394784·10} = e^{−7.9} = 3.7×10⁻⁴ — short ripples die fast.
23. `nb.figure` — **Fig. 7.16 remake (N182)** η(x) from `pond_ripples` at three times (stacked), the front at
    c_g,max·t…(the long waves) and the calm-centre edge at c_g,min·t marked by purple dashed lines. *see:* "the train
    stretches; the region inside ~0.18t stays flat." *read:* "the longest waves lead; the calm edge moves at the minimum
    group velocity, not at c_min = 23 cm/s." *change:* "…σ = 0 (no surface tension): no minimum, the ripples fill the
    centre down to the shortest waves."
24. `nb.note` — **N73 [B]** "Following one crest of a deep-water group: it runs from the rear to the front at twice the
    group speed, its wavelength grows as it goes, and it dies at the front while new crests are born at the rear (seen in
    A2 and in E5)."
25. `nb.explainer("group_velocity_packets", heading="Where does a wave's energy go?", why="Follow one crest with a marker
    while the envelope moves at half its speed; switch the dispersion to shallow or capillary and see c_g = c or c_g > c;
    the chord and tangent of ω(k) move with the dot. Crests passing through a group can only be seen moving.",
    tries=["Play the *swell set*: count how many crests pass through one group.", "Choose *shallow*: the packet keeps
    its shape — no dispersion.", "Choose *capillary*: crests now appear at the front and vanish at the back.", "Open the
    Derivation (D20) and toggle the dropped ω″ term: the packet spreads."])`
26. `nb.md` — "> ⚠️ **Common confusion:** 'the wave speed' is two speeds. c = ω/k is how fast crests move; c_g = dω/dk is
    how fast the envelope, the energy and any signal move. For deep-water gravity waves c_g = c/2; for capillary waves
    c_g = 3c/2; only non-dispersive waves (shallow water, sound) have c_g = c."
27. `nb.md` — **What would change if…** "…the medium changed slowly along the way (depth decreasing toward a beach)? k
    and ω can no longer both stay constant: C10 shows that ω is the one that is carried unchanged along a path moving at
    c_g."
28. `nb.core("C10", "Kinematic wave theory: $\\frac{\\partial\\omega}{\\partial t}+c_g\\frac{\\partial\\omega}{\\partial x}=0$
    (7.79), rays and refraction", question="How do the wavenumber and the frequency of a wave train change when the depth
    changes slowly — and why do waves turn toward the shore?")`
29. `nb.md` — **The problem in plain words:** "Walk along a beach on any day: the breakers arrive nearly parallel to the
    shore, whatever direction the wind blew offshore. Waves round an island even reach its sheltered side. We want the
    rules that carry a wave train's wavelength, frequency and direction across slowly changing depth — the ray (WKB) theory
    that Ch. 13 uses for internal and planetary waves."
30. `nb.md` — **The idea** (ASCII):
    ```
    phase θ(x,t):  k = ∂θ/∂x (crests per metre),  ω = −∂θ/∂t (crests per second passing a point)
    crests are conserved  ⇒  ∂k/∂t + ∂ω/∂x = 0
    uniform depth:  k is carried at c_g (rays are straight lines in the x–t plane)
    depth H(x):     ω is carried at c_g, k changes;  the shoreward end of a crest slows  ⇒  the crest turns
    ```
31. `nb.note` — **N74 [B], N75 [B]** "A slowly varying train $\eta=a(x,t)\cos[\theta(x,t)]$ (7.72) has a local wavenumber
    and frequency $k(x,t)\equiv\frac{\partial}{\partial x}\theta(x,t)$ and $\omega(x,t)\equiv-\frac{\partial}{\partial t}
    \theta(x,t)$ (7.73) (for θ = kx − ωt they are the usual k and ω)." + `nb.code`: a chirped train `theta = lambda x, t:
    0.5*x + 0.002*x**2 - 2.0*t`; `print(ch07.local_wavenumber_frequency(theta, 10.0, 0.0))` → (0.54, 2.0).
32. `nb.primer("first-order wave equation and characteristics", "The equation ∂q/∂t + c ∂q/∂x = 0 says q does not change
    for an observer moving at speed c: along each line dx/dt = c in the x–t plane (a characteristic) q is constant. It is
    the material derivative of Ch. 3 with the velocity replaced by c.", code="import numpy as np   # numbers\nq0 = lambda
    x: np.exp(-x**2)   # a bump at t = 0\nc, t = 2.0, 3.0\nq = lambda x, t: q0(x - c*t)   # the solution: the bump shifted
    by ct\nprint(q(6.0, 3.0), q0(0.0))   # 1.0 1.0 — the value rode along x = ct")` (**P174**).
33. `nb.derivation("D22", …)` — Part F D22 (8 steps), ref "7.75". (Gloss in step 3's *why*: mixed partials commute,
    Ch. 4 P121.)
34. `nb.note` — **N76 [B], N77 [B]** "Crest conservation $\partial k/\partial t+\partial\omega/\partial x=0$ (7.74) and, with
    a dispersion relation, $\frac{\partial k}{\partial t}+c_g\frac{\partial k}{\partial x}=0$ (7.75): **wavenumbers are
    carried at the group velocity** — an observer moving at c_g always sees the same wavelength." + `nb.figure` **Fig.
    7.17 remake (N183)**: x–t diagram from the FFT packet of C09: crest lines (thin orange, slope c) and rays (thick
    purple, slope c_g), the packet's edges widening. *see:* "crest lines cross the rays and end at the front." *read:*
    "along a purple line k and ω stay fixed; along an orange line a crest keeps its phase." *change:* "…shallow water: both
    families have the same slope — crests and groups move together."
35. `nb.note` — **N78 [B], N79 [B], N80 [B]** "In slowly varying depth $\omega=\sqrt{gk\tanh[kH(x)]}$ has the form
    $\omega=\omega(k,x)$ (7.76), a local group velocity is $\partial\omega(k,x)/\partial k=c_g$ (7.77), and multiplying by
    ∂k/∂t gives $c_g\frac{\partial k}{\partial t}=\frac{\partial\omega}{\partial k}\frac{\partial k}{\partial t}=\frac
    {\partial\omega}{\partial t}$ (7.78)."
36. `nb.derivation("D23", …)` — Part F D23 (8 steps), ref "7.79".
37. `nb.primer("Snell's law for waves", "When a wave crosses a medium whose properties change only in one direction (say
    x), the wavenumber component along the other direction (y) cannot change: the crests must match along every line of
    constant x. With |k| changing, the angle α between k and the x-direction obeys |k| sin α = constant — the optics rule
    n₁ sin α₁ = n₂ sin α₂.", code="import numpy as np   # numbers\nk1, a1 = 0.07, np.radians(30)   # offshore wavenumber
    and angle\nk2 = 0.18   # a larger k closer to shore (shallower water)\nprint(np.degrees(np.arcsin(k1*np.sin(a1)/k2)))
    # 11.2° — the crest has turned toward the shore")` (**P175**).
38. `nb.derivation("D24", …)` — Part F D24 (9 steps), ref "7.79".
39. `nb.note` — **N47 [B]** "**Refraction on a sloping beach** (Fig. 7.8): ω is fixed along the path, $c=\sqrt{gH}$ and λ
    shrink as H falls, so $k\sin\alpha$ constant (D24) forces α toward zero — crests end up parallel to the depth
    contours." **N48 [B]** "Round a circular island with a gently sloping beach (Fig. 7.9) the same bending brings crests
    onto the shadow side too — optics' lens, made of water depth."
40. `nb.worked_example("swell turning on a 1:50 beach", "T = 8 s, arriving at 30° to the shore normal in 20 m of water.
    1. ω = 2π/8 = 0.785 rad/s (fixed along the ray). 2. At H = 20 m solve (7.28): k = 0.0708 rad/m (λ = 88.8 m). 3. At H = 2
    m: k = 0.181 rad/m (λ = 34.7 m). 4. Snell: $\\sin\\alpha_2=0.0708\\times\\sin30°/0.181=0.195$ → α₂ = 11.3°. 5. The crest
    turned by almost 19° while crossing 900 m of beach (from 20 m to 2 m depth on a 1:50 slope).")`
41. `nb.code` — `print(ch07.crest_conservation_residual(theta, 10.0, 0.0))` (a chirp satisfying (7.74): residual ~1e-9) ·
    `r = ch07.snell_ray_plane_beach(np.linspace(1000, 100, 10), np.radians(30), 1000.0, 8.0, 0.02)` · `print(np.degrees(
    r["alpha"]), r["snell"])` · `s = ch07.refraction_state(8.0, np.radians(30), 20.0, 2.0); print(s["alpha_deg"])` ·
    `ray = ch07.ray_trace(None, (1000.0, 0.0), (-0.0708*np.cos(np.radians(30)), 0.0708*np.sin(np.radians(30))), (0,
    250), H_fn=lambda X: 0.02*X[0])` · `print(ray["omega_drift"])`. *expect:* α falls from 30° at 20 m to 11.3° at 2 m; `snell`
    constant to 1e-10; α = 11.27° at H = 2 m; ω drift along the traced ray < 1e-8 relative.
42. `nb.check_agree` — **from scratch:** fixed-step RK4 (P95) on $dx/dt=\partial\omega/\partial k_x$, $dk/dt=-\partial
    \omega/\partial x$ with central differences for the plane beach; compare the angle at x = 100 m with `snell_ray_plane_
    beach` and ω along the path with its start: `assert np.allclose(alpha_mine, r_alpha, rtol=1e-6)`.
43. `nb.figure` — three panels (10 × 3.6 in): **(a) Fig. 7.8 remake (N174)** plan view of a 1:50 beach, depth contours
    (grey), 7 rays from offshore at 30° (purple) and crest lines drawn every wavelength (orange) turning parallel to the
    shore; **(b) Fig. 7.9 remake (N175)** rays round a circular island with H(r) rising linearly from 0 at the shore to 30
    m at r = 3 km (our bathymetry), some rays reaching the shadow side; **(c) Fig. 7.18 remake (N184)** x–t plot of four rays
    over the slope (curved, ω printed along each). Title "Frequency is carried along rays; the wavelength and direction
    change". *see:* "rays bending toward shallow water; curved rays in x–t." *read:* "each ray keeps its ω (printed values
    agree to 8 digits); k grows as the depth falls, so the ray bends toward the shore." *change:* "…a longer period
    (T = 14 s): the swell feels the bottom sooner, so the turning starts further offshore."
44. `nb.plotly` — **IF5:** slider over the incidence angle α₀ (0…60°, 13 steps): rays and crest lines over the plane beach
    from `snell_ray_plane_beach`.
45. `nb.explainer("wave_rays_refraction", heading="Why do waves arrive parallel to the beach?", why="Set the angle and the
    period and watch rays bend over a beach, a ridge or round an island; the side view prints ω flat along the selected
    ray while k, c and c_g change. The conservation is seen, not asserted.", tries=["Preset *swell at 30°*: read α at the
    2 m contour (11°).", "Switch to *island*: find a ray that reaches the shadow side.", "Choose *homogeneous*: the rays are
    straight — nothing to refract.", "Step the Derivation (D23) to 'multiply by c_g'."])`
46. `nb.md` — **What would change if…** "…the medium varied in time too (a tide changing the depth)? Then ω is no longer
    conserved along a ray: dω/dt = ∂ω/∂t at fixed k — the Doppler-like frequency shift used for waves on tidal currents
    and for internal waves in Ch. 13."

### A.6 §7.6 Nonlinear Waves in Shallow and Deep Water — R10, C11 (+N81–N85, N93–N97, N185, N186, N189 · D25 D26 · A5 · E7), C12 (+N86–N92, N187, N188 · D27 · A4 · E2)
1. `nb.section("7.6", "Nonlinear Waves in Shallow and Deep Water", intro="**What is this section about?** Everything so
   far assumed ka ≪ 1. Now the amplitude matters. In shallow water crests outrun troughs, fronts steepen and break into
   a hydraulic jump — unless dispersion balances the steepening and a solitary wave travels unchanged. In deep water the
   crests sharpen (Stokes waves) and the orbits of C05 fail to close: the water drifts forward (Stokes drift).")`
2. `nb.recap("R10", "The Froude number", "$\\mathrm{Fr}\\equiv u/\\sqrt{gH}=u/c$ (4.104) compares the flow speed with the
   shallow-water wave speed of C04 (7.49): Fr > 1 supercritical (waves cannot travel upstream), Fr < 1 subcritical.",
   where="Ch. 4 §4.11")` + `nb.code`: `print(similarity.froude_number(2.97, 0.1, g=9.81))` → 3.00.
3. `nb.core("C11", "The hydraulic jump $\\frac{H_2}{H_1}=\\frac12\\big(-1+\\sqrt{1+8\\mathrm{Fr}_1^2}\\big)$ (7.81) — and
   the other fate, the solitary wave", question="Given the incoming depth and speed, how high is the jump, where does the
   energy go — and when does a steep shallow-water wave not break at all?")`
4. `nb.md` — **The problem in plain words:** "At the foot of a dam spillway a thin, fast sheet of water suddenly rises into
   a deep, slow, churning roller; a tidal bore climbs up a river mouth as a moving step. Both are hydraulic jumps. We want
   the downstream depth from the upstream depth and speed, and the price paid in energy. Some steep waves, though, never
   break: a single hump can travel for kilometres unchanged (Russell's solitary wave, 1844)."
5. `nb.note` — **N81 [B]** "**Why shallow waves steepen.** Relative to water moving at u, a small wavelet travels at the
   local $c'=\sqrt{gH'}$; in the fixed frame at $c=c'+u$. Under a crest H′ and u are both larger, so crests outrun troughs
   and the front steepens (Fig. 7.19). The book stays qualitative; our labelled extension uses the exact simple wave of the
   shallow-water equations, $c=3\sqrt{g(H+\eta)}-2\sqrt{gH}$, which breaks at $t_b=-1/\min(\partial c/\partial x)$ at t = 0
   (characteristics: Ch. 15)." + `nb.code`: `x = np.linspace(0, 200, 2001); eta0 = 0.2*np.exp(-((x - 50)/10)**2)`; `s =
   ch07.simple_wave_evolve(eta0, x, np.linspace(0, 12, 7), 1.0)`; `print(s["t_break"])`; `print(ch07.nonlinear_wavelet_
   speed(np.array([-0.1, 0.0, 0.1]), 1.0, model="book"), ch07.nonlinear_wavelet_speed(np.array([-0.1, 0.0, 0.1]), 1.0))`.
   *expect:* t_break ≈ 13 s for this 20 cm hump on 1 m (printed); wavelet speeds ≈ [2.66, 3.13, 3.60] m/s (book form) and
   ≈ [2.65, 3.13, 3.59] (simple wave) — they agree to first order in η/H.
6. `nb.figure` — **Fig. 7.19 remake (N185)** frames of `simple_wave_evolve` at t = 0, 4, 8 s and at t_b (7 × 3 in), the
   front's slope growing; the characteristics (thin grey lines x = x₀ + c(η₀)t) crossing at t_b. *see:* "the front of the
   hump steepens and becomes vertical at t_b." *read:* "points higher up move faster: the curve would fold over after
   t_b — real water breaks instead." *change:* "…the hump halved: t_b roughly doubles (∝ 1/amplitude)."
7. `nb.note` — **N82 [B]** "After breaking, the front becomes a hydraulic jump: a nearly steady step joining a shallow
   supercritical stream (depth H₁, speed √(gH₁) < u₁) to a deeper subcritical one. It may stand still (spillway, a kitchen
   sink's circular jump) or move (a tidal bore). In the frame of a moving bore the same relations hold." + `nb.figure`
   **Fig. 7.20 remake (N186)**: `cv_box(ax, 0.1, 0.377)` — the channel, the jump roller, the dashed control volume with
   faces 1 and 2, u₁ and u₂ arrows, the hydrostatic face forces ½ρgH₁², ½ρgH₂² (amber arrows). *see:* "a CV around the
   jump with two vertical faces." *read:* "only four things cross the CV's accounts: two momentum fluxes and two pressure
   forces (the bed is frictionless, the roller is inside)." *change:* "…the observer ran upstream at the bore speed: a
   moving bore becomes this stationary jump (Galilean frames, Ch. 3)."
8. `nb.derivation("D25", …)` — Part F D25 (10 steps), ref "7.81". (Gloss in step 9's *why*: the quadratic formula; keep
   the positive root — Vieta, Ch. 3 P71: the roots multiply to −2Fr₁², so exactly one is positive.)
9. `nb.note` — **N83 [B]** "The momentum balance of D25 in the book's form: $Q^2\big(\frac1{H_2}-\frac1{H_1}\big)=\frac12g
   (H_1^2-H_2^2)$ (7.80)." equation (7.80).
10. `nb.derivation("D26", …)` — Part F D26 (8 steps), ref "7.81".
11. `nb.note` — **N84 [B]** "The mechanical energy per unit mass of a surface particle drops across the jump:
    $E_2-E_1=-(H_2-H_1)\frac{g(H_2-H_1)^2}{4H_1H_2}<0$ whenever H₂ > H₁; a 'jump down' (Fr₁ < 1) would create mechanical
    energy — the second law forbids it. The lost energy goes into the turbulent roller and ends as heat."
12. `nb.worked_example("a spillway jump", "H₁ = 0.1 m, Fr₁ = 3. 1. $u_1=\\mathrm{Fr}_1\\sqrt{gH_1}=3\\times0.990=2.97$ m/s,
    $Q=u_1H_1=0.297$ m²/s. 2. $H_2/H_1=\\tfrac12(-1+\\sqrt{1+72})=\\tfrac12(-1+8.544)=3.772$ → H₂ = 0.377 m. 3. $u_2=Q/H_2=
    0.788$ m/s, $\\mathrm{Fr}_2=0.788/\\sqrt{9.81\\times0.377}=0.41$ (subcritical). 4. Head loss $\\frac{(H_2-H_1)^3}{4H_1H_2}
    =\\frac{0.2772^3}{4\\times0.1\\times0.377}=0.141$ m — the energy of 14 cm of fall, ≈ ρgQ × 0.141 = 412 W per metre of
    width turned into heat.")`
13. `nb.code` — `j = ch07.hydraulic_jump(0.1, Fr1=3.0); print(j)` · `print(ch07.jump_momentum_residual(0.1, j["H2"],
    j["Q"]))` · `print(ch07.jump_state(0.1, 3.0))` · `print(ch07.hydraulic_jump(1.0, Fr1=0.7))` · `print(ch07.jump_state(1.0,
    3.0**0.5, frame="moving")["bore_speed"], ch04.bore_speed(1.0, 2.0))` (a bore of depth 2 m running into still water
    1 m deep, H₂/H₁ = 2 ⇔ Fr₁ = √3: the relative speed). *expect:* H₂ = 0.3772 m, u₂ = 0.788 m/s, Fr₂ = 0.410, dE =
    −1.385 J/kg, head loss 0.1412 m; residual ≈ 1e-16; momentum terms 882.9 and 234.1 N/m, pressure forces 49.05 and
    697.9 N/m (balance to 1e-12); Fr₁ = 0.7 gives H₂/H₁ = 0.609 with dE > 0 flagged `allowed=False`; bore speed 5.42 m/s both
    ways. (Builder: take the bore case from `ch04.bore_speed`'s own argument convention, stated in its docstring.)
14. `nb.check_agree` — **from scratch:** `r = np.roots([1.0, 1.0, -2*3.0**2])` (the quadratic of D25 step 8) · `r_pos =
    r[r > 0][0]` · `assert np.allclose(r_pos, ch07.hydraulic_jump(0.1, Fr1=3.0)["ratio"], rtol=1e-12)`. Markdown: "np.roots
    returns both roots, 3.772 and −4.772; physics keeps the positive one."
15. `nb.figure` — (7 × 3.4 in) twin axes: H₂/H₁ (black) and head loss/H₁ (rose) vs Fr₁ ∈ [0.3, 6]; the region Fr₁ < 1 greyed
    and labelled "would create energy (forbidden)"; dots at Fr₁ = 1.3 (undular jump), 3 (our example), 5 (spillway).
    Title "Momentum sets the height, energy sets the direction". *see:* "the depth ratio grows roughly like √2·Fr₁; the
    loss rises steeply." *read:* "at Fr₁ = 1 nothing happens (H₂ = H₁, no loss); below 1 the formula still gives a
    number but the loss would be negative." *change:* "…H₁ doubled at the same Fr₁: every ratio stays, the head loss
    doubles (it scales with H₁)."
16. `nb.note` — **N85 [C]** "In a dispersive medium the different wavelengths separate, which works against steepening; the
    two can balance and give waves of permanent form."
17. `nb.note` — **N93 [B]** "For waves with λ/H between about 10 and 20, Korteweg and de Vries (1895) showed
    $\frac{\partial\eta}{\partial t}+c_0\frac{\partial\eta}{\partial x}+\frac32c_0\frac\eta H\frac{\partial\eta}{\partial x}
    +\frac16c_0H^2\frac{\partial^3\eta}{\partial x^3}=0$, $c_0=\sqrt{gH}$ (7.87): term by term, shallow-water propagation,
    nonlinear steepening (N81), and the first correction for the depth not being quite shallow (weak dispersion)." equation
    (7.87). **N94 [B]** "Without the nonlinear term, η = a cos(kx − ωt) gives $c=c_0(1-\tfrac16k^2H^2)$ — the first two
    Taylor terms of $c=\sqrt{(g/k)\tanh kH}$ (7.29)." + `nb.figure` (5 × 3.4 in) log–log |c_KdV/c − 1| vs kH (0.01…1) from
    `kdv_linear_phase_speed` and `phase_speed`, with a slope-4 guide and the fitted slope printed (expect 4.00; the error
    is ≈ (19/360)(kH)⁴). *see/read/change* one line each.
18. `nb.note` — **N95 [B]** "Which fate wins? The ratio of the nonlinear to the dispersive term of (7.87) is
    $\frac{\eta}{H}\frac{\partial\eta}{\partial x}\Big/H^2\frac{\partial^3\eta}{\partial x^3}\sim\frac{a\lambda^2}{H^3}$ (the
    Ursell number; the book's text says '(7.88)' but means the terms of (7.87)). Above about 16 the front steepens into a
    jump; below, dispersion can balance it." + `nb.code`: `for a, lam in ((0.1, 10.0), (0.2, 10.0), (0.05, 30.0)):
    print(a, lam, ch07.ursell_number(a, lam, 1.0))` → 10, 20, 45.
19. `nb.note` — **N96 [C]** "The periodic permanent-form solutions of (7.87) are cnoidal waves (Jacobi's cn functions,
    drawn optionally with `ch07.cnoidal_wave`)." **N97 [B]** "The single-hump solution is the solitary wave
    $\eta=a\,\mathrm{sech}^2\big[\big(\frac{3a}{4H^3}\big)^{1/2}(x-ct)\big]$, $c=c_0\big(1+\frac a{2H}\big)$ (7.88): taller is
    faster (check by substitution, Exercise 7.15; gloss: $\tanh'=\mathrm{sech}^2$, $\mathrm{sech}'=-\mathrm{sech}\tanh$)."
    equation (7.88) + `nb.code`: `print(ch07.kdv_residual_sympy())` → 0 · `print(np.sqrt(G*1.0)*(1 + 0.2/2))` → 3.445 m/s
    for a = 0.2 m on H = 1 m.
20. `nb.figure` — **Fig. 7.23 remake (N189)** two panels: (a) a cnoidal wave train (`cnoidal_wave`, optional; else a
    solitary-wave train) and (b) `solitary_wave` for a = 0.2 m, H = 1 m with its speed; both on the bottom line. *see:*
    "sharp crests and long flat troughs; a lone hump." *read:* "steepening (the ηη_x term) and spreading (the η_xxx term)
    cancel exactly for this shape." *change:* "…a doubled: the hump narrows by √2 and speeds up to c₀(1 + a/H)."
21. `nb.animation` — **A5** (`player="frames"`; `frames=30 if not FAST else 16`): (a) the simple-wave hump steepening up to
    t_b (frames from N81's cell); (b) `ch07.kdv_solve` of a 0.1 m Gaussian hump of 10 m width on H = 1 m (N = 512, FAST 256,
    cached in `outputs/ch07/kdv_run.npz`) splitting into solitons, tallest (fastest) in front; `ch07.kdv_invariants` printed
    under the player (mass to 1e-10, energy to ~1e-6). *see / read / change* after it: "the same kind of hump either breaks
    (no dispersion) or sorts itself into solitons (dispersion on)."
22. `nb.explainer("hydraulic_jump", heading="How high does the water jump — and why only upward?", why="The Fr₁ slider moves
    the downstream depth, the momentum bars and the energy bar together; pushing Fr₁ below 1 turns the energy change
    positive — the second law choosing the root is felt, not just stated.", tries=["Preset *spillway*: read H₂ and the
    head loss.", "Drag Fr₁ to 0.7: the status turns red — the energy bar is positive.", "Switch to the moving-bore frame:
    the same depths, a different observer.", "Choose the solitary-wave fate and raise a/H: the Ursell number crosses 16."])`
23. `nb.md` — **What would change if…** "…the layer were a denser fluid under a lighter one (cold air under warm air,
    salty water under fresh)? g becomes the reduced gravity g′ of C14 and the same jump appears as an internal hydraulic
    jump in the lee of mountains (Ch. 13)."
24. `nb.core("C12", "Stokes waves and Stokes drift $\\bar u_L=a^2\\omega ke^{2kz_0}$ (7.85)", question="If every orbit is a
    closed circle, why does floating debris travel with the waves?")`
25. `nb.md` — **The problem in plain words:** "A swimmer floating beyond the surf slowly drifts toward the beach although
    'waves only move up and down'; oil, plastic and larvae at sea are carried downwind by the waves themselves. We want
    that mean drift: its size, how fast it dies with depth, and why a current meter fixed in place measures no mean
    current at all."
26. `nb.note` — **N86 [B], N87 [B]** "In 1847 Stokes found periodic waves of finite amplitude in deep water:
    $\eta=a\cos k(x-ct)+\tfrac12ka^2\cos2k(x-ct)+\tfrac38k^2a^3\cos3k(x-ct)+\ldots$ (7.82) with
    $c=\sqrt{\frac gk\big(1+k^2a^2+\ldots\big)}$ (7.83): the harmonics sharpen the crests and flatten the troughs, and the
    speed grows with amplitude (4.4 % at ka = 0.3). Stated here; a sympy cell checks the second-harmonic coefficient ½
    and the O(k²a²) speed correction (`ch07.stokes_expansion_sympy()` → alpha = 1/2, gamma = 1)." + `nb.code` (check
    cell): `print(ch07.stokes_expansion_sympy())`; `print(ch07.stokes_wave_speed(1.0, 0.3)/ch07.phase_speed(1.0))` → 1.0440.
27. `nb.note` — **N88 [B]** "The steepest possible Stokes wave has a crest angle of 120° and height/wavelength 0.1411
    (`ch07.STOKES_LIMIT_STEEPNESS`, Dyachenko et al. 2016) — amplitude a ≈ 0.07λ; beyond it the crest spills (white caps).
    The truncated series (7.82) cannot reach that corner; we show profiles only up to ka = 0.3."
28. `nb.figure` — **Fig. 7.21 remake (N187)** `stokes_wave_profile` at ka = 0.3 (black) vs the linear cosine (blue dashed)
    over two wavelengths. *see:* "higher, sharper crests; shallower, broader troughs." *read:* "the second harmonic adds at
    the crests and subtracts at the troughs." *change:* "…ka = 0.1: the two curves almost coincide (the correction is ∝ ka)."
29. `nb.note` — **N89 [B]** "**Stokes drift.** At finite amplitude a particle is a little further forward at the top of its
    loop, where the forward velocity is larger, than it is backward at the bottom, where it is smaller: the loop does not
    close and the particle creeps forward (Fig. 7.22)." Gloss (in the markdown): "**Lagrangian vs Eulerian mean** — the
    average of the velocity *following a particle* (ū_L) versus the average *at a fixed point* (ū): they differ whenever the
    particle samples a non-uniform field (Ch. 3 §3.2)."
30. `nb.note` — **N90 [B]** "Keep the first Taylor terms of u about the mean position (Ch. 3 P98):
    $\frac{dx_p}{dt}=u(x_0,z_0,t)+\xi\big(\frac{\partial u}{\partial x}\big)_{x_0,z_0}+\zeta\big(\frac{\partial u}{\partial z}
    \big)_{x_0,z_0}+\ldots$ (7.84a) and the same for w (7.84b)." equations (7.84a), (7.84b).
31. `nb.derivation("D27", …)` — Part F D27 (9 steps), ref "7.85".
32. `nb.note` — **N91 [B]** "At any depth the same moves with cosh and sinh give
    $\bar u_L=a^2\omega k\frac{\cosh(2k(z_0+H))}{2\sinh^2(kH)}$ (7.86) (Exercise 7.14); the vertical drift is zero
    (Exercise 7.13)." equation (7.86). **N92 [B]** "At a point that stays under water the mean velocity is exactly zero:
    irrotationality gives $u=u\vert_{z=-H}+\int_{-H}^z\frac{\partial w}{\partial x}dz$, whose wavelength average vanishes for a
    periodic wave. All the transport is Lagrangian." + `nb.code`: `print(ch07.eulerian_mean_u(-0.5, 0.05, 1.0, 5.0),
    ch07.stokes_drift(-0.5, 0.05, 1.0, 5.0))` → 0 (to 1e-12) vs a positive drift.
33. `nb.worked_example("drift under a 1 m swell", "a = 1 m, T = 8 s, deep water. 1. ω = 0.785 rad/s, k = ω²/g = 0.0629
    rad/m. 2. Surface drift $a^2\\omega k=1\\times0.785\\times0.0629=0.0494$ m/s ≈ 4.9 cm/s ≈ 4 km per day. 3. It decays like
    $e^{2kz_0}$ — twice as fast as the orbits: halved every $\\ln2/2k=5.5$ m. 4. Compare the orbital speed aω = 0.785 m/s:
    the drift is ka = 6 % of it — a second-order effect.")`
34. `nb.code` — `k = 1.0; a = 0.05; w = ch07.omega_gravity(k)` · `for z0 in (0.0, -0.5, -1.0): print(z0,
    ch07.stokes_drift(z0, a, k), ch07.stokes_drift_numeric(z0, a, k, periods=20 if not FAST else 10))` · `print(ch07.
    stokes_drift(0.0, 1.0, (2*np.pi/8)**2/G))`. *expect:* analytic 7.83e-3, 2.88e-3, 1.06e-3 m/s at the three depths;
    numeric within ~1 % (differences O(ka)²); 0.0494 m/s for the swell. *explain:* gloss "`solve_ivp` with DOP853 and
    rtol = 1e-10: the drift is a few % of the orbit, so the integrator must be far more accurate than that (Ch. 3 P94)".
35. `nb.check_agree` — **from scratch:** integrate the exact path line (7.33) over 20 periods from (0, −0.5) with
    `ch07.particle_path(model="exact")`, fit the mean slope of x_p(t) (period-averaged positions), `assert
    np.allclose(slope, ch07.stokes_drift(-0.5, 0.05, 1.0), rtol=2e-2)`. Markdown: "The particles really do creep forward at
    (7.85) — the analytic drift is the average of the exact motion, not an extra assumption."
36. `nb.figure` — (6 × 3.6 in) ū_L(z₀) for deep water (analytic (7.85), purple) and for kH = 1 ((7.86), purple dashed), with
    numeric dots from `stokes_drift_numeric` at 6 depths; the orbit radius ae^{kz₀} (teal, scaled) for comparison of decay
    rates. Title "The drift lives close to the surface". *see:* "the purple curve falls twice as fast as the teal one." *read:*
    "e^{2kz₀} vs e^{kz₀}: at one wavelength's depth (z₀ = −λ) the orbit keeps e^{−2π} ≈ 0.2 % and the drift e^{−4π} ≈ 4×10⁻⁶."
    *change:* "…finite depth (kH = 1): the drift does not vanish at the bottom — cosh 2k(z₀ + H) ≥ 1."
37. `nb.animation` — **A4 (N188)** exact path lines of a steep wave (ka = 0.25, deep): 9 particles starting on a vertical
    line (`ch07.dyed_line`), their open loops traced, the dyed line redrawn each period leaning forward
    (`frames=80 if not FAST else 40`). *see / read / change:* "the top of the dye line runs ahead; after a few periods the
    line is curved like Fig. 7.22 — the shape of e^{2kz₀}."
38. `nb.explainer("particle_orbits", heading="What does the water under a wave actually do?", why="Tracer particles loop
    under the moving surface; drag the depth to squash circles into ellipses, then switch the model from linear to exact
    and watch the loops open and the dyed column lean forward. Motion and one model toggle are the whole lesson.",
    tries=["Start in *deep water*, click a particle at 1 m depth and read its A, B in the inspector.", "Drag kH to 0.3:
    the circles become flat ellipses; the bottom particle only slides.", "Choose *steep deep wave (exact)*: after five
    periods read the drift in the end-of-run card and compare with a²ωk.", "Open the Derivation (D27) at 'average the
    two terms'."])`
39. `nb.md` — "> ⚠️ **Common confusion:** 'the mean current is zero' and 'the water drifts' are both true. A fixed meter
    averages over whatever water passes it (zero); a float averages following one parcel (a²ωk e^{2kz₀}). Ch. 13 adds the
    Earth's rotation, and the Stokes drift then drives the Stokes–Coriolis force and Langmuir circulation in the ocean's
    surface layer."
40. `nb.md` — **What would change if…** "…the waves broke? Breaking transfers momentum to the water as a real (Eulerian)
    current — the longshore currents and undertow of beaches — which the irrotational theory here cannot describe."

### A.7 §7.7 Waves on a Density Interface — R11, C13 (+N98–N108, N190–N192 · D28 · P176), C14 (+N109–N131, N193, N194 · D29 D30 D31 · IF6 · E8)
1. `nb.section("7.7", "Waves on a Density Interface", intro="**What is this section about?** Replace the air above the
   water by a slightly lighter liquid — fresh water over salt water, warm water over cold. Gravity still restores the
   interface, but only through the small density difference: the waves become slow and tall. With a free surface on
   top, a layer of fluid can oscillate in two ways — the surface way and the interface way — the barotropic and baroclinic
   modes that layered ocean and atmosphere models are built on.")`
2. `nb.recap("R11", "Laplace in each fluid", "Each fluid is irrotational and incompressible on its own side, so
   $\\frac{\\partial^2\\phi_1}{\\partial x^2}+\\frac{\\partial^2\\phi_1}{\\partial z^2}=0$ and the same for φ₂ (7.90): one
   Laplace equation per layer.", where="Ch. 6 §6.2")`
3. `nb.core("C13", "Interfacial waves $\\omega=\\sqrt{gk\\frac{\\rho_2-\\rho_1}{\\rho_2+\\rho_1}}=\\varepsilon\\sqrt{gk}$
   (7.95)", question="How does a wave on a weak density step differ from a wave on the sea surface?")`
4. `nb.note` — **N98 [C]** "Set-up: a lighter fluid ρ₁ above a heavier one ρ₂, both infinitely deep, small slopes, no
   interfacial tension (Fig. 7.24) — a fjord's fresh upper layer, a sun-warmed surface layer, oil on water."
5. `nb.md` — **The problem in plain words:** "In a fjord, a layer of river water floats on sea water; in summer a lake's
   warm surface layer floats on cold water. Their interface carries waves you cannot see from a boat but a thermometer
   string can: slow, and metres tall. We want their speed and why they are so different from surface waves."
6. `nb.md` — **The idea**: "Gravity pulls a raised piece of interface down with the weight of the *difference* ρ₂ − ρ₁,
   but the inertia to be moved is both fluids, ρ₂ + ρ₁:
   ```
   surface wave:    ω² = gk · (ρ − 0)/(ρ + 0)            = gk
   interface wave:  ω² = gk · (ρ₂ − ρ₁)/(ρ₂ + ρ₁)  = ε² gk    ε ≈ 0.03 for Δρ/ρ = 0.002
   ```
   **A small density step makes a slow wave — and the two fluids slide past each other.**"
7. `nb.primer("complex amplitudes", "Write a real wave as the real part of a complex one, $\\zeta=\\mathrm{Re}\\{a\\,e^{i(kx-
   \\omega t)}\\}$, and drop the Re during linear algebra: ∂/∂x becomes multiplication by ik and ∂/∂t by −iω, and a
   complex amplitude b = |b|e^{iφ} carries both size and phase shift (Ch. 1 P45 Euler's formula; Ch. 6 P153). ⚠️ Take real
   parts before multiplying two fields (energy, flux: C16).", code="import numpy as np   # numbers\nk, w, t = 1.0, 2.0,
   0.3\nx = np.linspace(0, 1, 3)\nzeta = np.real(1.0*np.exp(1j*(k*x - w*t)))   # a cos(kx − ωt) with a = 1\ndzdt =
   np.real(-1j*w*1.0*np.exp(1j*(k*x - w*t)))   # ∂/∂t ↦ −iω\nprint(np.allclose(dzdt, w*np.sin(k*x - w*t)))   # True")`
   (**P176**).
8. `nb.note` — **N99 [B]** "From here on the book writes $\zeta(x,t)=a\exp[i(kx-\omega t)]$ (7.89) (N100) with Re{}
   dropped; complex constants A, B, C, b carry phase shifts. `W.real_field(amp, phase)` returns the real part when we plot."
   **N100 [B]** equation (7.89).
9. `nb.note` — **N101 [B], N102 [B], N103 [B], N104 [B]** "The problem: (7.90) in each fluid; $\phi_1\to0$ as $z\to\infty$
   (7.91); $\phi_2\to0$ as $z\to-\infty$ (7.92); both fluids move with the interface, $\frac{\partial\phi_1}{\partial z}=
   \frac{\partial\phi_2}{\partial z}=\frac{\partial\zeta}{\partial t}$ at z = 0 (7.93); and the pressure is continuous (no
   interfacial tension), $\rho_1\frac{\partial\phi_1}{\partial t}+\rho_1g\zeta=\rho_2\frac{\partial\phi_2}{\partial t}+
   \rho_2g\zeta$ at z = 0 (7.94) — C02's conditions, now two-sided." equations (7.91)–(7.94).
10. `nb.derivation("D28", …)` — Part F D28 (10 steps), ref "7.95".
11. `nb.worked_example("a thermocline wave", "ρ₁ = 1000, ρ₂ = 1002 kg/m³, λ = 100 m. 1. $\\varepsilon^2=2/2002=
    9.99\\times10^{-4}$, ε = 0.0316. 2. A surface wave of this length: $\\omega=\\sqrt{gk}=0.785$ rad/s, T = 8.0 s. 3. The
    interfacial wave: ω = 0.0316 × 0.785 = 0.0248 rad/s, T = 253 s ≈ 4.2 min — 32 times slower. 4. Same energy as a 1 m
    surface wave needs $a=\\sqrt{\\rho/\\Delta\\rho}\\times1\\,\\text{m}=\\sqrt{500}=22$ m.")`
12. `nb.code` — `k = 2*np.pi/100` · `print(ch07.interface_omega(k, 1000.0, 1002.0), ch07.omega_gravity(k))` · `f =
    ch07.interface_fields(np.array([0.0, 25.0]), np.array([0.0, 0.0]), 0.0, 1.0, k, 1000.0, 1002.0)`; `print(f["u1"],
    f["u2"], f["gamma_sheet"])` · `print(ch07.interface_residuals(10.0, 3.0, 1.0, k, 1000.0, 1002.0))` · `print(ch07.
    interface_energy(1.0, k, 1000.0, 1002.0), ch07.interface_energy(1.0, k, 1000.0, 1002.0, method="quad"))` · `print(
    ch07.interface_omega(k, 1002.0, 1000.0))` (heavier on top). *expect:* 0.02481 vs 0.7851 rad/s; u₁ = −u₂ = −0.0248 m/s at
    the crest (x = 0), sheet strength 0.0496 m/s; residuals of (7.90)–(7.94) ≲ 1e-10; E_k = E_p = 4.905 J/m², E = 9.81 J/m²
    (¼Δρga² each) both ways; NaN with a warning — the unstable Rayleigh–Taylor case of Ch. 11.
13. `nb.check_agree` — **from scratch:** solve the two kinematic conditions (7.93) for A and B with `np.linalg.solve([[-k,
    0], [0, k]], [-1j*w*a, -1j*w*a])`, then put them into (7.94) and solve for ω²: `assert np.allclose(w_mine,
    ch07.interface_omega(k, 1000.0, 1002.0), rtol=1e-12)` (Ch. 2 P57 `np.linalg.solve` works with complex numbers too).
14. `nb.note` — **N105 [B], N106 [B]** "The kinetic and potential energies are equal, $E_k=E_p=\tfrac14(\rho_2-\rho_1)ga^2$
    (Exercise 7.18; the column swap of Fig. 7.25 now lifts ρ₂ and lowers ρ₁), so $E=E_k+E_p=\tfrac12(\rho_2-\rho_1)ga^2$
    (7.96): for the same energy an interfacial wave is $\sqrt{\rho_2/(\rho_2-\rho_1)}$ times taller than a surface wave
    (22× for the thermocline above)." equation (7.96). **N191 [C]** "Fig. 7.25 is the column-swap sketch of C06 with a
    second fluid (not redrawn)."
15. `nb.note` — **N107 [B]** "The horizontal velocities $u_1=\partial\phi_1/\partial x=-\omega ae^{-kz}e^{i(kx-\omega t)}$ and
    $u_2=\omega ae^{kz}e^{i(kx-\omega t)}$ are opposite: the interface is a **vortex sheet** of strength $u_2-u_1=2\omega a
    \cos(kx-\omega t)$ (the jump of Ch. 5's vortex sheets) — the seed of the Kelvin–Helmholtz instability (Ch. 11). A
    continuous stratification spreads this vorticity through the fluid, so internal waves (C15) are not irrotational." +
    `nb.code`: `print(ch05.vortex_sheet_strength(f["u2"][0], f["u1"][0]))` = `f["gamma_sheet"][0]` (sign convention stated
    in the ch05 docstring).
16. `nb.figure` — **Fig. 7.24 remake (N190)** (7 × 3 in): the interface ζ (navy), arrows of u₁ above (light) and u₂ below
    (navy) decaying away from the interface, the sheet marked with ± signs of its strength. *see:* "above and below the
    interface the water moves in opposite directions." *read:* "the arrows' size ∝ e^{−k|z|}: the motion is confined within
    about a wavelength of the interface." *change:* "…ρ₁ → 0 (air over water): the upper arrows vanish and (7.95) returns
    $\omega=\sqrt{gk}$, the deep-water wave of C04."
17. `nb.note` — **N108 [C], N192 [C]** "'Dead water' (Fig. 7.26): a ship in a fjord with a thin fresh layer makes
    interfacial waves and feels an unexpected drag — Bjerknes' explanation of the old sailors' puzzle (two-layer flows
    return in Ch. 13)."
18. `nb.md` — **What would change if…** "…the upper layer had a free surface a finite distance above? A second wave
    (on the surface) enters and the two interact — C14's barotropic and baroclinic modes."
19. `nb.core("C14", "Barotropic and baroclinic modes; reduced gravity $c=\\sqrt{g'H}$, $g'=g\\frac{\\rho_2-\\rho_1}{\\rho_2}$
    (7.116)–(7.117)", question="What are the two ways a layer over deep water can oscillate — and why is one of them almost
    invisible at the surface?")`
20. `nb.note` — **N109 [C]** "Set-up (Fig. 7.27): an upper layer of thickness H with a free surface, over an infinitely
    deep lower layer; the origin at the mean free surface, the interface at z = −H. Two modes: surface and interface in
    phase or in antiphase."
21. `nb.md` — **The problem in plain words:** "The summer thermocline of a lake sits 10–50 m down; tides and winds move both
    the surface and the thermocline. Ocean models are built from exactly two kinds of motion: the barotropic one, which
    moves the whole column and shows at the surface, and the baroclinic one, which heaves the thermocline by metres while
    the surface barely twitches. We want both from one calculation."
22. `nb.md` — **The idea** (ASCII):
    ```
    barotropic:  surface ~~~   interface ~ (in phase, smaller by e^{−kH})     ω² = gk          fast
    baroclinic:  surface ~ (tiny, antiphase)   interface ~~~~~ (large)        ω² ≈ g′kH·k      slow
    ```
23. `nb.note` — **N110–N118 [B]** "The problem: $\phi_2\to0$ at $z\to-\infty$ (7.97); $\frac{\partial\phi_1}{\partial z}=
    \frac{\partial\eta}{\partial t}$ at z = 0 (7.98) (the book prints ∂φ₁/dz; read ∂z); $\frac{\partial\phi_1}{\partial t}+g\eta
    =0$ at z = 0 (7.99); $\frac{\partial\phi_1}{\partial z}=\frac{\partial\phi_2}{\partial z}=\frac{\partial\zeta}{\partial t}$
    at z = −H (7.100); $\rho_1\frac{\partial\phi_1}{\partial t}+\rho_1g\zeta=\rho_2\frac{\partial\phi_2}{\partial t}+\rho_2g
    \zeta$ at z = −H (7.101); with $\eta=ae^{i(kx-\omega t)}$ (7.102) (a real), $\zeta=be^{i(kx-\omega t)}$ (7.103) (b complex:
    a phase difference is allowed), $\phi_1=(Ae^{kz}+Be^{-kz})e^{i(kx-\omega t)}$ (7.104) and $\phi_2=Ce^{kz}e^{i(kx-\omega
    t)}$ (7.105). The book prints (7.105) with $e^{i(kz-\omega t)}$; that cannot satisfy (7.100) for all x — read
    $e^{i(kx-\omega t)}$ (`ch07.two_layer_sympy()["printed_7105_residual"]` ≠ 0 shows it)."
24. `nb.derivation("D29", …)` — Part F D29 (12 steps, ★★★, with `check_src`), ref "7.109".
25. `nb.note` — **N119–N122 [B]** "The constants: $A=-\frac{ia}2\big(\frac\omega k+\frac g\omega\big)$ (7.106),
    $B=\frac{ia}2\big(\frac\omega k-\frac g\omega\big)$ (7.107), $C=-\frac{ia}2\big(\frac\omega k+\frac g\omega\big)-\frac{ia}2
    \big(\frac\omega k-\frac g\omega\big)e^{2kH}$ (7.108) and $b=\frac a2\big(1+\frac{gk}{\omega^2}\big)e^{-kH}+\frac a2\big(1-
    \frac{gk}{\omega^2}\big)e^{kH}$ (7.109)." equations.
26. `nb.derivation("D30", …)` — Part F D30 (12 steps, ★★★, with `check_src`), ref "7.110".
27. `nb.note` — **N123 [B]** "$\big(\frac{\omega^2}{gk}-1\big)\big\{\frac{\omega^2}{gk}[\rho_1\sinh kH+\rho_2\cosh kH]-(\rho_2-
    \rho_1)\sinh kH\big\}=0$ (7.110) (Exercise 7.19): two factors, two modes."
28. `nb.derivation("D31", …)` — Part F D31 (11 steps), ref "7.117".
29. `nb.note` — **N124 [B], N125 [B]** "First root $\omega^2=gk$ (7.111) — a deep-water surface wave — with $b=ae^{-kH}$
    (7.112): the interface moves in phase, smaller by $e^{-kH}$: the **barotropic** mode (surfaces of constant pressure
    and density coincide)." **C14's second root** "$\omega^2=\frac{gk(\rho_2-\rho_1)\sinh kH}{\rho_2\cosh kH+\rho_1\sinh kH}$
    (7.113), which becomes (7.95) as kH → ∞." **N126 [B]** "$\eta=-\zeta\big(\frac{\rho_2-\rho_1}{\rho_1}\big)e^{-kH}$
    (7.114): surface and interface in antiphase and the surface tiny — the **baroclinic** (internal) mode; u reverses
    across the interface." **N127 [B]** "Long waves (kH ≪ 1): $\omega^2=kg\big(\frac{\rho_2-\rho_1}{\rho_2}\big)kH$ (7.115)."
    **N128 [B]** "so $c=[g'H]^{1/2}$, $g'=g\big(\frac{\rho_2-\rho_1}{\rho_2}\big)$ (7.116, 7.117): the shallow-water speed
    with g replaced by the **reduced gravity**." **N129 [B]** "$\eta=-\zeta\big(\frac{\rho_2-\rho_1}{\rho_1}\big)$ (7.118)."
    **N130 [B]** "$p'=-\rho_1\frac{\partial\phi_1}{\partial t}=i\rho_1\omega(A+B)e^{i(kx-\omega t)}=\rho_1g\eta$ (7.119): the
    upper layer is hydrostatic for long waves — shallow means hydrostatic again (C04)."
30. `nb.md` — "> ⚠️ **Which ρ in g′?** The book's (7.117) divides by the **lower** density ρ₂; Ch. 4's `reduced_gravity`
    (P131) divided by ρ₁. Ocean (1025 over 1027 kg/m³): 0.01910 vs 0.01914 m/s² (0.2 % apart); oil (800) over water (1000):
    1.96 vs 2.45 m/s² (25 %). Our code: `W.reduced_gravity_book(rho1, rho2, ref="lower")` for the book, `ref="upper"` for
    Ch. 4's. Ch. 13 uses a reference ρ₀ — for the ocean all three agree to 0.2 %." + `nb.code`: print both for the two
    cases.
31. `nb.worked_example("a thermocline 50 m down", "H = 50 m, ρ₁ = 1000, ρ₂ = 1002 kg/m³. 1. $g'=9.81\\times2/1002=0.0196$
    m/s². 2. $c=\\sqrt{g'H}=\\sqrt{0.98}=0.99$ m/s — against $\\sqrt{gH}=22.1$ m/s for the barotropic long wave. 3. A 10 m
    heave of the thermocline shows at the surface as $\\eta=-10\\times2/1000=-0.02$ m: 2 cm, in antiphase. 4. For λ = 1 km
    (kH = 0.31, not quite long) (7.113) gives T = 19.5 min, c = 0.85 m/s.")`
32. `nb.code` — `for lam in (100.0, 1000.0, 5000.0): k = 2*np.pi/lam; print(lam, ch07.two_layer_free_surface_omega(k, 50.0,
    1000.0, 1002.0), ch07.two_layer_state(k, 50.0, 1000.0, 1002.0))` · `print(ch07.two_layer_long_wave_speed(50.0, 1000.0,
    1002.0))` · `m = ch07.two_layer_modes(2*np.pi/1000, 50.0, 1000.0, 1002.0, mode="baroclinic"); print(m["eta_over_zeta"])` ·
    `print(ch07.two_layer_rigid_lid_omega(2*np.pi/1000, 50.0, 200.0, 1000.0, 1002.0))` (N131, Exercise 7.20 form).
    *expect:* baroclinic periods 4.22, 19.5, 86.9 min, barotropic 8.0, 25.3, 56.6 s; c_long = 0.990 m/s; η/ζ = −0.00146 at
    λ = 1 km (→ −0.002 for long waves).
33. `nb.check_agree` — **from scratch:** (7.110) is a quadratic in s = ω²/gk: expand the product to
    $(\rho_1\sinh kH+\rho_2\cosh kH)s^2-[\rho_1\sinh kH+\rho_2\cosh kH+(\rho_2-\rho_1)\sinh kH]s+(\rho_2-\rho_1)\sinh kH=0$,
    `np.roots` it, take ω = √(gk s) for both roots: `assert np.allclose(sorted(w_mine), sorted(ch07.
    two_layer_free_surface_omega(k, 50.0, 1000.0, 1002.0)), rtol=1e-10)`.
34. `nb.figure` — **Fig. 7.27 remake (N193)** two panels from `two_layer_modes` at λ = 500 m (amplitudes exaggerated, true
    ratios printed): barotropic (surface and interface in phase, u arrows the same way in both layers) and baroclinic
    (tiny surface in antiphase, big interface, u opposite above and below). *see:* "in phase vs antiphase." *read:* "the
    printed η/ζ: +e^{kH} ≈ 1.9 for the barotropic mode (the interface is the smaller one), −0.0011 for the baroclinic mode."
    *change:* "…Δρ → 0: the baroclinic surface signal vanishes and its speed goes to zero like √Δρ."
35. `nb.figure` — **Fig. 7.28 remake (N194)** + **N131 [B]** "When both layers are thin compared with λ and Δρ/ρ ≪ 1
    (Boussinesq), the barotropic mode has the same u at every depth and the baroclinic mode has uniform, opposite u in the
    two layers (two layers between rigid lids: Exercise 7.20)." u(z) profiles from `two_layer_modes(long_wave=True)`.
36. `nb.figure` — ω(k) of both modes (log–log, k ∈ [10⁻⁴, 1] rad/m, H = 50 m, Δρ = 2 kg/m³) with the limits √(gk) and
    ε√(gk) (7.95) (dashed) and √(g′H)·k (dotted). *see:* "two branches far apart; the baroclinic one joins ε√(gk) for short
    waves." *read:* "for kH ≫ 1 the surface is too far away to matter — the interface forgets it (7.113) → (7.95)."
    *change:* "…H doubled: the baroclinic long-wave line moves up by √2."
37. `nb.plotly` — **IF6:** `slider_figure` over Δρ/ρ₂ (log, 1e-4…0.5, 25 steps): ω(k) of both modes with the (7.95) limit.
38. `nb.explainer("two_layer_modes", heading="What is a baroclinic mode?", why="The density-contrast and thickness sliders
    move both dispersion curves and both mode shapes at once; the reader sees the interface amplitude dwarf the surface one
    as Δρ → 0 and recovers the two-deep-fluids case as kH → ∞.", tries=["Preset *ocean thermocline*: compare the two
    periods in the Explain tab.", "Switch between the modes: watch u reverse across the interface in the baroclinic one.",
    "Drag Δρ/ρ up to 0.2 (oil over water) and read how different the two g′ conventions are.", "Step the Derivation (D30) to
    'factor out (s − 1)'."])`
39. `nb.md` — **What would change if…** "…the Earth rotated? The baroclinic long-wave speed √(g′H) ≈ 1 m/s divided by the
    Coriolis parameter f ≈ 10⁻⁴ s⁻¹ gives the internal (baroclinic) Rossby radius ≈ 10 km — the size of ocean eddies
    (Ch. 13)."

### A.8 §7.8 Internal Waves in a Continuously Stratified Fluid — R12–R20, C15 (+N132–N142, N145–N147, N196 · D32 D33 D34 · P177 · IF8), C16 (+N143, N144, N148–N166, N195, N197–N199 · D35 D36 D37 · P178 · A6 · E9), S01
1. `nb.section("7.8", "Internal Waves in a Continuously Stratified Fluid", intro="**What is this section about?** Now the
   density changes smoothly with height, as in most of the ocean and atmosphere. Displaced parcels bob at the buoyancy
   frequency N, and waves can travel in any direction — with astonishing rules: the frequency depends only on the direction
   of the wave, never on its wavelength; the water moves along the crests; and the energy travels at right angles to the
   crests.")`
2. `nb.recap("R12", "Assumptions", "Boussinesq (density differences matter only in the weight), inviscid, small
   amplitude (advection dropped), and frequencies well above the Coriolis frequency (rotation joins in Ch. 13).",
   where="Ch. 4 §4.9")`
3. `nb.recap("R13", "Density and continuity equations", "$D\\rho/Dt=0$ (4.9) (each parcel keeps its density) and
   $\\partial u/\\partial x+\\partial v/\\partial y+\\partial w/\\partial z=0$ (4.10) (incompressible flow).", where="Ch. 4
   §4.2")` + `nb.code`: `print(ch07.internal_wave_fields(0.3, -0.2, 1.0, 1.0, 1.0, 0.01, 1e-3, residuals=True)[
   "continuity"])` → ≈ 1e-12 (the wave of C16 is divergence-free).
4. `nb.recap("R14", "Linear Boussinesq momentum", "$\\frac{\\partial u}{\\partial t}=-\\frac1{\\rho_0}\\frac{\\partial p}
   {\\partial x}$, $\\frac{\\partial v}{\\partial t}=-\\frac1{\\rho_0}\\frac{\\partial p}{\\partial y}$, $\\frac{\\partial w}
   {\\partial t}=-\\frac1{\\rho_0}\\frac{\\partial p}{\\partial z}-\\frac{\\rho g}{\\rho_0}$ (7.120, 7.121, 7.122): ρ₀ in the
   inertia, the true ρ only in the weight.", where="Ch. 4 §4.9")`
5. `nb.recap("R15", "Why Dρ/Dt = 0", "Heat and salt do not diffuse on wave time scales, and the density depends on
   temperature and salinity but not on pressure (δρ/ρ = −αδT, +βδS), so a parcel keeps its density.", where="Ch. 1
   §1.10, Ch. 4 §4.2")`
6. `nb.recap("R16", "The resting state", "Before the waves, the fluid is hydrostatic: $0=-\\frac1{\\rho_0}\\frac{d\\bar p}
   {dz}-\\frac{\\bar\\rho g}{\\rho_0}$ (7.123).", where="Ch. 1 §1.7")`
7. `nb.recap("R17", "Background plus perturbation", "$p=\\bar p(z)+p'$, $\\rho=\\bar\\rho(z)+\\rho'$ (7.124) — ⚠️ here p′ is
   measured from the stratified background p̄(z), not from −ρgz as in (7.30).", where="Ch. 4 §4.9")`
8. `nb.recap("R18", "The buoyancy frequency", "$N^2\\equiv-\\frac g{\\rho_0}\\frac{d\\bar\\rho}{dz}$ (7.127) — (1.29) with
   no adiabatic gradient: a displaced parcel oscillates at N.", where="Ch. 1 §1.10")` + `nb.code`: `print(stratification.
   brunt_vaisala_sq(1025.0, -1.045e-2, 0.0, g=9.81))` → N² ≈ 1.0e-4 s⁻² (N = 0.01 rad/s: a period of 10.5 min, a typical
   thermocline).
9. `nb.recap("R19", "Potential density", "Compressibility is handled by using the potential density (the density brought
   adiabatically to a reference pressure, oceanographers' 'sigma-theta') in N, so the incompressible equations still
   apply.", where="Ch. 1 §1.10")`
10. `nb.recap("R20", "Internal waves are rotational", "Where density surfaces tilt against pressure surfaces the baroclinic
    torque ∇ρ × ∇p/ρ² makes vorticity (Ch. 5 §5.4), so internal waves have no velocity potential — C13's vortex sheet is
    the thin-layer limit.", where="Ch. 5 §5.4")`
11. `nb.core("C15", "Internal waves: $\\omega=N\\cos\\theta$ (7.139)", question="If the buoyancy frequency N is the only
    frequency the fluid knows, what sets the frequency of a wave?")`
12. `nb.md` — **The problem in plain words:** "Fill a tank with salt water that gets fresher toward the top and wiggle a
    rod: the disturbance does not spread in rings like on a pond but in slanted beams. In the ocean, tides flowing over
    ridges radiate such beams through the thermocline; in the atmosphere, wind over mountains makes lee waves. We want the
    rule that fixes their frequency."
13. `nb.md` — **The idea** (ASCII):
    ```
    parcel pushed straight up:           it bobs at N (Ch. 1)
    parcel pushed along a slope at θ:    only the vertical part of gravity's pull restores it
                                         ⇒ it bobs at N cos θ, whatever the wavelength
    so:  ω = N cos θ,   0 < ω < N,   θ = angle of K above the horizontal = angle of the motion from the vertical
    ```
14. `nb.note` — **N132 [B]** "Put the split into the density equation: $\frac{\partial}{\partial t}(\bar\rho+\rho')+u\frac
    {\partial}{\partial x}(\bar\rho+\rho')+v\frac{\partial}{\partial y}(\bar\rho+\rho')+w\frac{\partial}{\partial z}(\bar\rho+
    \rho')=0$ (7.125)." equation. Gloss (markdown, before D32): "**Linearisation about a base state** keeps products of one
    small quantity with a background quantity (w dρ̄/dz) and drops products of two small quantities (u ∂ρ′/∂x) — Ch. 1 P68."
15. `nb.derivation("D32", …)` — Part F D32 (10 steps), ref "7.131".
16. `nb.note` — **N133 [B], N134 [B], N135 [B]** "$\frac{\partial\rho'}{\partial t}+w\frac{d\bar\rho}{dz}=0$ (7.126): density
    changes at a point only because the vertical motion carries the background up or down; the perturbation momentum
    equations $\frac{\partial u}{\partial t}=-\frac1{\rho_0}\frac{\partial p'}{\partial x}$, $\frac{\partial v}{\partial t}=
    -\frac1{\rho_0}\frac{\partial p'}{\partial y}$, $\frac{\partial w}{\partial t}=-\frac1{\rho_0}\frac{\partial p'}{\partial z}
    -\frac{\rho'g}{\rho_0}$ (7.128, 7.129, 7.130) look like (7.120)–(7.122) with ρ′, p′; and $\frac{\partial\rho'}{\partial t}-
    \frac{N^2\rho_0}{g}w=0$ (7.131)."
17. `nb.primer("operator elimination for linear PDEs", "With several linear PDEs in several unknowns, apply derivatives
    (∂/∂t, ∂/∂z, the horizontal Laplacian ∇_H²) to whole equations and substitute one into another until a single unknown
    is left — like eliminating variables from linear equations. Allowed because the derivatives of smooth fields commute
    (Schwarz, Ch. 4 P121) and the coefficients do not depend on the variable being differentiated.", code="import sympy as
    sp   # symbolic check of one elimination move\nx, t = sp.symbols('x t'); u, p = sp.Function('u')(x, t),
    sp.Function('p')(x, t)\neq1 = sp.Eq(sp.diff(u, t), -sp.diff(p, x))   # ∂u/∂t = −∂p/∂x\nprint(sp.simplify(sp.diff(eq1.lhs,
    x) - sp.diff(sp.diff(u, x), t)))   # 0: ∂/∂x of ∂u/∂t is ∂/∂t of ∂u/∂x")` (**P177**).
18. `nb.derivation("D33", …)` — Part F D33 (12 steps, ★★★, with `check_src`), ref "7.134".
19. `nb.note` — **N136 [B], N137 [B], N138 [B]** "D33 passes through $\frac1{\rho_0}\nabla_H^2p'=\frac{\partial^2w}{\partial
    z\,\partial t}$ (7.132) ($\nabla_H^2\equiv\partial^2/\partial x^2+\partial^2/\partial y^2$), $\frac1{\rho_0}\frac{\partial^2p'}
    {\partial t\,\partial z}=-\frac{\partial^2w}{\partial t^2}-N^2w$ (7.133), and ends at the w-equation
    $\frac{\partial^2}{\partial t^2}\nabla^2w+N^2\nabla_H^2w=0$ (7.134) — which holds for N(z) too."
20. `nb.note` — **N139 [B], N140 [B], R20 link** "The medium has no preferred horizontal direction but a preferred vertical,
    so ω depends on all of K: $\omega=\omega(k,l,m)=\omega(\mathbf K)$ (7.135) — anisotropic. Try a plane wave
    $w=w_0e^{i(kx+ly+mz-\omega t)}=w_0e^{i(\mathbf K\cdot\mathbf x-\omega t)}$ (7.136)." Gloss (in the markdown): "**Angle of a
    vector from its components**: cos θ = |k|/K with `np.arccos`, degrees with `np.degrees`."
21. `nb.derivation("D34", …)` — Part F D34 (8 steps), ref "7.139".
22. `nb.note` — **N141 [B], N142 [B]** "$\omega^2=\frac{k^2+l^2}{k^2+l^2+m^2}N^2$ (7.137); in the x–z plane $\omega=\frac{kN}
    {\sqrt{k^2+m^2}}=\frac{kN}{K}$ (7.138) — written for k > 0; for a wave travelling toward −x use |k| (our code does)."
23. `nb.note` — **N145 [B]** "Limits: θ = 0 (K horizontal, m = 0) gives ω = N — vertical columns bobbing, Ch. 1's parcel;
    θ → π/2 gives ω → 0, where the wave solution says nothing and the full equations must be used (next note)." +
    `nb.code`: `t = np.linspace(0, 1000, 5); f = ch07.internal_wave_fields(0.0, 0.0, t, 0.01, 0.0, 0.01, 1e-3)` (m = 0) ·
    `z = stratification.parcel_displacement(t, f["zeta_particle"][0], 1e-4)` · `assert np.allclose(f["zeta_particle"], z,
    atol=1e-12)` — the ω = N limit is Ch. 1's parcel.
24. `nb.note` — **N146 [B]** "A possible steady solution of (4.10) and (7.128)–(7.131): $w=p'=\rho'=0$ with u, v any
    horizontally non-divergent field, $\frac{\partial u}{\partial x}+\frac{\partial v}{\partial y}=0$ (7.142), at each level
    separately — strong stratification allows flat, layered flows (pancake eddies, cloud sheets seen from aircraft, Ch. 13)."
    + `nb.code`: `print(ch07.layered_flow_check(lambda x, y: -np.sin(y), lambda x, y: np.sin(x), 0.3, 0.7))` → all residuals
    0. **N147 [C], N196 [C]** "In strong stratification a 2-D body blocks a whole horizontal layer of fluid ahead of it
    (Fig. 7.30: the fluid cannot go over or under); our schematic draws streamlines stopping between the two tangent planes.
    Orographic blocking in Ch. 13."
25. `nb.worked_example("a thermocline internal wave", "N = 0.01 rad/s. 1. Buoyancy period 2π/N = 628 s = 10.5 min —
    no internal wave can be faster. 2. A wave with K at θ = 45° to the horizontal: ω = N cos 45° = 0.00707 rad/s, T = 14.8
    min. 3. The same direction with a 10× longer wavelength: still 14.8 min — only the direction counts. 4. For T = 12 h
    (a tide): cos θ = (2π/43 200)/0.01 = 0.0145, θ = 89.2° — K almost vertical, the energy beams almost horizontal.")`
26. `nb.code` — `print(ch07.internal_wave_omega(1.0, 1.0, 0.01), ch07.internal_wave_omega(0.1, 0.1, 0.01),
    ch07.internal_wave_omega(-1.0, 1.0, 0.01))` · `print(np.degrees(ch07.beam_angle(0.00707, 0.01)))` · `s =
    ch07.boussinesq_linear_sympy(); print({k: v for k, v in s.items() if k.startswith("r")})`. *expect:* 0.0070711 three
    times (size and the sign of k do not matter); 45.0°; all residuals 0.
27. `nb.check_agree` — **from scratch:** the finite-difference residual of (7.134) for a plane wave w = cos(kx + mz − ωt) at a
    point (h = 1e-3, second differences in t, x, z): with ω from (7.138) it is ≈ 0 (≲ 1e-9 relative), with the wrong ω = N it
    is ≈ N²k² (≈ 1e-4 relative 1): `assert abs(res_right) < 1e-6*abs(res_wrong)`; and `np.allclose(res_right,
    ch07.w_equation_residual(w_fn, 0.3, -0.2, 1.0, 0.01), atol=1e-12)`.
28. `nb.figure` — two panels (9 × 3.6 in): (a) ω/N vs θ (cos curve, 0…90°) with dots at 0°, 45°, 89°; (b) the (k, m) plane
    with contours of ω/N = 0.2, 0.5, 0.8 (straight rays from the origin) and three K arrows of different length on the 45°
    ray. Title "Frequency is set by direction, not wavelength". *see:* "straight lines of constant ω through the origin."
    *read:* "moving along a ray (changing |K|) keeps ω; turning K changes it." *change:* "…a surface wave's ω(K) instead:
    contours would be circles (ω depends only on |K|) — the opposite geometry."
29. `nb.plotly` — **IF8:** slider over ω/N (0.05…0.95, 19 steps): the four beam directions and the c (orange) and c_g
    (purple) arrows of each beam (from `ch07.internal_wave_state`), the St Andrew's cross skeleton. Title "The angle is set
    by the frequency".
30. `nb.md` — "> ⚠️ **Common confusion:** θ in $\omega=N\cos\theta$ (7.139) is the angle of **K** above the horizontal.
    The same θ is the angle of the particle motion and of the energy beams from the **vertical**. The book's Fig. 7.33
    caption says '45° with the horizontal' — true for the beams only because 45° is symmetric. We always print both."
31. `nb.md` — **What would change if…** "…the fluid rotated at f (Ch. 13)? Then $\omega^2=N^2\cos^2\theta+f^2\sin^2\theta$:
    internal waves live between f and N — the inertia–gravity waves of the ocean and atmosphere."
32. `nb.core("C16", "Transverse waves with $\\mathbf c_g\\cdot\\mathbf c=0$ (7.146): beams and $\\mathbf F=\\mathbf c_gE$
    (7.159)", question="Why does an internal wave's energy leave at right angles to the way its crests move?")`
33. `nb.md` — **The problem in plain words:** "Oscillate a cylinder up and down in a stratified tank and four beams appear —
    the St Andrew's cross; inside each beam the crests (dark and light stripes) slide *across* the beam, not along it. In
    the ocean, energy from tides over a ridge goes down into the abyss in such beams while the phase lines move up. We
    want to understand this and prove where the energy goes."
34. `nb.md` — **The idea** (ASCII):
    ```
    incompressible plane wave   ⇒  K·u = 0          the water moves along the crests (transverse)
    ω = N|k|/K depends on direction only  ⇒  ∇_K ω ⊥ K   ⇒  c_g ⊥ c
    phase up  ⇔  energy down;   horizontal parts of c and c_g point the same way
    ```
35. `nb.note` — **N143 [B]** "Write the velocity as $u=u_0e^{i(kx+ly+mz-\omega t)}$ (7.140) and similarly v, w." equation.
36. `nb.derivation("D35", …)` — Part F D35 (5 steps), ref "7.141".
37. `nb.note` — **N144 [B]** "$\mathbf K\cdot\mathbf u=0$ (7.141): particle motion is perpendicular to K, i.e. along the crests
    — a transverse (shear) wave. Surface waves are not transverse: their fields decay in z instead of oscillating."
38. `nb.note` — **N148 [B]** "The group velocity is the gradient of ω in wavenumber space (N71), component by component:
    $\mathbf c=(\omega/K)\mathbf e_K$ and $\mathbf c_g=\mathbf e_x\frac{\partial\omega}{\partial k}+\mathbf e_y\frac{\partial
    \omega}{\partial l}+\mathbf e_z\frac{\partial\omega}{\partial m}$ (7.8, 7.143)." Gloss: "**gradient in wavenumber space**:
    treat (k, l, m) as coordinates and take partial derivatives (Ch. 1 P25)."
39. `nb.derivation("D36", …)` — Part F D36 (9 steps), ref "7.146".
40. `nb.note` — **N149 [B]** "$\mathbf c=\frac{\omega}{K^2}(k\mathbf e_x+m\mathbf e_z)$ and $\mathbf c_g=\frac{Nm}{K^3}(m\mathbf
    e_x-k\mathbf e_z)$ (7.144, 7.145) (for k > 0; for k < 0, as drawn in Fig. 7.29, the code uses ∇_K of N|k|/K)."
41. `nb.figure` — **Figs. 7.29 and 7.31 remake (N195, N197)** two panels: (a) the book's geometry with k < 0: K and c up-left
    at θ, c_g down-left at θ from the vertical, crest lines and the particle-motion double arrow along them; (b) c and c_g as
    the two legs of a right triangle with horizontal hypotenuse for ω/N = 0.3, 0.71, 0.9. *see:* "c and c_g at a right angle
    with opposite vertical parts." *read:* "the hypotenuse is horizontal: c_x and c_gx have the same sign, c_z = −c_gz."
    *change:* "…ω → N: θ → 0, K horizontal, c_g → 0 — vertical columns bob and nothing propagates."
42. `nb.note` — **N150 [B]** "**St Andrew's cross.** A source oscillating at ω < N can only radiate waves with cos θ = ω/N,
    so its energy leaves along four beams at θ from the vertical (Fig. 7.33: ω = 0.71N → 44.8°), phase lines moving across
    each beam." + `nb.figure` **Fig. 7.33 remake (N199)** `ch07.st_andrews_cross` field on a 400² grid (FAST 200²) for
    ω/N = 0.71 with the four c_g (purple) and c (orange) arrows and both angles labelled (to the vertical and to the
    horizontal) — "our superposition of four Gaussian beams, an illustration of the geometry, not the photograph". *see:*
    "an X of striped beams." *read:* "stripes are crests (constant ρ′); energy runs along the beam, phase across it."
    *change:* "…ω/N = 0.3: the beams tilt toward the horizontal (72.5° from the vertical)."
43. `nb.note` — **N151 [C]** "In the real ocean N depends on z (below 0.01 rad/s); the plane-wave results hold locally where
    N changes little over a vertical wavelength 2π/m, and rays (C10's method) bend as N changes — WKB theory in Ch. 13."
44. `nb.note` — **N152 [B]** "Multiply (7.128)–(7.130) by ρ₀u, ρ₀v, ρ₀w and add (the Ch. 4 energy move, (4.56)):
    $\frac{\partial}{\partial t}\big[\tfrac12\rho_0(u^2+v^2+w^2)\big]+g\rho'w+\nabla\cdot(p'\mathbf u)=0$ (7.147): kinetic
    energy changes by buoyancy work gρ′w and by the divergence of the energy flux p′u." + `nb.code`: `print(ch07.
    internal_energy_budget_residual(0.3, -0.2, 1.0, 1.0, 1.0, 0.01, 1e-3))` → ≈ 1e-12.
45. `nb.note` — **N153 [B], N154 [B], N155 [B]** "With (7.131) the buoyancy work is a time derivative:
    $\frac{\partial E_p}{\partial t}=g\rho'w=\frac{\partial}{\partial t}\big[\frac{g^2\rho'^2}{2\rho_0N^2}\big]$ (7.148); with
    w = ∂ζ/∂t, $\rho'=\frac{N^2\rho_0\zeta}{g}$ (7.149), so $E_p=\frac{g^2\rho'^2}{2\rho_0N^2}=\tfrac12N^2\rho_0\zeta^2$ (7.150)
    per unit **volume** — the available potential energy of Ch. 13. Number: N = 0.01 rad/s, ζ = 10 m → 5 J/m³."
46. `nb.note` — **N156 [C], N157 [B]** "Consistency with C13: for two deep fluids the potential energy per unit area was
    $\tfrac14(\rho_2-\rho_1)ga^2$ (7.151); writing a density jump as $N^2=-\frac g{\rho_0}\frac{d\bar\rho}{dz}=\frac g{\rho_0}
    (\rho_2-\rho_1)\delta(z)$ (7.152) and integrating (7.150) over z recovers it (gloss: a 1-D Dirac delta is the limit of a
    narrow step's derivative — Ch. 6 P149 in one dimension)." + `nb.code`: `for eps in (1.0, 0.3, 0.1, 0.03): print(eps,
    ch07.internal_pe_interface_limit(1.0, 1000.0, 1002.0, eps=eps))` → approaching ¼ × 2 × 9.81 × 1 = 4.905 J/m² as ε → 0
    (the builder prints the error at each ε and its observed order).
47. `nb.primer("mean of a product of real parts", "For two fields written as real parts of complex amplitudes, the
    average over a period of their product is $\\langle\\mathrm{Re}(Ae^{i\\theta})\\,\\mathrm{Re}(Be^{i\\theta})\\rangle=
    \\tfrac12\\mathrm{Re}(AB^*)$ (B* the complex conjugate, Ch. 3 P81). Two fields 90° out of phase (B = iA) give zero:
    they carry no mean product.", code="import numpy as np   # numbers\nA, B = 2.0, 3.0*np.exp(1j*0.5)   # amplitudes
    with a phase shift\nth = np.linspace(0, 2*np.pi, 100001)\nmean = np.mean(np.real(A*np.exp(1j*th))*np.real(B*np.exp(1j*th)))\n
    print(mean, 0.5*np.real(A*np.conj(B)))   # both 2.633")` (**P178**).
48. `nb.note` — **N158 [B]** "Assume $[u,w,p',\rho']=[\hat u,\hat w,\hat p,\hat\rho]e^{i(kx+mz-\omega t)}$ with real ŵ." (D37
    step 1)
49. `nb.derivation("D37", …)` — Part F D37 (12 steps), ref "7.159".
50. `nb.note` — **N159 [B]** "The polarization relations $p'=-\frac{\omega m\rho_0}{k^2}\hat we^{i(kx+mz-\omega t)}$,
    $\rho'=\frac{iN^2\rho_0}{\omega g}\hat we^{i(kx+mz-\omega t)}$, $u=-\frac mk\hat we^{i(kx+mz-\omega t)}$ (7.153): ρ′ is 90°
    out of phase with w (the factor i)." **N160–N163 [B]** "$E_k=\tfrac12\rho_0\overline{(u^2+w^2)}=\tfrac14\rho_0\big(\frac
    {m^2}{k^2}+1\big)\hat w^2$ (7.154), $E_p=\frac{g^2\overline{\rho'^2}}{2\rho_0N^2}=\frac{N^2\rho_0}{4\omega^2}\hat w^2$ (7.155),
    equal (7.156) because $\omega^2=k^2N^2/(k^2+m^2)$, so $E=E_k+E_p=\tfrac12\rho_0\big(\frac{m^2}{k^2}+1\big)\hat w^2$ (7.157)
    (the same ½Re(AB*) move twice)." **N164–N166 [B]** "$\mathbf F=\overline{p'\mathbf u}=\mathbf e_x\overline{p'u}+\mathbf e_z
    \overline{p'w}=\frac{\rho_0\omega m\hat w^2}{2k^2}\big(\mathbf e_x\frac mk-\mathbf e_z\big)$ (7.158), and
    $\mathbf c_gE=\frac{Nm}{K^3}(m\mathbf e_x-k\mathbf e_z)\big[\frac{\rho_0}2\big(\frac{m^2}{k^2}+1\big)\hat w^2\big]$ reduces to
    it with (7.138): $\mathbf F=\mathbf c_gE$ (7.159) — the same law as (7.71)."
51. `nb.md` — "> ⚠️ **Common confusion (units):** in (7.159) F is per unit **area** and E per unit **volume**; in (7.71) F was
    per metre of crest (the whole depth) and E per unit horizontal area. Same law, different bookkeeping."
52. `nb.worked_example("c and c_g at 45°", "N = 1 rad/s, K = 1 rad/m at θ = 45° (k = m = 0.707 rad/m). 1. ω = N cos 45° =
    0.707 rad/s. 2. $\\mathbf c=\\frac{\\omega}{K^2}(k,m)=0.707\\times(0.707,0.707)=(0.5,0.5)$ m/s: up and to the right. 3.
    $\\mathbf c_g=\\frac{Nm}{K^3}(m,-k)=0.707\\times(0.707,-0.707)=(0.5,-0.5)$ m/s: down and to the right. 4. c·c_g =
    0.25 − 0.25 = 0: a right angle; |c| = |c_g| = 0.707 m/s at this angle only. 5. ω/N = 0.71 → θ = arccos 0.71 = 44.8°.")`
53. `nb.code` — `v = ch07.internal_wave_velocities(0.7071, 0.7071, 1.0); print(v)` · `print(ch07.internal_wave_velocities(-0.7071,
    0.7071, 1.0))` (Fig. 7.29's k < 0) · `e = ch07.internal_wave_energy(0.7071, 0.7071, 1.0, 0.01); print(e)` · `assert
    np.allclose(e["F"], e["cgE"]) and np.isclose(e["Ek"], e["Ep"])` · `print(ch07.internal_wave_state(0.71, 1.0))`. *expect:*
    c = (0.5, 0.5), c_g = (0.5, −0.5), dot 0; for k < 0: c = (−0.5, 0.5) (up-left), c_g = (−0.5, −0.5) (down-left, as drawn in Fig. 7.29; the printed
    (7.145) $\frac{Nm}{K^3}(m\mathbf e_x-k\mathbf e_z)$ would give (0.5, 0.5), pointing up-right — wrong); E_k = E_p = 0.05 J/m³ (ρ₀ = 1000, ŵ = 0.01 m/s), F = c_gE; beam 44.77° from the vertical.
54. `nb.check_agree` — **from scratch:** central-difference gradient of ω(k, m) = N|k|/√(k² + m²) at (−0.7071, 0.7071), then
    the dot product with c: `assert np.allclose(cg_mine, ch07.internal_wave_velocities(-0.7071, 0.7071, 1.0)["cg"],
    rtol=1e-7)` and `assert abs(np.dot(cg_mine, c)) < 1e-12`.
55. `nb.animation` — **A6 (N198)** a 2-D internal-wave packet from `ch07.linear_evolve_2d` (256², FAST 128²; K₀ at 45°
    up-right, N = 1 rad/s, Gaussian envelope), `frames=60 if not FAST else 30`: phase lines move up-right while the packet
    moves down-right, sliding along its own crests; purple and orange arrows of c_g and c. *see / read / change*: "the
    packet never moves the way its stripes move; at ω → N it would stand still."
56. `nb.explainer("internal_wave_beams", heading="Why does energy leave at right angles to the crests?", why="The ω/N slider
    tilts the beams while the crests visibly slide across them; drag K in the wavenumber plane: ω stays fixed when only |K|
    changes, and c_g stays perpendicular to c.", tries=["Preset *45° cross*: read both angles in the status.", "Drag the K
    arrow outward along its ray: ω does not change.", "Choose *near N*: the beams stand almost vertical and the energy
    hardly moves.", "Set the book's Fig. 7.29 geometry (k < 0) and check that c_g points down-left."])`
57. `nb.md` — **What would change if…** "…N varied with depth, as in the real ocean (N largest in the thermocline)? The beams
    curve, like the rays of C10, and can reflect where ω = N — internal waves get trapped in the thermocline (Ch. 13)."
58. `nb.pointer("**Exercises 7.1–7.20** are not reproduced (S01). Ideas from them used in this notebook: 7.3 → N02 (any
    initial shape); 7.4 → N30 (ψ); 7.5–7.6 → N63 (basin modes); 7.9–7.10 → N69, N72 (capillary c_g, minimum c_g); 7.11 → N72
    (viscous decay); 7.12 → D20 (Gaussian packet); 7.13–7.14 → N91 (drift at any depth); 7.15 → N97 (soliton check); 7.18 →
    N105 (interfacial energies); 7.19 → D30 (two-layer dispersion); 7.20 → N131 (rigid lids).")`
59. `nb.summary(clicked=[16 bullets, one per CORE: C01 "A wave is a moving phase: crests move at c = ω/k while the water
    only oscillates." · C02 "For gentle waves the surface conditions can be moved to the flat level z = 0; the error is of
    order ka." · C03 "Laplace plus the bottom, kinematic and dynamic conditions leave one frequency per wavenumber:
    ω² = gk tanh kH." · C04 "Long waves outrun short ones; depth matters only for waves longer than about twice the depth,
    and all long waves travel at √(gH)." · C05 "Under a wave the water goes round closed ellipses — circles in deep water,
    flat ellipses in shallow water." · C06 "A wave holds ½ρga² per square metre, half kinetic and half potential, and
    carries it at a speed that is not the crest speed." · C07 "Surface tension adds a second restoring force; together with
    gravity it makes a slowest wave, 23 cm/s at 1.7 cm." · C08 "Two opposite waves make a standing wave; walls keep only the
    patterns whose nodes fit, so a lake rings at its own periods." · C09 "Groups and energy travel at c_g = dω/dk — half the
    crest speed in deep water — and F = E c_g." · C10 "Along a ray moving at c_g the frequency never changes; the wavelength
    shrinks in shallow water, so crests turn parallel to the shore." · C11 "A hydraulic jump's height comes from momentum
    alone; energy is lost, which is why only fast shallow flow can jump." · C12 "At finite amplitude the orbits do not close:
    the water drifts forward at a²ωk e^{2kz₀} while the mean at a fixed point stays zero." · C13 "A small density step gives
    slow, tall waves, ε√(gk), with the two fluids sliding past each other." · C14 "A layer over deep water rings in two
    modes: barotropic (fast, surface) and baroclinic (slow, √(g′H), almost invisible at the surface)." · C15 "In a stratified
    fluid the frequency is set by direction only: ω = N cos θ, never above N." · C16 "Internal waves are transverse and
    their energy moves at right angles to the crests: phase up, energy down, F = c_g E."], feeds_forward=["Ch. 11: the
    interfacial vortex sheet with a shear flow is the Kelvin–Helmholtz instability; the σ-term of (7.56) stabilises short
    waves.", "Ch. 13: √(gH) and hydrostatic shallow water carry tides, tsunamis and Kelvin waves; rotation adds f to ω²
    (Poincaré and inertia–gravity waves); √(g′H)/f is the baroclinic Rossby radius; rays and WKB follow waves through N(z);
    equipartition fails in geostrophic adjustment.", "Ch. 15: characteristics and shocks are the compressible cousins of
    N81 and the hydraulic jump."], left_out=["Exercise text (S01).", "Nonlinear internal waves and wave breaking (research
    literature).", "Viscous damping in detail (Ch. 8 tools; named in N72)."])`

### A.9 Placement check (every curation id has exactly one home)
| IDs | Block | Section |
|---|---|---|
| C01, N01–N11, N167, D01 | C01 | §7.1 |
| R01–R06 | before C02 | §7.2 |
| C02, N12–N17, N168, D02–D04 | C02 | §7.2 |
| C03, N18–N25, D05, D06 | C03 | §7.2 |
| R07 | before C04 | §7.2 |
| C04, N26, N37, N38, N39, N42, N43, N46, N173, D07–D09, E1 | C04 | §7.2 |
| R08 | before C05 | §7.2 |
| C05, N27–N30, N40, N41, N44, N45, N169–N171, D10, D11, A1 | C05 | §7.2 |
| C06, N31–N36, N172, D12–D14 | C06 | §7.2 (+ pointer to C10 for N47, N48, N174, N175) |
| R09 | before C07 | §7.3 |
| C07, N49–N57, N176, D15, D16, E3 | C07 | §7.3 |
| C08, N58–N63, N177, N178, D17, D18, A3, E4 | C08 | §7.4 |
| C09, N64–N73, N179–N182, D19–D21, A2, E5 | C09 | §7.5 |
| C10, N47, N48, N74–N80, N174, N175, N183, N184, D22–D24, E6 | C10 | §7.5 |
| R10 | before C11 | §7.6 |
| C11, N81–N85, N93–N97, N185, N186, N189, D25, D26, A5, E7 | C11 | §7.6 |
| C12, N86–N92, N187, N188, D27, A4, E2 | C12 | §7.6 |
| R11 | before C13 | §7.7 |
| C13, N98–N108, N190–N192, D28 | C13 | §7.7 |
| C14, N109–N131, N193, N194, D29–D31, E8 | C14 | §7.7 |
| R12–R20 | before C15 | §7.8 |
| C15, N132–N142, N145–N147, N196, D32–D34 | C15 | §7.8 |
| C16, N143, N144, N148–N166, N195, N197–N199, D35–D37, A6, E9 | C16 | §7.8 |
| S01 | pointer after C16 | §7.8 (end) |
Totals: 16 CORE blocks, 199 NOTE ids, 20 RECAP calls, 1 SKIP pointer, 37 derivations, 9 explainers, 6 animations, 8
plotly figures, 1 live cell, 14 primers (P165–P178).

---

## Part B — explainer storyboards

Common to all ten: created with `tools/new_viz.py`; `<meta name="viz:chapter" content="ch07">`; tabs Walkthrough ·
Explore · Explain · Derivation · Equations · Code · Check; every displayed number is computed by a JS function that
mirrors a `ch07` callable and is proved by `selftest()` parity rows (`py:` expressions use only `ch07.…`, `np.pi`,
`np.inf`, numbers, strings, lists and dicts — Part C convention 3; results indexed down to one number). Explain is
"Explanation & interpretation" in numbered sections built with `Viz.work.step / line / box / table / hint / interpret`,
modelled on `forced_damped_vibrations.html` and `amplitude_phase_second_order_II_3.html`: **0** what the views show and
what each colour means · **1…n** every displayed quantity from the controls ("formula = substituted = result — why",
results boxed) · a section per view hidden on phones, or a hint · the values at the current time (live) · **Reading the
current setting** (regime-dependent). Derivation steps are copied from Part F (same `did` titles, same step count; phones
shorten *why* to its first sentence; plain-text *why* and *watch* never contain raw TeX — write e^(kz), not
`e^{kz}`). Every tour, Explain, Derivation, notes, status, equation and quiz text that names a book equation **writes it
out** next to its number (convention 9; `tools/eq_refs.py` → 0). Drafts below write equations in Unicode for readability;
builders set each in TeX (`$…$`, backslashes doubled in JS strings). Colours as in convention 8 (surface blue, particles
teal, phase/crests orange, group/energy purple, pressure amber, surface tension rose, lower layer navy). Walkthrough texts
≤ 45 words. A view hidden on portrait phones never carries a step's key number (repeat it in a visible title or readout).
**Time:** explainers whose periods span decades run the transport in **periods** (`tau` = t/T) and display the real
time t = τT with `Viz.fmtTime`. **Hyperbolic functions** use the overflow-safe ratios of Part C (a local `ratioCS(k, z,
H)` = e^{kz}(1 + e^{−2k(z+H)})/(1 − e^{−2kH}), `H = Infinity` handled). Each explainer fits 360×640 … 1920×1080 and the
1000×700 notebook frame with no scrolling (fit plan per explainer). Parallel builders use private scratch subfolders.

### E1 · dispersion_relation
- **Title:** "Why do long waves outrun short ones?" · **Summary:** "Drag the wavelength and the depth: the crest marker in
  the tank, the dot on the c(λ) curve and the pressure profile move together between the shallow line √(gH) and the deep
  curve √(g/k)." · **CORE:** C03, C04 (also C01's c = ω/k (7.4), N26 (7.31), N37, N38 (7.45), N39, N42 (7.48), N43
  (7.49), N46 (7.52)) · **Reference:** `angular_frequency_explorer_1.html` (system animation + curves on one clock,
  presets, "Right now" notes with a highlighted table of real cases).
- **meta:** `viz:order 1` · `viz:sections 7.1 7.2` · `viz:equations 7.4 7.28 7.29 7.31 7.45 7.48 7.49 7.52` · `viz:fluidpy
  ch07.dispersion_state ch07.omega_gravity ch07.phase_speed ch07.group_velocity ch07.pressure_response ch07.depth_regime
  ch07.wavenumber_from_omega` · `viz:derivations D01 D05 D06 D07 D08 D09`.
- **Physics (JS ↔ Python):** `omega(k, H)` = √(gk tanh kH) (deep branch √(gk) when H = Infinity) ↔ `ch07.omega_gravity(k,
  H)`; `cph(k, H)` ↔ `ch07.phase_speed(k, H)`; `cg(k, H)` = (c/2)(1 + 2kH/sinh 2kH) ↔ `ch07.group_velocity(k, H)`;
  `presp(k, z, H)` = cosh k(z + H)/cosh kH via exponentials ↔ `ch07.pressure_response(k, z, H)`; `regime(k, H)` returns
  {regime, deep_error = 1 − √(tanh kH), shallow_error = 1 − √(tanh kH/kH)} ↔ `ch07.depth_regime(k, H)`; `kFromT(T, H)`
  (`Viz.num.brentq` between the deep root ω²/g·0.5 and max(ω²/g, ω/√(gH))·2) ↔ `ch07.wavenumber_from_omega(2π/T, H)`;
  `state(lam, H)` bundles them ↔ `ch07.dispersion_state(lam, H)`. g = 9.81, ρ = 1000.
- **Views** (rows [1.15, 1]):
  1. `tank` "The wave in the tank" (row 0, flex 1.4) — x ∈ [0, 2λ] (axis in m, or km above 5 km), z from −min(H, 1.2λ) to
     +0.25 of the drawn depth (a vertical break mark ≈ when H > 1.2λ, labelled "bottom at −H = … m, too deep to draw");
     amplitude drawn at a fixed 8 % of the view height (true a only in the inspector). Surface η = a cos(kx − ωt) (blue),
     bottom (hatched), pressure amplitude shading |p′|/ρga = cosh ratio as an amber alpha ramp under the surface, an orange
     crest marker on one crest with an orange arrow c, five teal orbit ghosts at the surface…−H/2 (circle/ellipse size from
     (7.36)). Animates with the transport. Pointer: click → probe (x, z) (black ✚) for the inspector.
  2. `curve` "Phase speed c(λ)" (row 0, flex 1) — log–log, λ from 0.5 m to 10⁷ m, c from 0.5 to 400 m/s; the current curve
     (black), the deep line √(gλ/2π) (muted dashed), the shallow line √(gH) (muted dotted), shaded bands where each limit is
     within 2 % / 3 % (light blue for deep, light sand for shallow), the current dot (orange) with a vertical guide; chips in
     the title switch to ω(k) (k from 10⁻⁶ to 10 rad/m, the asymptotes √(gk) and k√(gH)). Pointer: click on the curve →
     set λ.
  3. `profile` "Pressure under the crest" (row 1, `hidePortrait: true`) — p′/ρga from 0 to 1 (x) against z/H from −1 to 0
     (y): the cosh ratio (amber), the hydrostatic line 1 (muted dotted), e^{kz} (muted dashed); a dot at the bottom value.
  Portrait: `tank` + `curve` stacked; the bottom-pressure fraction is repeated in the `tank` title ("bottom sees 53 %").
- **Controls (≤ 5 visible):** `lam` "Wavelength $\lambda$" log 0.5 m … 10⁷ m, default 156.1 m, help "crest to crest" ·
  `H` "Depth $H$" log 0.5 m … 5000 m, default 4000 m, help "still-water depth" · `Tset` "Keep the period $T$ fixed"
  toggle (optional, default off; when on, changing H re-solves λ from T with `kFromT`, the way swell shoals) · `a`
  "Amplitude $a$" 0.1 … 5 m, default 1 m (optional; affects only the pressure numbers) · `curve` chips "c(λ) · ω(k)".
- **Transport:** `tau` 0 → 4 periods, rate 0.5 period/s, `end: 'loop'`; readout t = τT.
- **Presets:** "wind swell" {lam: 156.1, H: 4000} · "swell on a 10 m beach" {Tset: true, H: 10} (λ → 92.4 m) · "tsunami"
  {lam: 200000, H: 4000} · "tide" {lam: 8.86e6, H: 4000} · "lab tank" {lam: 1, H: 0.5} · "kH = 2 (deep threshold)"
  {lam: 31.42, H: 10} · "H = 0.07λ (shallow threshold)" {lam: 142.9, H: 10}.
- **Status:** "🌊 deep: kH = 4.03 > 2 — the bottom is invisible, c = √(g/k) to 0.03 %" · "↔ intermediate: kH = 1.26 —
  both g/k and H matter" · "🏖️ shallow: H/λ = 0.02 < 0.07 — c = √(gH) = 198 m/s for every λ" (thresholds from
  `regime`).
- **Readouts:** "Frequency $\omega$" (rad/s) · "Period $T$" (s, `Viz.fmtTime`) · "Phase speed $c$" (m/s) · "Group speed
  $c_g$" (m/s) · "Bottom $p'$ / surface" (%).
- **Depth features:** Explain + Code + Derivation · **linked views** (3) · **transport** · **presets** · **status** ·
  **inspector** (click the tank: "p′ = ρga cosh k(z + H)/cosh kH cos(kx − ωt) = 1000 × 9.81 × 1 × cosh(0.126 × 5)/cosh
  1.257 × cos(…) = 9810 × 0.634 × 0.72 = **4478** Pa") · **notes** ("Right now": a table of real waves — lab tank,
  swell, swell on a shelf, tsunami, tide — with λ, H, kH, c; the nearest row highlighted).
- **Explain** ("Explanation & interpretation"):
  0. *What the views show* — "**The tank**: the blue line is the surface, the orange dot rides one crest at the phase
     speed, the amber shading is how strongly the wave's pressure is felt at each depth (dark = full), the teal loops are
     the paths of the water (C05). **Phase speed c(λ)**: the black curve is (7.29) for your depth; dashed = deep water
     (7.45), dotted = shallow water (7.49); the orange dot is your wave. **Pressure under the crest** (hidden on phones)
     shows p′ against depth."
  1. *Wavenumber and depth ratio* — "k = 2π/λ = 2π/**live lam** = **live k** rad/m; kH = **live kH** (H/λ = **live
     HoL**)." (why: every result below depends on k and H only through kH and g/k.)
  2. *Frequency from the dispersion relation* — "(7.28) ω = √(gk tanh kH) = √(9.81 × **k** × tanh **kH**) = √(9.81 × **k**
     × **tanh**) = **ω** rad/s" boxed; "T = 2π/ω = **T** s".
  3. *Phase speed* — "(7.4) c = ω/k = **ω**/**k** = **c** m/s" boxed; "the same as (7.29) √((g/k) tanh kH)"; "deep limit
     √(g/k) = **cdeep** m/s (off by **deep_error** %), shallow limit √(gH) = **cshal** m/s (off by **shallow_error** %)".
  4. *Group speed* — "(7.69) c_g = (c/2)(1 + 2kH/sinh 2kH) = **cg** m/s = **ratio** × c (C09: energy moves at c_g)".
  5. *Pressure below* — "(7.31) at the bottom p′/ρga = 1/cosh kH = 1/cosh **kH** = **pbot** — a bottom gauge sees
     **pbot%** of the surface pressure ρga = **pa** Pa."
  6. *How long to cross 1000 km* — "crests: 10⁶/c = **tc**; energy: 10⁶/c_g = **tcg** (the energy arrives later when
     c_g < c)."
  7. *Right now* — "t = **live t** (τ = **live tau** periods); the tracked crest is at x = **live xc** m."
  8. *Reading the current setting* — deep: "Deep water: the wave is shorter than about twice the depth (kH > 2), so the
     bottom never enters: c = √(g/k), longer waves are faster (dispersion), and the pressure has died to e^(kz) well above
     the bottom. A bottom gauge would not notice this wave." · intermediate: "Between the limits: both the wavelength and
     the depth set the speed; this is where swell is when it 'feels the bottom' near the coast and starts to slow and
     steepen." · shallow: "Shallow water: the wave is at least 14 times longer than the depth, so every wavelength travels
     at √(gH) — no dispersion, the pressure is hydrostatic (ρgη at every depth), and a bottom gauge sees the full signal.
     Tsunamis and tides live here."
- **Derivation tab:** **D01** (5 steps) `view: 'tank'`, goal page `set` preset "wind swell"; step 4 **live** "Δx = (ω/k)Δt
  = **c** × 1 s = **c** m" with `watch` "the orange dot moves one c per second"; interpret `s => "With your numbers a crest
  covers ${c} m every second, ${lam} m every period."`. **D05** (11 steps) `view: 'profile'` (phones: `tank`), goal `set
  {lam: 50, H: 10}`; step 8 (regroup into cosh k(z + H)) `set` highlights the amber profile with `watch` "the profile is
  the cosh ratio"; step 11 **live** "φ amplitude aω/(k sinh kH) = …". **D06** (7 steps) `view: 'curve'`, step 4 (cancel the
  cosine) highlights the orange dot, step 6 **live** "ω² = 9.81 × k × tanh kH = …"; **D07** (5 steps) `view: 'curve'`, step 5
  `watch` "the black curve rises to the right: longer is faster". **D08** (7 steps) `set {lam: 20, H: 10}` (kH = 3.1), step
  3 **live** "tanh 3.14 = 0.996: c within 0.2 % of √(g/k)"; **D09** (7 steps) `set {lam: 300, H: 10}` (kH = 0.21), step 3
  **live** "c/√(gH) = 1 − (kH)²/6 = …".
- **Code:**
  ```python
  lam, H = {{lam}}, {{H}}                    # wavelength and depth [m]
  k = 2*np.pi/lam                            # wavenumber = {{k}} rad/m
  omega = ch07.omega_gravity(k, H)           # (7.28) ω = √(gk tanh kH) = {{w}} rad/s
  T = 2*np.pi/omega                          # period = {{T}} s
  c = ch07.phase_speed(k, H)                 # (7.29) c = ω/k = {{c}} m/s
  cg = ch07.group_velocity(k, H)             # (7.69) c_g = {{cg}} m/s
  print(ch07.depth_regime(k, H))             # {{regime}}, deep err {{de}}, shallow err {{se}}
  print(ch07.pressure_response(k, -H, H))    # (7.31) bottom/surface = {{pb}}
  ```
- **Walkthrough (7 steps):** 1. "Two waves, two speeds" — "A tsunami crosses the Pacific in a day; swell needs a week.
  Same water, same gravity. What sets a wave's speed?" `set` wind swell, `play: true` · 2. "What moves" — "The orange dot
  rides one crest: it moves at c = ω/k (7.4), 15.6 m/s here. The water only goes round the teal loops." `readouts: ['c']`,
  `derive: {id: 'D01', step: 4}` · 3. "The rule" — "The surface and the bottom allow only ω = √(gk tanh kH) (7.28). Its
  slope from the origin is c." `eq: 'disp'`, `code: {lines: [2, 3]}`, `derive: {id: 'D06', step: 6}` · 4. "Deep water" —
  "Here kH = 4: tanh kH ≈ 1, so c = √(g/k) (7.45) — longer waves go faster and the depth is invisible." `set` wind swell,
  highlight `view:curve` · 5. "Shallow water" — "Drag λ to 200 km: now c = √(gH) (7.49) = 198 m/s for every long wave, and
  the pressure reaches the bottom unchanged." `set` tsunami, `controls: ['lam']` · 6. "Feel the bottom" — "Keep T = 10 s and
  bring the swell onto a 10 m beach: λ shrinks to 92 m and c to 9.2 m/s." `set` swell on a 10 m beach, `notes: true` · 7.
  "Your turn" — "Predict: at which λ does a 20 m-deep bay stop being 'deep' for waves? Drag λ, then check the status."
  `controls: ['lam', 'H']`.
- **Equations:** `c` "Phase speed" ref 'Eq. (7.4)' `c=\omega/k=\lambda\nu` · `disp` "Dispersion relation" ref 'Eq. (7.28)'
  `\omega=\sqrt{gk\tanh kH}` live "= √(9.81 × k × tanh kH) = … rad/s" · `cph` "Phase speed" ref 'Eq. (7.29)'
  `c=\sqrt{\frac gk\tanh kH}` · `deep` "Deep water" ref 'Eq. (7.45)' `c=\sqrt{g/k}=\sqrt{g\lambda/2\pi}` · `shal` "Shallow
  water" ref 'Eq. (7.49)' `c=\sqrt{gH}` · `press` "Pressure" ref 'Eq. (7.31)' `p'=\rho ga\frac{\cosh k(z+H)}{\cosh
  kH}\cos(kx-\omega t)` · symbols λ, k, ω, T, c, H, g, a, p′ with units.
- **Check yourself:** (1) "Double λ in deep water. By what factor does c change?" — "√2: c = √(gλ/2π) (7.45)." `set
  {lam: 156, H: 4000}` · (2) "Double H for a 200 km tsunami. By what factor does c change?" — "√2 again, but now because
  c = √(gH) (7.49): long waves do not care about λ." · (3) "Why does the bottom gauge see 100 % of the tsunami but almost
  nothing of the swell?" — "Pressure decays like cosh k(z + H)/cosh kH (7.31): for kH ≪ 1 the ratio is ≈ 1 (hydrostatic),
  for kH = 4 it is 1/cosh 4 ≈ 4 %." · (4) "With T fixed at 10 s, what happens to λ as the swell moves into shallower
  water?" — "It shrinks (156 → 137 → 92 m at 4000, 30, 10 m): ω is fixed, so k must grow as tanh kH falls." `set
  {Tset: true}`.
- **Selftest parity rows:** `{name: 'omega k=0.1 H=10', js: omega(0.1, 10), py: 'ch07.omega_gravity(0.1, 10.0)', rtol:
  1e-12}` · `{name: 'c swell', js: state(156.1, 4000).c, py: 'ch07.dispersion_state(156.1, 4000.0)["c"]', rtol: 1e-12}` ·
  `{name: 'cg lam50 H10', js: cg(2*Math.PI/50, 10), py: 'ch07.group_velocity(2*np.pi/50, 10.0)', rtol: 1e-12}` · `{name:
  'bottom pressure', js: presp(2*Math.PI/50, -10, 10), py: 'ch07.pressure_response(2*np.pi/50, -10.0, 10.0)', rtol:
  1e-12}` · `{name: 'k from T=10 H=30', js: kFromT(10, 30), py: 'ch07.wavenumber_from_omega(2*np.pi/10, 30.0)', rtol:
  1e-9}` · `{name: 'deep error kH=2', js: regime(0.2, 10).deep_error, py: 'ch07.depth_regime(0.2, 10.0)["deep_error"]',
  rtol: 1e-12}` · `{name: 'deep H=inf', js: omega(1, Infinity), py: 'ch07.omega_gravity(1.0, np.inf)', rtol: 1e-12}` ·
  invariant `{name: 'c = omega/k', js: cph(0.3, 7)*0.3 - omega(0.3, 7), expect: 0, atol: 1e-12}`.
- **Fit plan:** 360×640: status one line ("🌊 deep · kH 4.0"), `tank` (55 %) + `curve` (45 %), `profile` hidden (its
  number in the tank title); presets as a chip row that wraps; walkthrough card paged. 844×390: `tank` | `curve` side by
  side. 1000×700 and desktop: rows [1.15, 1], `profile` beside `curve` in row 1 with flex 0.6.

### E2 · particle_orbits
- **Title:** "What does the water under a wave actually do?" · **Summary:** "Tracer particles loop under a moving
  surface: circles in deep water, flat ellipses in shallow water; switch to exact path lines and the loops open — the
  Stokes drift — while a fixed current meter still reads zero." · **CORE:** C05, C12 (also N27–N30 ((7.33)–(7.37)), N40
  (7.46), N41 (7.47), N44 (7.50), N45 (7.51), N86 (7.82), N89–N92 ((7.84)–(7.86))) · **Reference:**
  `forced_damped_vibrations.html` (a phase-plane-like orbit view + a time graph with toggled components and the
  explanation panel).
- **meta:** `viz:order 2` · `viz:sections 7.2 7.6` · `viz:equations 7.33 7.34 7.35 7.36 7.37 7.46 7.50 7.84 7.85 7.86` ·
  `viz:fluidpy ch07.orbit_semi_axes ch07.orbit_linear ch07.particle_path ch07.stokes_drift ch07.orbit_state
  ch07.dyed_line` · `viz:derivations D10 D11 D27`.
- **Physics:** fixed wavelength λ = 20 m (k = 0.3142 rad/m), g = 9.81; H = kH/k; a = ka/k. `axes(z0)` → {A, B, focal}
  ↔ `ch07.orbit_semi_axes(z0, a, k, H)`; `orbitLin(x0, z0, t)` ↔ `ch07.orbit_linear`; `rhsExact(t, [x, z])` = (7.33)
  integrated by `Viz.num.rk4Step` with 200 steps per period ↔ `ch07.particle_path(…, model="exact")` (DOP853; parity at
  rtol 1e-6 after one period); `drift(z0)` = a²ωk cosh 2k(z₀ + H)/(2 sinh² kH) ((7.86); deep limit (7.85)) ↔
  `ch07.stokes_drift(z0, a, k, H)`; `psi(x, z, t)` (7.37); `state(z0)` ↔ `ch07.orbit_state(z0, a, k, H)`.
- **Modes (`model` chips):** "linear (7.34)" (closed orbits from the formula) · "exact (7.33)" (RK4 path lines; loops open).
- **Views** (rows [1.2, 1]):
  1. `column` "Water under the wave" (row 0, flex 1.5, equal aspect) — x ∈ [−λ/2, λ/2], z ∈ [−min(H, λ), a + 0.5 m];
     surface (blue), bottom (hatched, when within λ), 5 columns × 6 depths of tracer particles (teal dots) with their orbit
     ghosts (faint teal), a dyed vertical line through the middle column (rose) redrawn each frame, optional streamlines of
     (7.37) (thin grey, toggle); the selected particle ringed black with its (ξ, ζ) arrows. Pointer: click → select the
     nearest particle (its z₀ becomes `z0sel`).
  2. `axes` "Orbit size vs depth" (row 1, flex 1) — A(z₀) (teal solid), B(z₀) (teal dashed) against z₀ from −H (or −λ) to 0
     (y axis), and in exact mode ū_L(z₀) (purple, top axis in cm/s); a black dot at `z0sel`.
  3. `track` "x of the selected particle" (row 1, flex 1.2, `hidePortrait: true`) — x_p(t) − x₀ over the run: linear
     (teal, closed oscillation) and exact (purple, oscillation + drift), the drift line ū_L t (purple dashed); the
     "so far" part bold, the rest faint.
  Portrait: `column` + `axes`; the drift per period is repeated in the `column` title in exact mode.
- **Controls:** `kH` "Depth ratio $kH$" log 0.2 … 6, default 3 · `ka` "Steepness $ka$" 0.01 … 0.3, default 0.1 · `model`
  chips · `stream` "Show streamlines" toggle (optional) · `z0sel` "Selected depth $z_0$" −H … 0 (optional; also set by
  clicking).
- **Transport:** `tau` 0 → 10 periods, rate 0.5 period/s, `end: 'hold'` with an end-of-run card (`Viz.card`): "In 10
  periods (T = **T** s) the selected particle drifted **d** m — ū_L × 10T = **pred** m from (7.86); a fixed meter at its
  depth measured 0."
- **Presets:** "deep water" {kH: 3, ka: 0.1, model: 'linear'} · "intermediate" {kH: 1} · "shallow" {kH: 0.3} · "steep deep
  wave (exact)" {kH: 4, ka: 0.25, model: 'exact'} · "bottom particle" {kH: 0.5, z0sel: −H}.
- **Status:** "⭕ deep: circles, radius halves every 2.2 m" · "⬭ intermediate: B/A = 0.46 at z₀ = −H/2" · "▭ shallow:
  flat ellipses, width a/kH = 3.3a at every depth" · exact mode adds "➡ drift ū_L = 5.6 cm/s at the surface (ka = 0.1)".
- **Readouts:** "Semi-axis $A$" (m) · "Semi-axis $B$" (m) · "Drift $\bar u_L$" (cm/s) · "Period $T$" (s).
- **Depth features:** Explain + Code + Derivation · **linked views** (3) · **transport** + end-of-run card · **modes**
  (linear / exact) · **presets** · **status** · **inspector** (click a particle: "A = a cosh k(z₀ + H)/sinh kH = 0.318 ×
  cosh(0.314 × 7.55)/sinh 3 = **0.172** m; B = … ; ū_L = a²ωk cosh 2k(z₀ + H)/(2 sinh² kH) = … = **1.60** cm/s").
- **Explain:**
  0. *What the views show* — "**Water under the wave**: teal dots are water parcels, faint teal loops their orbits, the rose
     line a column of dye; blue is the surface. **Orbit size vs depth**: the half-widths A (solid) and half-heights B
     (dashed) of the loops; in exact mode the purple curve is the drift. **x of the selected particle** (hidden on phones):
     its to-and-fro motion, closed (teal) or creeping forward (purple)."
  1. *The wave* — "k = 2π/20 = 0.314 rad/m, H = kH/k = **H** m, a = ka/k = **a** m; (7.28) ω = √(gk tanh kH) = **ω** rad/s,
     T = **T** s."
  2. *The orbit at your depth* — "(7.36) A = a cosh k(z₀ + H)/sinh kH = **A** m, B = a sinh k(z₀ + H)/sinh kH = **B** m;
     B/A = tanh k(z₀ + H) = **BA**" boxed; "focal half-distance √(A² − B²) = a/sinh kH = **f** m — the same at every depth
     (cosh² − sinh² = 1)".
  3. *Deep and shallow shapes* — deep: "(7.46): circles of radius a e^(kz₀) = **r** m"; shallow: "(7.50): width a/kH =
     **w** m, height a(1 + z₀/H) = **h** m".
  4. *Sense of rotation* — "At the top of the loop (ζ = +B) the parcel moves forward (ξ̇ = +aω cosh/sinh > 0): clockwise
     for a wave travelling to the right."
  5. *The drift* — "(7.86) ū_L = a²ωk cosh 2k(z₀ + H)/(2 sinh² kH) = **uL** m/s = **uLcm** cm/s; per period ū_L T = **dp**
     m = **pct** % of the orbit width 2A." (deep: (7.85) a²ωk e^(2kz₀).)
  6. *The Eulerian mean* — "At a fixed point that stays under water the mean velocity is 0 (N92): the drift exists only
     following the particles."
  7. *Right now* — "t = **live t** (τ = **live tau**); the selected parcel is at ξ = **live xi** m, ζ = **live zeta** m;
     exact − linear x difference **live dx** m."
  8. *Reading the current setting* — deep linear: "Deep water, small amplitude: circles that shrink like e^(kz₀); at half a
     wavelength down only 4 % remains, so a diver there feels almost nothing." · shallow: "Shallow water: every parcel
     sloshes back and forth by the same a/kH, vertical motion dies linearly to zero at the bottom — the sea bed is swept to
     and fro (sediment transport)." · exact, ka > 0.15: "Finite amplitude: the loops do not close — the drift a²ωk e^(2kz₀)
     is **pct** % of the orbital speed and decays twice as fast as the orbits: only the top few metres are carried."
- **Derivation tab:** **D10** (9 steps) `view: 'column'`, goal `set` deep water linear; step 4 (freeze the argument) `watch`
  "the black arrow is the neglected ξ ∂u/∂x — one ka smaller than u"; step 7 **live** "ξ amplitude a cosh k(z₀ + H)/sinh kH
  = **A** m". **D11** (6 steps) `view: 'axes'` (phones `column`), step 4 **live** "B/A = tanh k(z₀ + H) = …", step 6 `watch`
  "the selected parcel runs clockwise". **D27** (9 steps) `view: 'track'` (phones `column`), goal `set` steep deep wave
  (exact); step 6 (average the two terms) **live** "⟨ξu_x + ζu_z⟩ = a²ωk e^(2kz₀) = **uL** m/s" with `watch` "the purple
  track creeps forward at this slope"; interpret "With your numbers the surface drift is **uLcm** cm/s, **n** km per day."
- **Code:**
  ```python
  k, kH, ka = 2*np.pi/20, {{kH}}, {{ka}}      # λ = 20 m
  H, a = kH/k, ka/k                           # depth {{H}} m, amplitude {{a}} m
  ax = ch07.orbit_semi_axes({{z0}}, a, k, H)  # (7.36) A = {{A}} m, B = {{B}} m
  T = 2*np.pi/ch07.omega_gravity(k, H)        # period {{T}} s
  p = ch07.particle_path(0.0, {{z0}}, np.linspace(0, 10*T, 2001),
                         a, k, H, model="{{model}}")   # (7.33) or (7.34)
  print(p["x"][-1] - p["x"][0])               # drift after 10 T: {{d}} m
  print(ch07.stokes_drift({{z0}}, a, k, H))   # (7.86) ū_L = {{uL}} m/s
  ```
- **Walkthrough (7 steps):** 1. "Does the water travel?" — "Floats bob as swell passes and come back. Watch the teal
  parcels: what path does each follow?" `set` deep water, `play: true` · 2. "Circles" — "In deep water each parcel goes
  round a circle of radius a e^(kz₀) (7.46): halved every 2.2 m here." `readouts: ['A']`, `derive: {id: 'D11', step: 4}` ·
  3. "Frozen argument" — "To first order the parcel feels the velocity at its mean position (7.34): the loops close."
  `derive: {id: 'D10', step: 4}`, `code: {lines: [3, 3]}` · 4. "Shallow water" — "Drag kH to 0.3: the circles squash into
  flat ellipses of equal width; the bottom parcel only slides (7.50)." `set` shallow, `controls: ['kH']` · 5. "Opening the
  loop" — "Switch to exact path lines at ka = 0.25: the loops no longer close and the dye line leans forward." `set` steep
  deep wave, `play: true` · 6. "The drift" — "The forward creep is ū_L = a²ωk e^(2kz₀) (7.85) — ka times the orbital speed,
  halved every 1.1 m." `derive: {id: 'D27', step: 6}`, `inspect: true` · 7. "Your turn" — "Predict the drift at z₀ = −2 m
  before clicking that parcel, then compare with the inspector." `controls: ['ka', 'kH']`.
- **Equations:** `path` "Path lines" ref 'Eq. (7.33)' · `orb` "Orbits" ref 'Eq. (7.36)' `\xi^2/A^2+\zeta^2/B^2=1` live A, B
  · `deep` "Deep circles" ref 'Eq. (7.46)' · `shal` "Shallow ellipses" ref 'Eq. (7.50)' · `psi` "Stream function" ref 'Eq.
  (7.37)' · `drift` "Stokes drift" ref 'Eq. (7.85)' `\bar u_L=a^2\omega ke^{2kz_0}` · `driftH` "Any depth" ref 'Eq. (7.86)'.
- **Check yourself:** (1) "In deep water, how deep must you go for the orbit to shrink to 1 % of its surface size?" —
  "e^(kz₀) = 0.01 ⇒ z₀ = −ln 100/k = −4.6/0.314 = −14.7 m (about 0.73 λ)." · (2) "In shallow water, which is larger at
  the bottom: the horizontal or the vertical excursion?" — "Horizontal: a/kH at every depth; the vertical one falls to 0
  at the bottom (7.50)." `set` shallow · (3) "Double ka in exact mode. By what factor does the drift grow?" — "4: ū_L ∝ a²
  (7.85)." `set` steep deep wave · (4) "Why does a current meter fixed 3 m down read zero mean while a float at 3 m
  drifts?" — "The meter averages whatever water passes it (N92: exactly 0 for an irrotational periodic wave); the float
  follows one parcel, which is further forward at the top of its loop."
- **Selftest parity rows** (inputs: the builder prints the preset's a, k, H and T to 10 significant digits into both the
  JS call and the `py` string — e.g. a = 0.3183098862, k = 0.3141592654, H = 9.549296586 for kH = 3, ka = 0.1):
  `{name: 'A at -2 m', js: axes(-2).A, py: 'ch07.orbit_semi_axes(-2.0, 0.3183098862, 0.3141592654, 9.549296586)["A"]',
  rtol: 1e-9}` · `{name: 'B at -2 m', js: axes(-2).B, py: '…["B"]', rtol:
  1e-9}` · `{name: 'drift surface deep', js: drift(0) (kH = 6), py: 'ch07.stokes_drift(0.0, 0.3183, 0.31416, 19.099)', rtol:
  1e-9}` · `{name: 'exact path x after 1 T', js: rk4Path(0, -1, 1).x, py: 'ch07.particle_path(0.0, -1.0, [0.0, T],
  a, k, H, model="exact")["x"][1]', rtol: 1e-6}` (ka = 0.25, kH = 4 → a = 0.7957747155, H = 12.73239545, T = 3.580… s
  printed to 10 digits) · `{name: 'orbit_state
  B/A', js: state(-3).B_over_A, py: 'ch07.orbit_state(-3.0, 0.3183, 0.31416, 3.183)["B_over_A"]', rtol: 1e-9}` · invariant
  `{name: 'focal constant', js: axes(-1).focal - axes(-5).focal, expect: 0, atol: 1e-12}`.
- **Fit plan:** 360×640: `column` (60 %) + `axes`; model chips and presets on one wrapping row; status one line. 844×390:
  `column` | `axes`. Desktop: rows [1.2, 1] with `axes` and `track` sharing row 1.

### E3 · capillary_gravity_waves
- **Title:** "Why is there a slowest ripple?" · **Summary:** "Sweep the wavelength over four decades and change the liquid:
  gravity's term g/k and surface tension's σk/ρ trade places, and between them the speed has a minimum — 23 cm/s at 1.7 cm
  on clean water." · **CORE:** C07 (also N50–N55 ((7.53)–(7.59)), N57 (7.60), R09 (1.5); links C09's c_g,min, N72) ·
  **Reference:** `overfitting_curves.html` (a minimal curve with an optimum marker and a regime verdict) with
  `fid_formula_lab.html`'s term bars.
- **meta:** `viz:order 3` · `viz:sections 7.3` · `viz:equations 7.53 7.54 7.55 7.56 7.57 7.58 7.60` · `viz:fluidpy
  ch07.capillary_state ch07.omega_capillary_gravity ch07.phase_speed ch07.capillary_minimum ch07.min_group_velocity
  ch07.capillary_surface_pressure` · `viz:derivations D15 D16`.
- **Physics:** `omegaCG(k, H, sig, rho)` ↔ `ch07.omega_capillary_gravity`; `c(k)` ↔ `ch07.phase_speed(k, H, sigma=, rho=)`;
  `cg(k)` ↔ `ch07.group_velocity(k, H, sigma=, rho=)`; `cmin(sig, rho)` = {c_min, lam_m, k_m} ↔ `ch07.capillary_minimum`;
  `cgmin(sig, rho)` (s = 2/√3 − 1) ↔ `ch07.min_group_velocity`; `crestP(a, k, sig)` = σak² ↔
  `ch07.capillary_surface_pressure(-a*k**2, sig)`; `state(lam)` ↔ `ch07.capillary_state(lam, sig, rho, H)`.
- **Modes (`liquid` chips):** water 20 °C (σ = 0.07274 N/m, ρ = 998.2) · soapy water (σ = 0.030, ρ = 1000) · mercury
  (σ = 0.485, ρ = 13 534) · ethanol (σ = 0.0223, ρ = 789) · custom (σ slider free).
- **Views** (rows [1, 1]):
  1. `patch` "The surface, magnified" (row 0, flex 1) — two wavelengths of η = a cos(kx − ωt) with a = λ/20 (drawn), a scale
     bar (1 mm … 1 m), at the crest two downward restoring arrows: gravity (blue, ∝ ρga) and surface tension (rose, ∝ σak²),
     their lengths in proportion; animated with the transport.
  2. `curve` "Phase speed c(λ)" (row 0, flex 1.3) — log–log, λ from 1 mm to 10 m, c from 5 cm/s to 5 m/s: c(λ) (black),
     capillary branch √(2πσ/ρλ) (rose dashed), deep gravity branch √(gλ/2π) (blue dashed), shallow plateau √(gH) (muted
     dotted), the minimum (black ring, labelled c_min, λ_m), the c_g curve (purple, thin) with its minimum (purple ring), the
     current dot (orange). Pointer: click → set λ.
  3. `bars` "Restoring terms" (row 1, `hidePortrait: true`) — term bars (the `terms` feature) g/k (blue) and σk/ρ (rose)
     summing to c²/tanh kH (black).
- **Controls:** `lam` "Wavelength $\lambda$" log 1 mm … 10 m, default 17.1 mm · `liquid` chips · `sig` "Surface tension
  $\sigma$" 0 … 0.5 N/m, step 0.001, default 0.0727 (optional; custom liquid) · `H` "Depth $H$" log 1 mm … 100 m, default
  1 m (optional).
- **Transport:** `tau` 0 → 3 periods, rate 0.5, loop.
- **Presets:** "λ = λ_m (the minimum)" {lam: lam_m(liquid)} · "raindrop ripple 5 mm" {lam: 0.005} · "pure capillary 1 mm"
  {lam: 0.001} · "wind wave 1 m" {lam: 1} · "σ = 0 (no tension)" {liquid: 'custom', sig: 0}.
- **Status:** "💧 capillary ripple: tension term 4.2× gravity's" (ratio > 2) · "⚖️ near the minimum: the two terms within a
  factor 2" · "🌊 gravity wave: tension only 1 % of the restoring force" (ratio < 0.5).
- **Readouts:** "Phase speed $c$" (cm/s) · "Minimum $c_{min}$" (cm/s) · "at $\lambda_m$" (cm) · "Crest pressure" (Pa, for
  a = λ/20).
- **Terms:** title "c² split (deep water)" unit m²/s²: g/k (blue) · σk/ρ (rose) · total (black).
- **Depth features:** Explain + Code + Derivation · **linked views** (3) · **term bars** · **presets** · **status** ·
  **modes** (liquids) · **notes** (a small table: liquid, σ, ρ, c_min, λ_m with the current liquid highlighted).
- **Explain:**
  0. *What the views show* — "**The surface, magnified**: the blue arrow at the crest is gravity's pull (the weight of the
     hump), the rose arrow the tension of the curved surface. **Phase speed c(λ)**: black = (7.57), rose dashed = pure
     capillary (7.60), blue dashed = pure gravity (7.45), purple = the group speed; the rings mark the minima. **Restoring
     terms** (hidden on phones; the same numbers are in section 2)."
  1. *Wavenumber* — "k = 2π/λ = **k** rad/m."
  2. *The two restoring terms* — "gravity g/k = 9.81/**k** = **G** m²/s²; tension σk/ρ = **sig** × **k**/**rho** = **S**
     m²/s²; ratio S/G = σk²/ρg = **r**" boxed.
  3. *Phase speed* — "(7.57) c = √((g/k + σk/ρ) tanh kH) = √((**G** + **S**) × **tanh**) = **c** m/s" boxed; "without σ:
     **c0** m/s — tension adds **pct** %."
  4. *The slowest wave for this liquid* — "(7.58) c_min = (4gσ/ρ)^(1/4) = (4 × 9.81 × **sig**/**rho**)^(1/4) = **cmin** m/s at
     λ_m = 2π√(σ/ρg) = **lm** cm" boxed.
  5. *Pressure under a crest* — "(7.54) p = −σ η_xx = σ a k² = **sig** × **a** × **k**² = **p** Pa (a = λ/20): the surface
     pushes the crest down."
  6. *Group speed* — "c_g = **cg** m/s (= c at λ_m: the tangent through the origin); its own minimum **cgmin** m/s at
     **lgm** cm sets the calm ring of a stone in a pond (C09)."
  7. *Right now* — "t = **live t**."
  8. *Reading the current setting* — capillary: "Surface tension wins: shorter ripples are faster, like waves on a drum
     skin; gravity hardly matters below about 4 mm." · crossover: "The two forces are comparable: you are near the slowest
     possible wave of this liquid — any other wavelength travels faster." · gravity: "Gravity wins: the tension correction
     is **pct** %; this is an ordinary water wave, longer is faster."
- **Derivation tab:** **D15** (10 steps) `view: 'patch'`, goal `set {lam: 0.005}`; step 3 (the sign under a crest) `watch`
  "the rose arrow points into the water: higher pressure below a crest"; step 8 **live** "η_xx = −k²η: σk²/ρ = **S·k** …";
  step 10 **live** "ω = √(k(g + σk²/ρ)) = **ω** rad/s". **D16** (7 steps) `view: 'curve'`, goal `set` λ = 1 cm; step 4 (set
  the slope to zero) `set {lam: lam_m}` with `watch` "the dot sits in the valley"; step 6 **live** "c_min² = 2√(gσ/ρ) =
  …"; interpret "For **liquid** the slowest wave is **cmin** m/s at **lm** cm."
- **Code:**
  ```python
  sigma, rho = {{sig}}, {{rho}}              # {{liquid}} [N/m], [kg/m^3]
  k = 2*np.pi/{{lam}}                         # wavenumber {{k}} rad/m
  s = ch07.capillary_state({{lam}}, sigma, rho, H={{H}})
  print(s["gravity_term"], s["tension_term"])   # g/k = {{G}}, σk/ρ = {{S}}
  print(s["c"], s["regime"])                  # (7.57) c = {{c}} m/s, {{regime}}
  m = ch07.capillary_minimum(sigma, rho)      # (7.58)
  print(m["c_min"], m["lam_m"])               # {{cmin}} m/s at {{lm}} m
  ```
- **Walkthrough (6 steps):** 1. "Racing ripples" — "Raindrop ripples outrun the slightly longer wind waves. Why would
  shorter waves be faster?" `set {lam: 0.005}`, `play: true` · 2. "Two forces" — "Gravity pulls a hump down with its
  weight; the curved surface pulls it flat like a drum skin (7.54): p = −σ ∂²η/∂x²." `terms: true`, `derive: {id: 'D15',
  step: 3}` · 3. "One new term" — "Everything of §7.2 carries over with g replaced by g + σk²/ρ (7.56)." `eq: 'disp'`,
  `derive: {id: 'D15', step: 9}` · 4. "The valley" — "g/k falls with k while σk/ρ grows: their sum has a minimum at
  λ_m = 2π√(σ/ρg) (7.58)." `set` λ = λ_m, `derive: {id: 'D16', step: 4}` · 5. "Other liquids" — "Mercury, soap, ethanol:
  the valley moves. Which liquid has the slowest ripple?" `controls: ['liquid']` · 6. "Your turn" — "Predict c for a 3 mm
  ripple on clean water before dragging λ there." `controls: ['lam']`.
- **Equations:** `pcap` ref 'Eq. (7.54)' `(p)_{z=\eta}=-\sigma\,\partial^2\eta/\partial x^2` · `disp` ref 'Eq. (7.56)'
  `\omega=\sqrt{k(g+\sigma k^2/\rho)\tanh kH}` · `cph` ref 'Eq. (7.57)' · `cmin` ref 'Eq. (7.58)'
  `c_{min}=[4g\sigma/\rho]^{1/4},\ \lambda_m=2\pi\sqrt{\sigma/\rho g}` · `pure` ref 'Eq. (7.60)' `c=\sqrt{2\pi\sigma/\rho\lambda}`.
- **Check yourself:** (1) "Halve σ (soap). How do λ_m and c_min change?" — "λ_m ∝ √σ → ×0.71; c_min ∝ σ^(1/4) → ×0.84
  (7.58)." · (2) "At λ_m the two bars are equal. Why?" — "d(c²)/dk = −g/k² + σ/ρ = 0 means g/k = σk/ρ (D16 step 4)." ·
  (3) "Which is faster on clean water: a 1 mm or a 5 cm wave?" — "The 1 mm ripple: 68 cm/s vs 30 cm/s — both sides of the
  valley rise." `set {lam: 0.001}` · (4) "Why does a stone leave a calm ring rather than a calm point?" — "Energy moves at
  c_g, whose minimum (≈ 18 cm/s) is below c_min: no wave can be slower than c_g,min, so the centre empties (C09)."
- **Selftest parity rows:** `{name: 'c at 1 cm water', js: c(2*Math.PI/0.01), py: 'ch07.phase_speed(2*np.pi/0.01, 1.0,
  sigma=0.07274, rho=998.2)', rtol: 1e-12}` · `{name: 'c_min water', js: cmin(0.07274, 998.2).c_min, py:
  'ch07.capillary_minimum(0.07274, 998.2)["c_min"]', rtol: 1e-12}` · `{name: 'lam_m mercury', js: cmin(0.485,
  13534).lam_m, py: 'ch07.capillary_minimum(0.485, 13534.0)["lam_m"]', rtol: 1e-12}` · `{name: 'cg_min water', js:
  cgmin(0.07274, 998.2).cg_min, py: 'ch07.min_group_velocity(0.07274, 998.2)["cg_min"]', rtol: 1e-10}` · `{name: 'crest
  pressure', js: crestP(0.001, 2*Math.PI/0.01, 0.07274), py: 'ch07.capillary_surface_pressure(-0.001*(2*np.pi/0.01)**2,
  0.07274)', rtol: 1e-12}` · `{name: 'tension term', js: state(0.005).tension_term, py: 'ch07.capillary_state(0.005,
  0.07274, 998.2, 1.0)["tension_term"]', rtol: 1e-12}` · invariant `{name: 'c=cg at minimum', js: cg(km) - c(km), expect:
  0, atol: 1e-9}`.
- **Fit plan:** 360×640: `curve` (60 %) above `patch`; bars hidden (numbers in Explain §2 and in the status); liquid chips
  wrap. 844×390: `patch` | `curve`. Desktop: rows [1, 1] with `bars` under `patch`.

### E4 · seiche_standing_waves
- **Title:** "Which waves can live in a lake?" · **Summary:** "Two equal waves travelling in opposite directions add into a
  pattern with fixed nodes; walls keep only the patterns that fit, so a basin rings at periods set by its length, depth
  and mode number." · **CORE:** C08 (also N58–N63 ((7.61)–(7.65)), N177, N178) · **Reference:**
  `angular_frequency_explorer_1.html` (system animation + a table of real cases with the current row highlighted).
- **meta:** `viz:order 4` · `viz:sections 7.4` · `viz:equations 7.61 7.62 7.63 7.64 7.65` · `viz:fluidpy
  ch07.standing_wave_fields ch07.seiche_modes ch07.seiche_state ch07.basin_modes` · `viz:derivations D17 D18`.
- **Physics:** `seiche(L, H, n)` ↔ `ch07.seiche_modes(L, H, n)` (k = (n + 1)π/L, ω from (7.65)); `state(L, H, n)` ↔
  `ch07.seiche_state(L, H, n)` (T_shallow = 2L/((n + 1)√(gH)), shallow_error); `fields(x, z, t)` ↔
  `ch07.standing_wave_fields(x, z, t, a, k, H)` (η = 2a cos kx cos ωt, (7.62), (7.63)) and the travelling parts from
  (7.2), (7.61).
- **Views** (rows [1.2, 1]):
  1. `basin` "The basin" (row 0, flex 1.5) — x ∈ [0, L], z ∈ [−H, 0] drawn with a vertical exaggeration printed in the
     title; walls (black), surface η (blue, exaggerated), streamlines of (7.62) (grey, 9 levels), velocity arrows (black,
     small), the right-going and left-going components as faint orange ghosts (when shown). Pointer: click → probe (x, z).
  2. `modes` "η and u along the basin" (row 1, flex 1) — η(x) envelope 2a cos kx (blue) and u(x) at the surface envelope ∝
     sin kx (black), nodes of η (blue rings) and of u (black rings), wall lines.
  3. `periods` "Period vs mode number" (row 1, flex 1, `hidePortrait: true`) — T_n (dots) for n = 0…6 and the shallow
     estimate 2L/((n + 1)√(gH)) (muted line), the current n ringed.
  Portrait: `basin` + `modes`; the period sits in the status.
- **Controls:** `L` "Basin length $L$" log 1 m … 500 km, default 1.5 m · `H` "Depth $H$" log 0.05 … 1000 m, default 0.2 m
  · `n` "Mode $n$" 0 … 5, step 1, default 0 · `comp` chips "right-going · left-going · both" (default both) · `a` optional.
- **Transport:** `tau` 0 → 4 periods, rate 0.5, loop.
- **Presets:** "bathtub" {L: 1.5, H: 0.2} · "swimming pool" {L: 25, H: 2} · "50 km lake" {L: 50000, H: 100} · "long
  shallow lake" {L: 388000, H: 19} (our inputs, Erie-like) · "harbour" {L: 500, H: 10}.
- **Status:** "🛁 T₀ = 2.20 s · shallow estimate 3 % short" · with comp ≠ both: "➡ one travelling wave: no nodes — the
  pattern moves".
- **Readouts:** "Wavelength $\lambda$" (m) · "Period $T$" (`fmtTime`) · "$kH$" · "Shallow error" (%).
- **Depth features:** Explain + Code + Derivation · **linked views** (3) · **transport** · **presets** · **inspector** (click
  the basin: "u = 2aω cosh k(z + H)/sinh kH sin kx sin ωt = … = **u** m/s; η = 2a cos kx cos ωt = …") · **notes** ("Right
  now": table bathtub / pool / harbour / 50 km lake / long lake with L, H, T₀ — current row highlighted) · **status**.
- **Explain:**
  0. *What the views show* — "**The basin**: blue surface (vertically exaggerated ×**ex**), grey streamlines of (7.62), black
     velocity arrows, faint orange the two travelling waves. **η and u**: the envelopes of height (blue) and surface
     velocity (black) — where one has a node the other has its maximum. **Period vs mode** (hidden on phones)."
  1. *Which wavelengths fit* — "u = 0 at both walls needs sin kL = 0: kL = (n + 1)π (7.64), λ = 2L/(n + 1) = 2 × **L**/**n1**
     = **lam** m" boxed.
  2. *The frequency* — "(7.65) ω = √((πg(n + 1)/L) tanh((n + 1)πH/L)) = √(**A** × tanh **kH**) = **w** rad/s; T = 2π/ω =
     **T**" boxed.
  3. *Shallow or not* — "kH = **kH**; shallow estimate T = 2L/((n + 1)√(gH)) = **Ts** (error **err** %)."
  4. *Why it stands still* — "(7.2) + (7.61): a cos(kx − ωt) + a cos(kx + ωt) = 2a cos kx cos ωt: the shape cos kx never
     moves."
  5. *Right now* — "t = **live t**, η at the left wall = **live eta0**, largest surface u = **live umax** m/s."
  6. *Reading the current setting* — "Mode **n** has **n+1** half-wavelengths in the basin: the water swings from one end to
     the other (n = 0) or in **n+1** opposite cells. " + shallow case: "Here kH is small, so T = 2L/((n + 1)√(gH)): long
     basins ring slowly — tens of minutes to hours for real lakes." · deep case: "Here the basin is deep compared with
     its length (kH > 2): the depth no longer matters, T depends on L alone."
- **Derivation tab:** **D17** (8 steps) `view: 'basin'`, goal `set {comp: 'right-going'}`; step 1 `set {comp: 'left-going'}`;
  step 3 `set {comp: 'both'}` with `watch` "the nodes stop moving"; step 7 **live** "u amplitude 2aω cosh kH/sinh kH = …".
  **D18** (5 steps) `view: 'modes'`, step 3 `set {n: 1}` with `watch` "two cells, λ = L"; step 5 **live** "ω = √(πg(n+1)/L
  tanh((n+1)πH/L)) = …"; interpret "Your basin rings at **T** in mode **n**."
- **Code:**
  ```python
  L, H, n = {{L}}, {{H}}, {{n}}              # basin length, depth [m], mode
  m = ch07.seiche_modes(L, H, n)             # (7.64)-(7.65)
  print(m["lam"], m["T"])                    # λ = {{lam}} m, T = {{T}} s
  s = ch07.seiche_state(L, H, n)
  print(s["T_shallow"], s["shallow_error"])  # 2L/((n+1)√(gH)) = {{Ts}} s, {{err}}
  f = ch07.standing_wave_fields({{x}}, {{z}}, {{t}}, 0.01, m["k"], H)
  print(f["eta"], f["u"])                    # η, u at the probe
  ```
- **Walkthrough (6 steps):** 1. "The slosh" — "Push the water in a bathtub and it rocks at a fixed rhythm. Where does
  that period come from?" `set` bathtub, `play: true` · 2. "One wave" — "A single wave travels: its crests march along."
  `set {comp: 'right-going'}` · 3. "Two waves" — "Add the same wave going left: cos(kx − ωt) + cos(kx + ωt) = 2 cos kx cos ωt
  — fixed nodes." `set {comp: 'both'}`, `derive: {id: 'D17', step: 5}` · 4. "Walls choose" — "A wall needs u = 0: only
  kL = (n + 1)π fits (7.64)." `derive: {id: 'D18', step: 3}`, `code: {lines: [2, 3]}` · 5. "A real lake" — "A 50 km lake,
  100 m deep: shallow, T₀ = 2L/√(gH) ≈ 53 min. Step n up." `set` 50 km lake, `controls: ['n']`, `notes: true` · 6. "Your
  turn" — "Predict the pool's lowest period, then check with the preset." `controls: ['L', 'H']`.
- **Equations:** `left` ref 'Eq. (7.61)' `\eta=a\cos[kx+\omega t]` · `stand` (unnumbered) `\eta=2a\cos kx\cos\omega t` ·
  `psi` ref 'Eq. (7.62)' · `u` ref 'Eq. (7.63)' · `lam` ref 'Eq. (7.64)' `\lambda=2L/(n+1)` · `freq` ref 'Eq. (7.65)'.
- **Check yourself:** (1) "Double L of a shallow lake. What happens to T₀?" — "It doubles: T = 2L/√(gH)." · (2) "At a node
  of η, what does the water do?" — "It moves fastest horizontally: u ∝ sin kx is largest where cos kx = 0 (7.63)." · (3)
  "Switch off one component. Why does the pattern move?" — "A single travelling wave has no fixed nodes; standing needs
  two equal opposite waves." · (4) "When does the depth stop mattering?" — "When the basin mode is 'deep' (kH > 2): tanh
  → 1 and ω = √(πg(n + 1)/L)." `set {L: 1, H: 1}`.
- **Selftest parity rows:** `{name: 'T bathtub', js: seiche(1.5, 0.2, 0).T, py: 'ch07.seiche_modes(1.5, 0.2, 0)["T"]', rtol:
  1e-12}` · `{name: 'T lake n=2', js: seiche(50000, 100, 2).T, py: 'ch07.seiche_modes(50000.0, 100.0, 2)["T"]', rtol:
  1e-12}` · `{name: 'shallow error pool', js: state(25, 2, 0).shallow_error, py: 'ch07.seiche_state(25.0, 2.0,
  0)["shallow_error"]', rtol: 1e-9}` · `{name: 'u probe', js: fields(0.3, -0.1, 0.4).u, py: 'ch07.standing_wave_fields(0.3,
  -0.1, 0.4, 0.01, 2.0943951, 0.2)["u"]', rtol: 1e-9}` · invariant `{name: 'node of eta', js: fields(L/(2*(n+1)), 0,
  0.3).eta, expect: 0, atol: 1e-12}`.
- **Fit plan:** 360×640: `basin` (55 %) + `modes`; comp chips and presets one wrapping row. 844×390: `basin` | `modes`.
  Desktop: rows [1.2, 1].

### E5 · group_velocity_packets
- **Title:** "Where does a wave's energy go?" · **Summary:** "Follow one crest while the envelope — and the energy —
  moves at c_g = dω/dk: half the crest speed in deep water, the same in shallow water, faster for capillary ripples; the
  chord and the tangent of ω(k) move with the dot." · **CORE:** C09 (also C06's E (7.42) and F (7.44), N65 (7.66), N66,
  N67 (7.68), N68–N70 ((7.69)–(7.71)), N72, N73; links C16's c ⟂ c_g) · **Reference:**
  `amplitude_phase_second_order_II_3.html` (system + graph windows linked by one state; a numbered derivation with live
  numbers in the explanation).
- **meta:** `viz:order 5` · `viz:sections 7.5` · `viz:equations 7.42 7.66 7.67 7.68 7.69 7.70 7.71` · `viz:fluidpy
  ch07.beat_wave ch07.gaussian_packet ch07.group_velocity ch07.packet_state ch07.pond_ripples ch07.min_group_velocity
  ch07.wave_energy_density` · `viz:derivations D19 D20 D21`.
- **Physics:** `omega(k)` for the chosen dispersion (deep √(gk); finite √(gk tanh kH); shallow k√(gH); capillary
  √(σk³/ρ)) ↔ `ch07.omega_capillary_gravity(k, H, sigma, rho, g)` with the matching data (capillary: g = 0); `cg(k)` ↔
  `ch07.group_velocity(k, H, g=, sigma=, rho=)`; `beat(x, t)` = 2a cos(½Δk x − ½Δω t) cos(k̄x − ω̄t) (and the printed
  ½Δω x ghost) ↔ `ch07.beat_wave(x, t, k1, k2, …, printed=)`; `packet(x, t, order)` = closed-form chirped Gaussian
  (order 2) or a(x − c_g t) cos(k₀x − ω₀t) (order 1) ↔ `ch07.gaussian_packet(x, t, a, k0, sigma_x, order=, H=, g=, sigma=,
  rho=)`; `pond(x, t)` = direct cosine sum of 256 modes ↔ `ch07.pond_ripples(x, t, width, n_modes=256, …)`; `pstate(k0)` ↔
  `ch07.packet_state(k0, H, g, sigma, rho, a)`.
- **Modes (`mode` chips):** "two waves (beats)" · "packet" · "stone in a pond".
- **Views** (rows [1.25, 1]):
  1. `wave` "η(x, t)" (row 0, flex 1.6) — x window following the group (width 12 group lengths in packet mode, 2 beat
     lengths in beats mode, 0…1 m radius in pond mode); the carrier (blue), the envelope ± (purple), one tracked crest
     (orange dot, re-seeded at the rear when it dies at the front — N73), a node or the envelope peak (purple dot); speed
     arrows c (orange) and c_g (purple) in the corner; in beats mode the printed-slip envelope (rose dashed) when the "show
     the book's printed form" toggle is on; in pond mode the calm-centre radius c_g,min t (purple dashed) and the front
     c_g,max t.
  2. `disp` "ω(k)" (row 1, flex 1) — the dispersion curve (black) around k₀ (range 0…3k₀), the chord from the origin
     (orange, slope c) and the tangent (purple, slope c_g) at the dot; in beats mode the two points (k₁, ω₁), (k₂, ω₂) and
     their chord Δω/Δk (purple dashed).
  3. `energy` "Energy density and flux" (row 1, flex 1, `hidePortrait: true`) — ½ρga(x)² along x (purple fill) and the flux
     arrow F = E c_g with its value; in pond mode the spectrum |A(k)| instead.
  Portrait: `wave` + `disp`; c, c_g and their ratio in the `wave` title.
- **Controls:** `mode` chips · `disp` chips "deep · finite H · shallow · capillary" · `lam0` "Carrier wavelength $\lambda_0$"
  log 1 mm … 500 m, default 156 m · `width` "Group length" 1 … 12 wavelengths, default 4 (packet: σ_x = width·λ₀/4; beats:
  Δk = k₀/width) · `H` "Depth $H$" log 0.5 … 500 m, default 10 m (optional; finite H and shallow) · `printed` "Show the
  printed ½Δω x" toggle (optional, beats mode) · `order` "Keep ω″ (spreading)" toggle (optional, packet mode; off = D20's
  first-order envelope).
- **Transport:** `tau` 0 → 20 carrier periods, rate 1, `end: 'hold'` with an end-of-run card: "In 20 T = **t** the crest
  travelled **xc** m, the envelope **xe** m (ratio **r** = c_g/c)."
- **Presets:** "swell set" {mode: 'packet', disp: 'deep', lam0: 156, width: 4} · "shallow (no dispersion)" {disp: 'shallow',
  lam0: 200, H: 5} · "pure capillary" {disp: 'capillary', lam0: 0.003} · "narrow packet (spreads)" {width: 1.5, order: true}
  · "wide packet" {width: 10} · "beats" {mode: 'beats', disp: 'deep'} · "stone in a pond" {mode: 'stone in a pond'}.
- **Status:** "🟠 crests overtake the group: c_g = 0.50 c — crests are born at the rear" · "🟣 together: c_g = c (shallow,
  no dispersion)" · "💧 the group overtakes the crests: c_g = 1.50 c — crests are born at the front".
- **Readouts:** "Phase speed $c$" · "Group speed $c_g$" · "$c_g/c$" · "Energy $E$" (J/m², a = 1 m) · "Flux $F$" (kW/m).
- **Depth features:** Explain + Code + Derivation · **linked views** (3) · **transport** + end-of-run card · **modes** ·
  **presets** · **status** · **notes** ("Right now": the crest you follow is **n** wavelengths from the rear of the group;
  it will die in **t** s).
- **Explain:**
  0. *What the views show* — "**η(x, t)**: blue carrier, purple envelope, the orange dot one crest, the purple dot the
     group's peak (or a node). **ω(k)**: orange chord = phase speed, purple tangent = group speed. **Energy** (hidden on
     phones): ½ρga² along the group and the flux F."
  1. *Carrier* — "k₀ = 2π/λ₀ = **k** rad/m; ω₀ = **w** rad/s from the chosen relation; c = ω₀/k₀ = **c** m/s" boxed.
  2. *Group speed* — "(7.67) c_g = dω/dk at k₀ = **cg** m/s" boxed; deep: "(7.70) c_g = c/2"; finite H: "(7.69) c_g = (c/2)(1 +
     2kH/sinh 2kH) = (c/2)(1 + **q**) = **cg**"; capillary: "ω = √(σk³/ρ) ⇒ c_g = 3c/2".
  3. *Beats (beats mode)* — "(7.66) with Δk = **dk**, Δω = **dw**: node speed Δω/Δk = **cgf** m/s → dω/dk as Δk → 0; the
     printed ½Δω x would give a frozen envelope."
  4. *Crests through the group* — "relative speed c − c_g = **dc** m/s; a crest crosses a group of length **Lg** m in **tx**
     s."
  5. *Energy and flux* — "(7.42) E = ½ρga² = **E** J/m² (a = 1 m); (7.71) F = E c_g = **F** W per metre of crest — not E c =
     **Ec**."
  6. *Arrival from a storm* — "over 1000 km the energy needs 10⁶/c_g = **tarr** (the crests would suggest **tc**)."
  7. *Spreading (packet mode)* — "the neglected ½ω″(k − k₀)² term matters after t ≈ 1/(|ω″|δk²) = **tspread** s; now t =
     **live t** s."
  8. *Stone in a pond (pond mode)* — "c_g,min = **cgmin** m/s at λ = **lgm** cm (Exercise 7.10 value): inside r = c_g,min t =
     **live rcalm** m the water is calm."
  9. *Reading the current setting* — deep: "Deep water: the group moves at half the crest speed. Follow a crest: it appears
     at the back, runs through, grows longer and dies at the front. Swell 'sets' arrive at c_g." · shallow: "No dispersion:
     every component moves at √(gH), the packet keeps its shape and c_g = c — a tsunami's energy arrives with its crests."
     · capillary: "Ripples: c_g > c, crests appear at the front and slide backwards through the group." · pond: "All
     wavelengths start together; each leaves at its own c_g: long waves lead, the shortest ripples trail, and nothing is
     slower than c_g,min."
- **Derivation tab:** **D19** (7 steps) `view: 'wave'`, goal `set {mode: 'beats'}`; step 4 `set {printed: true}` with `watch`
  "the rose envelope stays put — the printed form cannot be right"; step 6 **live** "Δω/Δk = **cgf** m/s vs dω/dk = **cg**".
  **D20** (12 steps) `view: 'wave'`, goal `set` swell set; step 5 (Taylor-expand ω) `view: 'disp'` with `watch` "the purple
  tangent is the linear term"; step 9 `set {order: false}` "the envelope keeps its shape"; step 11 `set {order: true, width:
  1.5}` `watch` "now it spreads: the dropped term"; step 12 **live** "envelope peak at c_g t = **xe** m". **D21** (11 steps)
  `view: 'disp'`, goal `set {disp: 'finite H', H: 10, lam0: 50}`; step 7 **live** "c_g/c = ½(1 + 2kH/sinh 2kH) = **ratio**";
  step 10 **live** "F = E c_g = **F** W/m".
- **Code:**
  ```python
  k0 = 2*np.pi/{{lam0}}                      # carrier wavenumber {{k}} rad/m
  disp = dict(H={{H}}, g={{g}}, sigma={{sig}}, rho=1000.0)   # {{disp}}
  c = ch07.phase_speed(k0, **disp)           # crest speed {{c}} m/s
  cg = ch07.group_velocity(k0, **disp)       # (7.67) dω/dk = {{cg}} m/s
  p = ch07.gaussian_packet({{x}}, {{t}}, 1.0, k0, {{sx}}, order={{order}}, **disp)
  print(p["envelope"], p["cg"])              # envelope at the probe
  print(ch07.wave_energy_density(1.0)*cg)    # (7.71) F = E c_g = {{F}} W/m
  ```
- **Walkthrough (7 steps):** 1. "Sets" — "Swell arrives in groups. Watch the orange crest: does it stay in its group?" `set`
  swell set, `play: true` · 2. "Two speeds" — "The crest runs at c = 15.6 m/s, the group at c_g = 7.8 m/s: crests are born
  at the back and die at the front." `readouts: ['c', 'cg']` · 3. "Beats" — "Two close waves add to (7.66): a fast carrier
  times a slow envelope moving at Δω/Δk." `set {mode: 'beats'}`, `derive: {id: 'D19', step: 6}` · 4. "Tangent vs chord" —
  "Crests move at the chord's slope ω/k, groups at the tangent's slope dω/dk (7.67)." highlight `view:disp` · 5. "Why the
  envelope moves at c_g" — "Taylor-expand ω about k₀: the linear term shifts the envelope by c_g t (7.68)." `derive: {id:
  'D20', step: 9}` · 6. "Energy" — "The energy flux is E c_g (7.71): energy rides with the group, not the crests." `derive:
  {id: 'D21', step: 10}`, `code: {lines: [4, 7]}` · 7. "Your turn" — "Predict: in shallow water, does the packet spread?
  Choose the preset and check." `controls: ['disp', 'lam0']`.
- **Equations:** `beats` ref 'Eq. (7.66)' · `cg` ref 'Eq. (7.67)' `c_g=\Delta\omega/\Delta k=d\omega/dk` · `pack` ref 'Eq.
  (7.68)' `\eta=a(x-c_gt)\cos(kx-\omega t)` · `cgw` ref 'Eq. (7.69)' · `lim` ref 'Eq. (7.70)' · `flux` ref 'Eq. (7.71)'
  `F=Ec_g` · `E` ref 'Eq. (7.42)' `E=\tfrac12\rho ga^2`.
- **Check yourself:** (1) "In deep water, how many crests pass through a group while it travels its own length?" — "The
  crests move c − c_g = c/2 faster, so in the time the group moves L, crests move 2L: they gain L — as many crests as the
  group holds." · (2) "Why does the shallow packet not spread?" — "ω = k√(gH) is a straight line: ω″ = 0 and c_g = c for
  all k." `set` shallow · (3) "Which arrives first from a storm: 10 s or 16 s swell?" — "16 s: c_g = gT/4π grows with T."
  · (4) "In the pond mode, why is there a calm circle?" — "c_g has a minimum ≈ 18 cm/s: no energy is slower, so the
  region inside c_g,min t empties."
- **Selftest parity rows:** `{name: 'cg deep lam156', js: cg(2*Math.PI/156), py: 'ch07.group_velocity(2*np.pi/156.0)', rtol:
  1e-12}` · `{name: 'cg H=10 lam50', js: cgFinite(2*Math.PI/50, 10), py: 'ch07.group_velocity(2*np.pi/50.0, 10.0)', rtol:
  1e-12}` · `{name: 'cg capillary ratio', js: cgCap(2*Math.PI/0.003)/cCap(2*Math.PI/0.003), py:
  'ch07.group_velocity(2*np.pi/0.003, np.inf, g=0.0, sigma=0.0727)/ch07.phase_speed(2*np.pi/0.003, np.inf, g=0.0,
  sigma=0.0727)', rtol: 1e-12}` · `{name: 'beat node speed', js: beatState(0.9, 1.1).cg_finite, py: 'ch07.beat_wave(0.0, 0.0,
  0.9, 1.1)["cg_finite"]', rtol: 1e-12}` · `{name: 'packet eta', js: packet(620, 200, 2).eta, py: 'ch07.gaussian_packet(620.0,
  200.0, 1.0, 0.0402768, 156.0, order=2)["eta"]', rtol: 1e-8}` · `{name: 'pond eta', js: pond(0.3, 1.0), py:
  'ch07.pond_ripples(0.3, 1.0, width=0.01, n_modes=256)', rtol: 1e-9}` · `{name: 'arrival h', js: pstate(2*Math.PI/156).arrival_h,
  py: 'ch07.packet_state(2*np.pi/156.0)["arrival_h"]', rtol: 1e-12}` · invariant `{name: 'F = E cg', js: F - E*cg, expect: 0,
  atol: 1e-9}`.
- **Fit plan:** 360×640: `wave` (60 %) + `disp`; mode and dispersion chips on one wrapping row (dispersion chips optional
  on phones, presets carry them). 844×390: `wave` | `disp`. Desktop: rows [1.25, 1], `energy` beside `disp`.

### E6 · wave_rays_refraction
- **Title:** "Why do waves arrive parallel to the beach?" · **Summary:** "Along a ray the frequency never changes while the
  wavelength shrinks with the depth, so the shoreward end of a crest slows and the crest swings round — over a beach, a
  ridge or an island." · **CORE:** C10 (also N47, N48, N74–N80 ((7.72)–(7.78)), N174, N175, N183, N184) · **Reference:**
  `stride_padding_playground.html` (classic-setting presets and the formula with numbers plugged in) with
  `angular_frequency_explorer_1.html`'s linked views.
- **meta:** `viz:order 6` · `viz:sections 7.2 7.5` · `viz:equations 7.28 7.73 7.74 7.75 7.79` · `viz:fluidpy
  ch07.snell_ray_plane_beach ch07.refraction_state ch07.ray_trace ch07.omega_gravity ch07.group_velocity` ·
  `viz:derivations D22 D23 D24`.
- **Physics:** `omega(kvec, x)` = √(g|k| tanh(|k|H(x))) with bathymetries: plane beach H = s·x (x offshore), ridge H =
  H₀ − h e^{−(y/w)²}, island H = min(H₀, s(r − R)) (r from the island centre), homogeneous H = H₀; `ray(x0, k0)` = RK4
  (`Viz.num.rk4Step`, 600 steps) on dx/dt = ∇_k ω, dk/dt = −∇_x ω (central differences, h = 1e-6 relative) ↔
  `ch07.ray_trace(None, x0, k0, t_span, H_fn=…)` (parity on ω drift and on the plane-beach angle); `snell(x)` (closed form:
  k(x) from `kFromOmega` (brentq) and α from k sin α = const) ↔ `ch07.snell_ray_plane_beach(x, alpha0, x0, T, slope)`;
  `stateAt(H)` ↔ `ch07.refraction_state(T, alpha0, H0, H)`.
- **Modes (`bathy` chips):** "plane beach" · "circular island" · "submarine ridge" · "homogeneous".
- **Views** (rows [1.3, 1]):
  1. `plan` "Plan view" (row 0, flex 1.5, equal aspect) — depth contours (grey, labelled in m), 9 rays (purple) launched from
     the offshore edge at α₀ to the x-axis normal, crest lines every wavelength drawn perpendicular to the local k (orange,
     thin), the shoreline or island (sand); the selected ray bold with a moving dot (transport). Pointer: click a ray to
     select it; click a point on it for the inspector.
  2. `along` "Along the selected ray" (row 1, flex 1) — against distance travelled: ω (black, flat — the conservation), k
     (orange), c (orange dashed), c_g (purple), α (teal, right axis in degrees); a vertical line at the moving dot.
  3. `xt` "x–t diagram" (row 1, flex 1, `hidePortrait: true`) — for the homogeneous chip: crest lines (thin orange, slope c)
     and rays (thick purple, slope c_g) (Fig. 7.17); otherwise the selected ray's x(t) (curved, Fig. 7.18) with ω printed.
  Portrait: `plan` + `along`; α at the shore and ω drift in the status.
- **Controls:** `bathy` chips · `alpha0` "Incidence $\alpha_0$" 0 … 70°, default 30° · `T` "Period $T$" 4 … 20 s, default
  8 s · `slope` "Beach slope" 1:10 … 1:200 (log, optional), default 1:50 · `H0` "Offshore depth" 5 … 200 m (optional),
  default 20 m.
- **Transport:** `s` 0 → 1 (fraction of the selected ray travelled), rate 0.1/s, `end: 'hold'`; the dot moves along the ray
  and the `along` cursor follows; real travel time printed via `fmtTime`.
- **Presets:** "swell at 30° on a 1:50 beach" {bathy: 'plane beach', alpha0: 30, T: 8, slope: 0.02, H0: 20} · "long period
  (turns sooner)" {T: 16} · "island" {bathy: 'circular island', alpha0: 0, T: 10} · "ridge (focusing)" {bathy: 'submarine
  ridge', alpha0: 0} · "homogeneous (straight rays)" {bathy: 'homogeneous'}.
- **Status:** "🧭 ω conserved to 3×10⁻¹⁰ along the ray · α 30.0° → 11.3° at the 2 m contour" (plane beach) · "🏝️ 4 of 9 rays
  reach the shadow side" (island) · "➖ uniform depth: straight rays" (homogeneous).
- **Readouts:** "$k$ at the dot" (rad/m) · "$\alpha$ at the dot" (°) · "$c$, $c_g$" (m/s) · "ω drift" (relative).
- **Depth features:** Explain + Code + Derivation · **linked views** (3) · **modes** (bathymetries) · **presets** ·
  **status** · **inspector** (click a ray point: "ω = 2π/T = 0.785 rad/s; solve ω² = gk tanh kH at H = 5 m → k = 0.118
  rad/m; Snell k sin α = k₀ sin α₀ = 0.0708 × 0.5 = 0.0354 ⇒ sin α = 0.0354/0.118 = 0.300 ⇒ α = **17.4°**") · **transport**.
- **Explain:**
  0. *What the views show* — "**Plan view**: grey depth contours, purple rays (paths of energy, moving at c_g), orange crest
     lines (perpendicular to k). **Along the selected ray**: black ω is flat — that is the law; k, c, c_g and α change.
     **x–t** (hidden on phones): crest lines vs rays."
  1. *Frequency is fixed* — "ω = 2π/T = **w** rad/s along every ray ((7.79): ∂ω/∂t + c_g ∂ω/∂x = 0 in a steady medium)."
  2. *Wavenumber from the depth* — "at the start H₀ = **H0** m: solve (7.28) ω² = gk tanh kH₀ → k₀ = **k0** rad/m (λ = **l0**
     m); at the dot H = **H** m → k = **k** rad/m (λ = **l** m)" boxed.
  3. *Snell* — "straight contours: the along-shore wavenumber k sin α is conserved (D24): sin α = k₀ sin α₀/k = **k0** ×
     **s0**/**k** = **sa** ⇒ α = **a** °" boxed.
  4. *Speeds* — "c = ω/k = **c** m/s, c_g = **cg** m/s (7.69); travel time so far **tt**."
  5. *Conservation check* — "ω drift along the traced ray: **drift** (RK4 in your browser)."
  6. *Island / ridge* — "(mode notes) rays bend toward the shallower region: a ridge focuses energy (bigger waves above it);
     an island's slope turns rays onto its lee side."
  7. *Reading the current setting* — plane beach: "The shoreward end of each crest is in shallower, slower water, so the
     crest pivots: by the 2 m contour it is within **a**° of the shoreline, whatever the offshore angle was." · long period:
     "Longer waves feel the bottom farther offshore (kH smaller), so they turn sooner and arrive straighter." ·
     homogeneous: "Nothing changes along the path: k, ω, c and c_g are all constant, rays are straight lines."
- **Derivation tab:** **D22** (8 steps) `view: 'xt'` (phones `along`), goal `set {bathy: 'homogeneous'}`; step 5 `watch`
  "crest lines end at the group's front"; step 8 **live** "c_g = **cg** m/s is the slope of the thick lines". **D23** (8
  steps) `view: 'along'`, goal `set` swell at 30°; step 6 `watch` "the black ω line stays flat"; step 8 **live** "ω drift
  **drift**". **D24** (9 steps) `view: 'plan'`, goal `set` swell at 30°; step 5 `watch` "the teal α curve falls"; step 8
  **live** "sin α = **k0** × 0.5/**k** = **sa**"; interpret "With T = **T** s the crest turns from **a0**° to **a**° by
  the **H** m contour."
- **Code:**
  ```python
  T, a0, H0, s = {{T}}, np.radians({{alpha0}}), {{H0}}, {{slope}}
  x = np.linspace(H0/s, 0.05/s, 200)          # offshore distance [m]
  r = ch07.snell_ray_plane_beach(x, a0, H0/s, T, s)
  print(np.degrees(r["alpha"][[0, -1]]))      # {{a0}}° → {{aend}}°
  print(r["snell"].std())                     # k sin α constant: {{sd}}
  st = ch07.refraction_state(T, a0, H0, {{H}})
  print(st["k"], st["alpha_deg"], st["cg"])   # at H = {{H}} m
  ```
- **Walkthrough (6 steps):** 1. "The puzzle" — "Waves reach the beach almost parallel to it on any day. Why?" `set` swell
  at 30°, `play: true` · 2. "Crests are counted" — "Crests are neither made nor destroyed: ∂k/∂t + ∂ω/∂x = 0 (7.74)."
  `derive: {id: 'D22', step: 4}` · 3. "Frequency rides the ray" — "In a steady medium ω is carried unchanged at c_g (7.79):
  the black line is flat." highlight `view:along`, `derive: {id: 'D23', step: 6}` · 4. "The crest pivots" — "k grows as H
  falls; along-shore k sin α is fixed, so α shrinks (Snell)." `derive: {id: 'D24', step: 8}`, `inspect: true` · 5.
  "Islands" — "On an island's slope the rays curl round to the lee side." `set` island · 6. "Your turn" — "Predict: does a
  16 s swell arrive straighter than an 8 s swell? Try it." `controls: ['T', 'alpha0']`.
- **Equations:** `loc` ref 'Eq. (7.73)' `k\equiv\partial\theta/\partial x,\ \omega\equiv-\partial\theta/\partial t` · `crest` ref
  'Eq. (7.74)' · `kadv` ref 'Eq. (7.75)' · `wray` ref 'Eq. (7.79)' `\partial\omega/\partial t+c_g\partial\omega/\partial x=0` ·
  `snell` (ours) `k\sin\alpha=\text{const}` · `disp` ref 'Eq. (7.28)'.
- **Check yourself:** (1) "What stays constant along a ray: k, ω or c?" — "ω (7.79); k and c change with depth." · (2)
  "Start at α₀ = 0. Do the rays bend on a plane beach?" — "No: k sin α = 0 stays 0; they go straight in." `set {alpha0: 0}`
  · (3) "Why does a ridge make bigger waves above it?" — "Rays bend toward the shallower ridge and crowd together: more
  energy flux per metre of crest." `set` ridge · (4) "In the x–t diagram, which lines are steeper — crest lines or rays?"
  — "Rays are steeper in x–t (slower, c_g < c) in deep water; crest lines cross them."
- **Selftest parity rows:** `{name: 'alpha at 2 m', js: stateAt(2).alpha_deg, py: 'ch07.refraction_state(8.0,
  0.5235987756, 20.0, 2.0)["alpha_deg"]', rtol: 1e-9}` · `{name: 'k at 5 m', js: stateAt(5).k, py:
  'ch07.refraction_state(8.0, 0.5235987756, 20.0, 5.0)["k"]', rtol: 1e-9}` · `{name: 'snell ray alpha', js:
  snell([1000, 100]).alpha[1], py: 'ch07.snell_ray_plane_beach([1000.0, 100.0], 0.5235987756, 1000.0, 8.0,
  0.02)["alpha"][1]', rtol: 1e-9}` · `{name: 'RK4 ray angle vs Snell', js: rayAngleAt(100), py: 'ch07.snell_ray_plane_beach([1000.0,
  100.0], 0.5235987756, 1000.0, 8.0, 0.02)["alpha"][1]', rtol: 1e-5}` · invariant `{name: 'omega drift', js: rayDrift(),
  expect: 0, atol: 1e-7}`.
- **Fit plan:** 360×640: `plan` (60 %) + `along`; bathymetry chips one row; slope/H0 optional. 844×390: `plan` | `along`.
  Desktop: rows [1.3, 1], `xt` beside `along`.

### E7 · hydraulic_jump
- **Title:** "How high does the water jump — and why only upward?" · **Summary:** "A fast shallow stream becomes a slow
  deep one with the same momentum flux: the Fr₁ slider sets H₂ from the momentum balance, and the energy bar shows the
  loss — which would turn into a gain, forbidden, for Fr₁ < 1." · **CORE:** C11 (also N81–N85, N93, N95, N97 ((7.87),
  (7.88)), R10 (4.104); links ch04 `control_volume_budgets`) · **Reference:** `fid_formula_lab.html` (term bars with a
  total) with ch04 `control_volume_budgets`' face-flux bars.
- **meta:** `viz:order 7` · `viz:sections 7.6` · `viz:equations 4.104 7.80 7.81 7.87 7.88` · `viz:fluidpy ch07.jump_state
  ch07.hydraulic_jump ch07.jump_momentum_residual ch07.ursell_number ch07.solitary_wave` · `viz:derivations D25 D26`.
- **Physics:** `jump(H1, Fr1)` = {ratio = ½(−1 + √(1 + 8Fr₁²)), H2, u1, u2, Q, Fr2, dE = −g(H₂ − H₁)³/(4H₁H₂), head_loss,
  power, mom_in, mom_out, p_in, p_out, residual, allowed, bore_speed} ↔ `ch07.jump_state(H1, Fr1, g, rho, frame)`;
  `ursell(a, lam, H)` ↔ `ch07.ursell_number`; `sol(x, t, a, H)` ↔ `ch07.solitary_wave(x, t, a, H)`.
- **Modes (`frame` chips):** "stationary jump" (the CV at rest; u₁ → from the left) · "bore into still water" (the same
  jump seen by an observer at rest: water ahead at rest, the step moving at the bore speed). **Fate chips** (`fate`):
  "jump" · "solitary wave" (replaces the channel view by a sech² hump on depth H with a/H slider and the Ursell number).
- **Views** (rows [1.2, 1]):
  1. `channel` "The channel" (row 0, flex 1.5) — side view of the bed, the upstream depth H₁ and downstream H₂ joined by a
     turbulent roller (drawn with a few curls), the dashed CV with faces 1 and 2, momentum-flux arrows ρQu (purple, lengths
     ∝ value) and face pressure forces ½ρgH² (amber) at each face, tracer dots moving at u₁ and u₂ (transport); in bore mode
     the step moves and the water ahead is still. In the solitary-wave fate: the hump moving at c₀(1 + a/2H).
  2. `curve` "Depth ratio and loss vs Fr₁" (row 1, flex 1.1) — H₂/H₁ (black) and head loss/H₁ (rose) for Fr₁ ∈ [0.3, 6]; the
     Fr₁ < 1 region greyed "would create energy"; the current dot. Pointer: click → set Fr₁.
  3. `budget` "Momentum budget" (row 1, flex 1, `hidePortrait: true`) — the `terms` bars: inflow ρQu₁, pressure ½ρgH₁²,
     outflow −ρQu₂, pressure −½ρgH₂², total ≈ 0 (black); and an energy bar ΔE (rose if negative, red "forbidden" if positive).
  Portrait: `channel` + `curve`; H₂ and the loss in the status.
- **Controls:** `H1` "Upstream depth $H_1$" 0.02 … 2 m, default 0.1 · `Fr1` "Froude number $\mathrm{Fr}_1$" 0.5 … 6, step
  0.01, default 3 · `frame` chips · `fate` chips · `aH` "Hump $a/H$" 0.02 … 0.6 (optional, solitary fate).
- **Transport:** `t` 0 → 10 s, rate 1, loop (tracers; the bore moving).
- **Presets:** "undular jump" {Fr1: 1.3} · "spillway" {Fr1: 5, H1: 0.1} · "tidal bore" {frame: 'bore into still water',
  H1: 1, Fr1: 1.5} · "Fr₁ = 1 (no jump)" {Fr1: 1} · "Fr₁ = 0.7 (forbidden)" {Fr1: 0.7} · "solitary wave a/H = 0.2" {fate:
  'solitary wave', aH: 0.2}.
- **Status:** "✅ jump: Fr₁ = 3.00 > 1 — H₂ = 0.377 m, 0.141 m of head lost" · "〰 undular: 1 < Fr₁ < 1.7 — a train of
  waves instead of a roller" · "⛔ Fr₁ < 1: this 'jump' would create energy — forbidden" · solitary: "🌊 solitary wave:
  c = c₀(1 + a/2H) = 3.45 m/s · Ursell a λ²/H³ = …".
- **Readouts:** "$H_2$" (m) · "$u_2$" (m/s) · "$\mathrm{Fr}_2$" · "Head loss" (m) · "Power lost" (W/m).
- **Terms:** title "Momentum through the CV (per m of width)" unit N/m: in ρQu₁ (purple) · p₁ = ½ρgH₁² (amber) · out −ρQu₂
  (purple, hatched) · −p₂ = −½ρgH₂² (amber, hatched) · total (black, ≈ 0).
- **Depth features:** Explain + Code + Derivation · **linked views** (3) · **term bars** · **presets** · **status** ·
  **modes** (frames; fates) · **transport**.
- **Explain:**
  0. *What the views show* — "**The channel**: the dashed box is the control volume; purple arrows carry momentum in and
     out, amber arrows are the hydrostatic pushes on the two faces. **Depth ratio and loss vs Fr₁**: the grey region is
     forbidden. **Momentum budget** (hidden on phones): the four bars add to zero."
  1. *Upstream state* — "u₁ = Fr₁√(gH₁) = **Fr1** × √(9.81 × **H1**) = **u1** m/s; Q = u₁H₁ = **Q** m²/s (R10: (4.104) Fr =
     u/√(gH))."
  2. *The depth ratio* — "(7.81) H₂/H₁ = ½(−1 + √(1 + 8Fr₁²)) = ½(−1 + √(1 + 8 × **Fr1sq**)) = **r**" boxed; "H₂ = **H2** m".
  3. *Downstream state* — "u₂ = Q/H₂ = **u2** m/s; Fr₂ = u₂/√(gH₂) = **Fr2** (always < 1 after a jump)."
  4. *Momentum check* — "(7.80): ρQ(u₂ − u₁) = **lhs** N/m and ½ρg(H₁² − H₂²) = **rhs** N/m — equal (residual **res**)."
  5. *Energy* — "E₂ − E₁ = −g(H₂ − H₁)³/(4H₁H₂) = **dE** J/kg; head loss **hl** m; power turned into heat ρgQ × head loss =
     **P** W per metre of width."
  6. *Bore frame (if chosen)* — "moving into still water of depth H₁ the step travels at **cb** m/s = u₁ of the stationary
     picture."
  7. *Solitary fate (if chosen)* — "(7.88) c = c₀(1 + a/2H) = √(gH)(1 + **aH**/2) = **cs** m/s; Ursell a λ²/H³ with λ ≈ the
     hump width ×4 = **U** (< 16: dispersion can balance steepening)."
  8. *Reading the current setting* — Fr₁ > 1.7: "A proper jump: the fast shallow flow cannot send waves upstream (Fr₁ > 1),
     so it piles up abruptly; momentum fixes how much, and the roller burns **pct** % of the upstream energy." · 1 < Fr₁ <
     1.7: "A weak (undular) jump: the loss is tiny and the step breaks into a train of waves instead of a roller." · Fr₁ <
     1: "The algebra gives an answer but the energy would rise: the second law rules it out — subcritical flow does not
     jump." · solitary: "Steepening (the ηη_x term of (7.87)) and dispersion (the η_xxx term) balance exactly: the hump
     travels unchanged; taller humps are faster."
- **Derivation tab:** **D25** (10 steps) `view: 'channel'`, goal `set` spillway; step 3 (the pressure forces) highlight the
  amber arrows `watch` "½ρgH² on each face"; step 4 highlight purple arrows; step 9 **live** "r = ½(−1 + √(1 + 8 × **Fr1²**))
  = **r**". **D26** (8 steps) `view: 'curve'` (phones `channel`), step 7 (the sign) `set {Fr1: 0.7}` with `watch` "the energy
  bar turns positive: forbidden"; step 8 **live** "E₂ − E₁ = **dE** J/kg"; interpret "At Fr₁ = **Fr1** the jump costs
  **hl** m of head."
- **Code:**
  ```python
  H1, Fr1 = {{H1}}, {{Fr1}}                   # upstream depth [m], Froude number
  j = ch07.jump_state(H1, Fr1)                # (7.80)-(7.81)
  print(j["H2"], j["u2"], j["Fr2"])           # {{H2}} m, {{u2}} m/s, {{Fr2}}
  print(j["mom_in"] + j["p_in"],              # momentum in + push in
        j["mom_out"] + j["p_out"])            # = out + push out: {{bal}} N/m
  print(j["dE"], j["head_loss"], j["allowed"])   # {{dE}} J/kg, {{hl}} m, {{ok}}
  print(ch07.ursell_number({{a}}, {{lam}}, H1))  # aλ²/H³ = {{U}}
  ```
- **Walkthrough (6 steps):** 1. "The roller" — "At a spillway's foot a thin fast sheet rises into a deep slow roller. How
  high?" `set` spillway, `play: true` · 2. "The box" — "Draw a box round the jump: momentum in, momentum out, and the two
  hydrostatic pushes." `terms: true`, `derive: {id: 'D25', step: 3}` · 3. "Momentum decides" — "Balancing them gives (7.81):
  H₂/H₁ = ½(−1 + √(1 + 8Fr₁²)) = 3.77 here." `derive: {id: 'D25', step: 9}`, `code: {lines: [2, 3]}` · 4. "Energy is lost" —
  "The head drops by (H₂ − H₁)³/(4H₁H₂) = 14 cm: burnt in the roller." `readouts: ['loss']` · 5. "Forbidden jumps" — "Drag
  Fr₁ below 1: the formula still answers, but energy would appear from nowhere." `set {Fr1: 0.7}`, `derive: {id: 'D26', step:
  7}` · 6. "Your turn" — "Predict H₂ for Fr₁ = 2 on H₁ = 0.5 m, then check." `controls: ['H1', 'Fr1']`.
- **Equations:** `Fr` ref 'Eq. (4.104)' `\mathrm{Fr}=u/\sqrt{gH}` · `mom` ref 'Eq. (7.80)' · `bel` ref 'Eq. (7.81)' · `loss`
  (unnumbered) `E_2-E_1=-g(H_2-H_1)^3/(4H_1H_2)` · `kdv` ref 'Eq. (7.87)' · `sol` ref 'Eq. (7.88)'.
- **Check yourself:** (1) "Double H₁ at fixed Fr₁. What happens to H₂/H₁ and to the head loss?" — "The ratio stays (it
  depends on Fr₁ only); the loss doubles (it scales with H₁)." · (2) "Why is Fr₂ always below 1?" — "Downstream the water
  is deeper and slower; a jump takes supercritical flow to subcritical, like a shock takes supersonic to subsonic (Ch.
  15)." · (3) "What happens at Fr₁ = 1 exactly?" — "H₂ = H₁ and no loss: no jump." `set {Fr1: 1}` · (4) "A bore 1 m high
  enters a river 1 m deep. How fast does it run?" — "H₂/H₁ = 2 ⇒ 2Fr₁² = 6, Fr₁ = √3: c = √3 × √(9.81) = 5.4 m/s."
- **Selftest parity rows:** `{name: 'H2 spillway', js: jump(0.1, 3).H2, py: 'ch07.jump_state(0.1, 3.0)["H2"]', rtol: 1e-12}` ·
  `{name: 'head loss', js: jump(0.1, 3).head_loss, py: 'ch07.jump_state(0.1, 3.0)["head_loss"]', rtol: 1e-12}` · `{name:
  'Fr2', js: jump(0.5, 2).Fr2, py: 'ch07.jump_state(0.5, 2.0)["Fr2"]', rtol: 1e-12}` · `{name: 'bore speed', js: jump(1,
  Math.sqrt(3)).bore_speed, py: 'ch07.jump_state(1.0, 1.7320508075688772, frame="moving")["bore_speed"]', rtol: 1e-12}` ·
  `{name: 'ursell', js: ursell(0.2, 10, 1), py: 'ch07.ursell_number(0.2, 10.0, 1.0)', rtol: 1e-12}` · `{name: 'soliton eta',
  js: sol(1.5, 0.0, 0.2, 1.0), py: 'ch07.solitary_wave(1.5, 0.0, 0.2, 1.0)', rtol: 1e-12}` · invariant `{name: 'momentum
  balance', js: jump(0.1, 3).residual, expect: 0, atol: 1e-9}`.
- **Fit plan:** 360×640: `channel` (55 %) + `curve`; budget hidden (its total in Explain §4), chips wrap. 844×390:
  `channel` | `curve`. Desktop: rows [1.2, 1], `budget` beside `curve`.

### E8 · two_layer_modes
- **Title:** "What is a baroclinic mode?" · **Summary:** "A layer over deep water rings in two ways: the barotropic mode
  (surface and interface in phase, fast) and the baroclinic mode (antiphase, slow, √(g′H) for long waves, almost
  invisible at the surface); with no free surface only the interfacial wave ε√(gk) remains." · **CORE:** C13, C14 (also
  N99 complex amplitudes, N105–N107 ((7.96), vortex sheet), N124–N131 ((7.111)–(7.119)), N128 g′ convention; links ch05
  `vortex_sheet_rollup`) · **Reference:** `amplitude_phase_second_order_II_3.html` (one state driving system and response
  windows; a numbered derivation with live numbers).
- **meta:** `viz:order 8` · `viz:sections 7.7` · `viz:equations 7.95 7.96 7.109 7.110 7.111 7.112 7.113 7.114 7.116 7.117`
  · `viz:fluidpy ch07.interface_omega ch07.two_layer_free_surface_omega ch07.two_layer_state ch07.two_layer_modes
  ch07.reduced_gravity_book` · `viz:derivations D28 D29 D30 D31`.
- **Physics:** `wInt(k)` ↔ `ch07.interface_omega(k, rho1, rho2)`; `wTwo(k, H)` → [ω_bt, ω_bc] ↔
  `ch07.two_layer_free_surface_omega(k, H, rho1, rho2)`; `ratio(k, H, mode)` η/ζ ((7.112): e^{kH} for barotropic;
  (7.114)) ↔ `ch07.two_layer_modes(k, H, rho1, rho2, mode=)["eta_over_zeta"]`; `gprime(ref)` ↔
  `ch07.reduced_gravity_book(rho1, rho2, ref=)`; `state` ↔ `ch07.two_layer_state`. ρ₂ = 1025 kg/m³ fixed; ρ₁ = ρ₂(1 −
  Δρ/ρ₂). Real fields: η = a cos θ, ζ = |b| cos(θ + arg b) (b real here: sign only), u₁ from φ₁ with A, B of (7.106)–(7.107),
  u₂ from C (7.108).
- **Modes (`mode` chips):** "barotropic" · "baroclinic" · "two deep fluids (no surface)".
- **Views** (rows [1.25, 1]):
  1. `side` "Side view" (row 0, flex 1.5) — one wavelength (or 2) with the free surface (blue) and the interface (navy line
     over a navy fill for the lower layer), both animated with exaggeration factors printed ("surface ×**fs**, interface
     ×**fi**"); arrows of u in both layers at 5 depths each (upper light, lower navy), reversing across the interface in the
     baroclinic mode; the sheet jump marked (± signs) at the interface.
  2. `disp` "ω(k) of both modes" (row 1, flex 1.1) — log–log k ∈ [10⁻⁵, 1] rad/m: barotropic (blue), baroclinic (navy), the
     limits √(gk) (dashed), ε√(gk) (7.95) (navy dashed), √(g′H)k (dotted); the current dot on the chosen mode.
  3. `ratio` "Surface signal and speed vs Δρ/ρ" (row 1, flex 1, `hidePortrait: true`) — |η/ζ| of the baroclinic mode and
     c_bc/c_bt against Δρ/ρ on log–log axes, the current value marked.
  Portrait: `side` + `disp`; η/ζ and both periods in the status.
- **Controls:** `mode` chips · `drho` "Density step $\Delta\rho/\rho_2$" log 10⁻⁴ … 0.5, default 0.002 · `H` "Upper layer
  $H$" log 1 … 500 m, default 50 m · `lam` "Wavelength $\lambda$" log 10 m … 100 km, default 1000 m · `a` optional.
- **Transport:** `tau` 0 → 3 periods of the chosen mode, rate 0.5, loop.
- **Presets:** "ocean thermocline" {H: 50, drho: 0.002, lam: 1000} · "fjord (fresh over salt)" {H: 5, drho: 0.02, lam:
  200} · "oil over water" {drho: 0.2, H: 1, lam: 20} · "kH → ∞ (recover (7.95))" {H: 500, lam: 100} · "long waves (g′H)"
  {lam: 50000}.
- **Status:** "🌊 barotropic: interface in phase, e^(−kH) = 0.73 of the surface · T = 25 s" · "↕ baroclinic: antiphase,
  surface signal ×0.0015 · T = 19.5 min · c = 0.85 m/s" · "⇅ two deep fluids: ω = ε√(gk), ε = 0.032".
- **Readouts:** "Period (barotropic)" · "Period (baroclinic)" · "$\eta/\zeta$" · "$g'$ (book, ρ₂)" · "$\sqrt{g'H}$".
- **Depth features:** Explain + Code + Derivation · **linked views** (3) · **transport** · **modes** · **presets** ·
  **status** · **notes** (the "which ρ?" table: g′ with ρ₂ and with ρ₁ and their % difference at the current Δρ).
- **Explain:**
  0. *What the views show* — "**Side view**: blue = free surface, navy = interface (drawn with separate magnifications, both
     printed), arrows = horizontal velocity in each layer. **ω(k)**: blue = barotropic, navy = baroclinic; dashed = the
     deep-water and two-deep-fluid limits, dotted = the long-wave limit. **Surface signal** (hidden on phones)."
  1. *Density numbers* — "ρ₂ = 1025, ρ₁ = **r1** kg/m³; ε² = (ρ₂ − ρ₁)/(ρ₂ + ρ₁) = **eps2** (7.95); g′ = g(ρ₂ − ρ₁)/ρ₂ =
     **gp** m/s² (7.117) (with ρ₁ below the line: **gp1**, **diff** % apart)" boxed.
  2. *Barotropic mode* — "(7.111) ω² = gk: ω = **wbt** rad/s, T = **Tbt**; (7.112) b = a e^(−kH) = **bt** a (in phase)."
  3. *Baroclinic mode* — "(7.113) ω² = gk(ρ₂ − ρ₁) sinh kH/(ρ₂ cosh kH + ρ₁ sinh kH) = **wbc**² → T = **Tbc**, c = **cbc** m/s"
     boxed; "(7.114) η/ζ = −((ρ₂ − ρ₁)/ρ₁) e^(−kH) = **etaz**".
  4. *Long-wave limit* — "kH = **kH**; (7.116) c = √(g′H) = √(**gp** × **H**) = **cl** m/s vs the barotropic √(gH) = **cs**
     m/s (ratio **q**)."
  5. *Two deep fluids* — "(7.95) ω = ε√(gk) = **wint** rad/s — the baroclinic branch for kH ≫ 1."
  6. *A 10 m internal wave* — "surface signal |η| = 10 × |η/ζ| = **sig10** m; energy (7.96) ½(ρ₂ − ρ₁)ga² with a = 10 m =
     **Eint** J/m² — equal to a surface wave of height **aeq** m."
  7. *Right now* — "t = **live t**; surface η = **live eta** m, interface ζ = **live zeta** m."
  8. *Reading the current setting* — barotropic: "The whole column moves together (in the Boussinesq limit even u is the
     same at every depth): the surface carries the signal — tides and storm surges are barotropic." · baroclinic: "The
     interface heaves **ratio** times more than the surface, in antiphase, and the layers slide past each other: a vortex
     sheet at the interface. Its speed √(g′H) ≈ **cl** m/s is ~**q** times slower than the surface wave — the
     thermocline's own slow dynamics (Ch. 13's reduced-gravity models)." · two deep fluids: "No surface: only the
     interfacial wave, a deep-water wave slowed by ε."
- **Derivation tab:** **D28** (10 steps) `view: 'side'`, goal `set {mode: 'two deep fluids'}`; step 5 **live** "A = iωa/k =
  i × **wa/k**"; step 9 **live** "ω² = gk(ρ₂ − ρ₁)/(ρ₂ + ρ₁) = …"; step 10 `set {drho: 0.5}` with `watch` "ρ₁/ρ₂ → 0 would give
  √(gk)". **D29** (12 steps) `view: 'side'`, goal `set` ocean thermocline; step 7 **live** "A = −(ia/2)(ω/k + g/ω) = …"; step
  12 **live** "b/a = **ba**". **D30** (12 steps) `view: 'disp'` (phones `side`), step 10 (factor out s − 1) highlights the blue
  branch, step 12 highlights the navy branch. **D31** (11 steps) `view: 'disp'`, step 3 `set {mode: 'barotropic'}`, step 5 `set
  {mode: 'baroclinic'}`, step 9 `set` long waves with **live** "c = √(g′H) = **cl** m/s".
- **Code:**
  ```python
  rho2 = 1025.0; rho1 = rho2*(1 - {{drho}})   # densities [kg/m^3]
  k, H = 2*np.pi/{{lam}}, {{H}}               # kH = {{kH}}
  w_bt, w_bc = ch07.two_layer_free_surface_omega(k, H, rho1, rho2)
  print(2*np.pi/w_bt, 2*np.pi/w_bc)           # (7.111), (7.113): {{Tbt}}, {{Tbc}} s
  m = ch07.two_layer_modes(k, H, rho1, rho2, mode="{{mode}}")
  print(m["eta_over_zeta"])                   # (7.112)/(7.114): {{etaz}}
  print(ch07.reduced_gravity_book(rho1, rho2),       # g' = {{gp}} (ρ₂)
        ch07.two_layer_long_wave_speed(H, rho1, rho2))  # √(g'H) = {{cl}} m/s
  ```
- **Walkthrough (7 steps):** 1. "Two ways to ring" — "A warm layer floats on cold water. Watch the surface and the
  thermocline: how can this system oscillate?" `set` ocean thermocline, `play: true` · 2. "One interface" — "Without the
  free surface: ω = ε√(gk) (7.95), a deep-water wave slowed by ε = 0.03." `set {mode: 'two deep fluids'}`, `derive: {id:
  'D28', step: 9}` · 3. "Add a surface" — "Now the pressure condition factors into two roots (7.110): two modes." `derive:
  {id: 'D30', step: 10}` · 4. "Barotropic" — "ω² = gk: surface and interface in phase, the interface smaller by e^(−kH)
  (7.112)." `set {mode: 'barotropic'}` · 5. "Baroclinic" — "The slow root: interface metres, surface centimetres, in
  antiphase (7.114)." `set {mode: 'baroclinic'}`, `derive: {id: 'D31', step: 5}` · 6. "Reduced gravity" — "Long waves: c =
  √(g′H) ≈ 1 m/s with g′ = g(ρ₂ − ρ₁)/ρ₂ (7.117) — 22 times slower than √(gH)." `set` long waves, `notes: true` · 7. "Your
  turn" — "Predict how the baroclinic period changes if Δρ doubles, then drag." `controls: ['drho', 'H']`.
- **Equations:** `int` ref 'Eq. (7.95)' · `E` ref 'Eq. (7.96)' · `b` ref 'Eq. (7.109)' · `disp2` ref 'Eq. (7.110)' · `bt` ref
  'Eq. (7.111)' `\omega^2=gk` · `bc` ref 'Eq. (7.113)' · `ratio` ref 'Eq. (7.114)' · `gp` ref 'Eq. (7.116, 7.117)'.
- **Check yourself:** (1) "Double Δρ in the long-wave preset. What happens to the baroclinic speed?" — "×√2: c = √(g′H),
  g′ ∝ Δρ." · (2) "Why is the barotropic period independent of Δρ?" — "Its root ω² = gk (7.111) contains no density: the
  surface wave hardly notices the weak step." · (3) "Make kH large (H = 500 m, λ = 100 m). Which formula does the
  baroclinic curve follow?" — "(7.95), ε√(gk): the interface no longer feels the far surface." · (4) "Oil over water: how
  different are the two g′ conventions?" — "g(ρ₂ − ρ₁)/ρ₂ = 1.96 vs /ρ₁ = 2.45 m/s² for 800 over 1000: 25 %." `set` oil over
  water.
- **Selftest parity rows:** `{name: 'omega interface', js: wInt(2*Math.PI/100), py: 'ch07.interface_omega(2*np.pi/100.0,
  1022.95, 1025.0)', rtol: 1e-12}` · `{name: 'omega baroclinic', js: wTwo(2*Math.PI/1000, 50)[1], py:
  'ch07.two_layer_free_surface_omega(2*np.pi/1000.0, 50.0, 1022.95, 1025.0)[1]', rtol: 1e-12}` · `{name: 'eta/zeta bc', js:
  ratio(2*Math.PI/1000, 50, 'baroclinic'), py: 'ch07.two_layer_modes(2*np.pi/1000.0, 50.0, 1022.95, 1025.0,
  mode="baroclinic")["eta_over_zeta"]', rtol: 1e-10}` · `{name: 'g prime book', js: gprime('lower'), py:
  'ch07.reduced_gravity_book(1022.95, 1025.0)', rtol: 1e-12}` · `{name: 'c long', js: cLong(50), py:
  'ch07.two_layer_long_wave_speed(50.0, 1022.95, 1025.0)', rtol: 1e-12}` · invariant `{name: 'kH→∞ limit', js: wTwo(1,
  1000)[1] - wInt(1), expect: 0, atol: 1e-9}`. (ρ₁ = 1025 × (1 − 0.002) = 1022.95 at the default.)
- **Fit plan:** 360×640: `side` (55 %) + `disp`; mode chips one row; `ratio` hidden (its numbers in the status). 844×390:
  `side` | `disp`. Desktop: rows [1.25, 1], `ratio` beside `disp`.

### E9 · internal_wave_beams
- **Title:** "Why does energy leave at right angles to the crests?" · **Summary:** "In a stratified fluid the frequency sets
  only the angle, ω = N cos θ: drag ω/N and the St Andrew's cross tilts, the crests slide across each beam while the
  energy runs along it, and in the wavenumber plane c_g stays perpendicular to c." · **CORE:** C15, C16 (also N141–N145
  ((7.137)–(7.141)), N146 layered flow, N148–N150 ((7.143)–(7.145)), N158–N166 ((7.153)–(7.159)); links ch01
  `parcel_stability`) · **Reference:** `angular_frequency_explorer_1.html` (rotating pointer + system + response curve on
  one clock) — the K arrow is its rotating pointer.
- **meta:** `viz:order 9` · `viz:sections 7.8` · `viz:equations 7.134 7.137 7.139 7.141 7.144 7.145 7.146 7.159` ·
  `viz:fluidpy ch07.internal_wave_state ch07.internal_wave_omega ch07.internal_wave_velocities ch07.st_andrews_cross
  ch07.internal_wave_energy ch07.beam_angle` · `viz:derivations D33 D34 D35 D36`.
- **Physics:** `wInt(k, m, N)` = N|k|/√(k² + m²) ↔ `ch07.internal_wave_omega(k, m, N)`; `vel(k, m, N)` = {c = (ω/K²)(k, m),
  cg = ∇_K of N|k|/K} ↔ `ch07.internal_wave_velocities(k, m, N)`; `state(ωN, N, K, ks, ms)` ↔ `ch07.internal_wave_state(
  omega_over_N, N, K, k_sign, m_sign)`; `cross(x, z, t)` (four Gaussian beams) ↔ `ch07.st_andrews_cross(x, z, t, omega, N,
  width)`; `energy(k, m, N, w0)` ↔ `ch07.internal_wave_energy`.
- **Modes (`view` chips):** "beams (St Andrew's cross)" · "single plane wave" · "packet".
- **Views** (rows [1.3, 1]):
  1. `tank` "The stratified tank" (row 0, flex 1.5, equal aspect) — faint horizontal density layers (background), the ρ′
     field of the chosen mode as a blue–red heatmap (`Viz.field.heatmap`): four beams from a small oscillating cylinder at
     the centre (beams mode), one plane wave (plane mode) or a Gaussian packet moving along c_g (packet mode); moving phase
     lines; 12 particle dots oscillating along the crests (teal double arrows); the beam angles labelled to the vertical and
     to the horizontal.
  2. `kplane` "Wavenumber plane" (row 0, flex 1) — axes k, m (rad/m); the K arrow (black, draggable tip), c (orange) along K,
     c_g (purple) at 90° with a right-angle mark; the circle |K| = const (grey) and rays of constant ω (grey, labelled
     ω/N = 0.25, 0.5, 0.75). Pointer: drag the tip → sets `theta` and `K` (and the signs of k, m by quadrant).
  3. `curve` "ω/N vs θ" (row 1, `hidePortrait: true`) — the cos curve for θ ∈ [0°, 90°], the current dot, and the parcel
     limit ω = N at θ = 0 marked "Ch. 1 parcel".
  Portrait: `tank` + `kplane`; the angles in the status.
- **Controls:** `wN` "Frequency $\omega/N$" 0.05 … 0.99, default 0.71 · `N` "Buoyancy frequency $N$" log 0.001 … 2 rad/s,
  default 1 · `view` chips · `K` "Wavenumber $\lvert\mathbf K\rvert$" 0.2 … 5 rad/m (optional; shows that ω does not change) ·
  `quad` chips "K up-right · up-left (Fig. 7.29)" (optional).
- **Transport:** `tau` 0 → 3 periods of ω, rate 0.5, loop.
- **Presets:** "45° cross" {wN: 0.71, N: 1} · "near N" {wN: 0.95} · "low frequency (flat beams)" {wN: 0.15} · "Fig. 7.29
  geometry (k < 0)" {view: 'single plane wave', quad: 'up-left', wN: 0.5} · "ocean thermocline" {N: 0.01, wN: 0.5}.
- **Status:** "✳ beams 44.8° from the vertical (45.2° from the horizontal) · phase ↑ ⇒ energy ↓" · near N: "↕ θ = 18° from
  the vertical: columns bob, energy barely moves" · low: "〰 beams 81° from the vertical: nearly flat, layered flow".
- **Readouts:** "$\theta$ (K from horizontal)" (°) · "Beam from vertical" (°) · "$\lvert\mathbf c\rvert$" (m/s) ·
  "$\lvert\mathbf c_g\rvert$" (m/s) · "$\mathbf c\cdot\mathbf c_g$".
- **Depth features:** Explain + Code + Derivation · **linked views** (3) · **transport** · **modes** · **presets** ·
  **status** · **inspector** (click the kplane: "ω = N|k|/K = 1 × 0.707/1 = **0.707** rad/s; c = (ω/K²)(k, m) = (0.5, 0.5)
  m/s; c_g = (Nm/K³)(m, −k) = (0.5, −0.5) m/s; c · c_g = 0.25 − 0.25 = **0**").
- **Explain:**
  0. *What the views show* — "**The tank**: blue/red = density perturbation ρ′ (crests and troughs), teal arrows = the
     water's motion along the crests, the beams carry energy away from the source. **Wavenumber plane**: black K, orange
     c (along K), purple c_g (at a right angle); grey rays are lines of constant ω. **ω/N vs θ** (hidden on phones)."
  1. *The angle from the frequency* — "(7.139) ω = N cos θ ⇒ θ = arccos(ω/N) = arccos **wN** = **th**° — the angle of K above
     the horizontal = the beam's angle from the vertical" boxed; "from the horizontal the beam makes **90−th**°."
  2. *The wavenumber components* — "k = K cos θ = **k**, m = K sin θ = **m** rad/m (signs from the quadrant)."
  3. *Phase velocity* — "(7.144) c = (ω/K²)(k, m) = (**cx**, **cz**) m/s, |c| = ω/K = **c** m/s."
  4. *Group velocity* — "(7.145) c_g = (Nm/K³)(m, −k) = (**gx**, **gz**) m/s, |c_g| = N sin θ/K = **cg** m/s" boxed.
  5. *At right angles* — "(7.146) c · c_g = **cx**×**gx** + **cz**×**gz** = **dot** — zero: the energy runs along the
     crests; vertical parts opposite: phase **up/down** ⇔ energy **down/up**."
  6. *Energy flux* — "(7.157) E = ½ρ₀(m²/k² + 1)ŵ² = **E** J/m³ (ŵ = 1 cm/s); (7.159) F = c_g E = (**Fx**, **Fz**) W/m²."
  7. *Periods* — "T = 2π/ω = **T** s; buoyancy period 2π/N = **TN** s (the fastest possible)."
  8. *Right now* — "t = **live t**."
  9. *Reading the current setting* — "ω/N = **wN**: " + near N: "K nearly horizontal: the water moves nearly vertically,
     like Ch. 1's bobbing parcel, and c_g → 0 — the energy hardly leaves." · middle: "The classic cross: phase lines move
     across the beams, energy moves along them; if the phase goes up, the energy goes down." · low: "K nearly vertical: the
     motion is almost horizontal and the beams almost flat — the layered, pancake-like flow of strongly stratified fluids
     (N146)."
- **Derivation tab:** **D33** (12 steps) `view: 'tank'`, goal `set {view: 'single plane wave'}`; steps 2, 6, 9 highlight the
  operator applied (a chip "∂/∂t", "∇_H²") in the step text; step 12 **live** "∂²/∂t²∇²w + N²∇_H²w with your wave: residual
  **res**". **D34** (8 steps) `view: 'kplane'`, step 4 **live** "ω²K² = N²k² ⇒ ω = N|k|/K = **w**"; step 7 `set {wN: 0.3}` then
  `{wN: 0.9}` with `watch` "K turns; its length does not matter". **D35** (5 steps) `view: 'tank'`, step 4 `watch` "the teal
  arrows lie along the stripes". **D36** (9 steps) `view: 'kplane'`, step 5 **live** "∂ω/∂k = Nm²/K³ = **gx**", step 6 **live**
  "∂ω/∂m = −Nkm/K³ = **gz**", step 9 `watch` "the right-angle mark"; interpret "At ω/N = **wN** the energy leaves at
  **th**° from the vertical, at **cg** m/s."
- **Code:**
  ```python
  N, wN, K = {{N}}, {{wN}}, {{K}}             # [rad/s], ω/N, |K| [rad/m]
  s = ch07.internal_wave_state(wN, N, K, k_sign={{ks}}, m_sign={{ms}})
  print(s["theta_K_deg"], s["beam_from_vertical_deg"])   # (7.139): {{th}}°
  v = ch07.internal_wave_velocities(s["k"], s["m"], N)
  print(v["c"], v["cg"], v["dot"])            # (7.144)-(7.146): dot = {{dot}}
  e = ch07.internal_wave_energy(s["k"], s["m"], N, 0.01)
  print(e["E"], e["F"], e["cgE"])             # (7.159) F = c_g E
  ```
- **Walkthrough (7 steps):** 1. "A strange cross" — "Wiggle a cylinder in a stratified tank: four beams, not rings. Why?"
  `set` 45° cross, `play: true` · 2. "One frequency, one angle" — "The w-equation (7.134) gives ω = N cos θ (7.139): the
  frequency fixes the direction of K, not its length." `derive: {id: 'D34', step: 4}`, `eq: 'disp'` · 3. "Length does not
  matter" — "Drag K along its ray: ω stays; turn it: ω changes." highlight `view:kplane`, `controls: ['K']` · 4. "Motion along
  the crests" — "Incompressibility gives K · u = 0 (7.141): the water slides along the stripes." `derive: {id: 'D35', step:
  4}` · 5. "Energy at right angles" — "c_g is the gradient of ω in the K plane — perpendicular to K (7.146)." `derive: {id:
  'D36', step: 9}`, `inspect: true` · 6. "Phase up, energy down" — "Watch a beam: stripes move up and out, energy down and
  out." `code: {lines: [4, 5]}` · 7. "Your turn" — "Predict the beam angle for ω = 0.3N, then drag." `controls: ['wN']`.
- **Equations:** `weq` ref 'Eq. (7.134)' `\partial_t^2\nabla^2w+N^2\nabla_H^2w=0` · `disp3` ref 'Eq. (7.137)' · `disp` ref 'Eq.
  (7.139)' `\omega=N\cos\theta` · `trans` ref 'Eq. (7.141)' `\mathbf K\cdot\mathbf u=0` · `ccg` ref 'Eq. (7.144, 7.145)' · `perp`
  ref 'Eq. (7.146)' `\mathbf c_g\cdot\mathbf c=0` · `flux` ref 'Eq. (7.159)' `\mathbf F=\mathbf c_gE`.
- **Check yourself:** (1) "At which ω/N do the beams make 60° with the vertical?" — "cos 60° = 0.5: ω = 0.5N." `set {wN:
  0.5}` · (2) "Halve |K| at fixed direction. What happens to ω, |c| and |c_g|?" — "ω stays; |c| = ω/K and |c_g| = N sin θ/K
  both double." · (3) "Why can't an internal wave have ω > N?" — "cos θ ≤ 1: N is the fastest (the bobbing parcel of Ch. 1)."
  · (4) "In Fig. 7.29's geometry (K up-left), where does the energy go?" — "Down-left: c_g = ∇_K(N|k|/K); the printed
  (7.145) $\mathbf c_g=\frac{Nm}{K^3}(m\mathbf e_x-k\mathbf e_z)$, written for k > 0, would point up-right — wrong." `set` Fig. 7.29 geometry.
- **Selftest parity rows:** `{name: 'omega 45', js: wInt(0.70710678, 0.70710678, 1), py: 'ch07.internal_wave_omega(0.70710678,
  0.70710678, 1.0)', rtol: 1e-12}` · `{name: 'cg x (k<0)', js: vel(-0.70710678, 0.70710678, 1).cg[0], py:
  'ch07.internal_wave_velocities(-0.70710678, 0.70710678, 1.0)["cg"][0]', rtol: 1e-12}` · `{name: 'beam angle', js:
  state(0.71, 1, 1, 1, 1).beam_from_vertical_deg, py: 'ch07.internal_wave_state(0.71, 1.0)["beam_from_vertical_deg"]', rtol:
  1e-12}` · `{name: 'flux Fz', js: energy(0.6, 0.8, 1, 0.01).F[1], py: 'ch07.internal_wave_energy(0.6, 0.8, 1.0,
  0.01)["F"][1]', rtol: 1e-12}` · `{name: 'cross field', js: cross(0.3, 0.2, 1.0), py: 'ch07.st_andrews_cross(0.3, 0.2, 1.0,
  0.71, 1.0, 0.1)', rtol: 1e-9}` · invariant `{name: 'c dot cg', js: dot(vel(0.3, -1.2, 0.5)), expect: 0, atol: 1e-14}`.
- **Fit plan:** 360×640: `tank` (55 %) + `kplane`; `curve` hidden (θ in the status). 844×390: `tank` | `kplane`. Desktop:
  rows [1.3, 1], `curve` under `kplane`.

### B1 · linearised_free_surface
(Backup — built only if one of E1–E9 fails review.)
- **Title:** "How can a condition on an unknown surface be applied at z = 0?" · **Summary:** "Drag the steepness ka: the
  exact kinematic and dynamic conditions, evaluated with the linear solution, miss by an amount that grows like (ka)²,
  while the linearised ones hold exactly; a floating leaf follows ∂η/∂t." · **CORE:** C02 (also N13–N17 ((7.14)–(7.21)),
  R04 (7.13), R05 (7.19), R06 (7.20)) · **Reference:** `overfitting_curves.html` (a minimal two-slider figure with a
  verdict) and the ch04 backup `kinematic_free_surface` storyboard.
- **meta:** `viz:order 10` · `viz:sections 7.2` · `viz:equations 7.13 7.16 7.17 7.18 7.19 7.20 7.21` · `viz:fluidpy
  ch07.free_surface_residuals ch07.free_surface_residual_scan ch07.surface_normal` · `viz:derivations D02 D03 D04`.
- **Physics:** `resid(ka, kH)` ↔ `ch07.free_surface_residuals(x, t, a, k, H)` (central differences on the linear fields);
  `scan(kaList)` ↔ `ch07.free_surface_residual_scan(ka_values, kH)`; `normal(etax)` ↔ `ch07.surface_normal`.
- **Views:** (1) `surface` "The surface and the flat level" — η (blue), the dashed z = 0, a leaf (green) riding the surface
  with its velocity arrow, the normal n (black) and U_s (grey) at the leaf; (2) `bars` "Residuals" — term bars: kinematic
  exact vs linear, dynamic exact vs linear (black, log scale); (3) `loglog` "Residual vs ka" (`hidePortrait`) — slope-2 and
  slope-1 lines with the current ka marked.
- **Controls:** `ka` 0.001 … 0.4 (log), `kH` 0.2 … 5 (log), transport over one period. **Presets:** "swell ka = 0.06",
  "steep ka = 0.3", "shallow kH = 0.3". **Status:** "✅ linear theory good to **p** % (ka = …)" / "⚠️ steep: dropped terms
  **p** % of the kept ones". **Depth features:** linked views (3), term bars, presets, status, transport.
- **Explain:** 0. views · 1. kept term aω = …, dropped slope × velocity ka·aω = … (7.16)→(7.17) · 2. Taylor transfer ηφ_zz =
  … (7.17)→(7.18) · 3. dynamic: ½|∇φ|² vs φ_t, gη at z = η (7.19)–(7.21) · 4. relative error ka = … · 5. reading the current
  setting (swell: invisible; steep: linear theory only a first guess; shallow: a/H is the stricter condition).
- **Derivation tab:** D02 (8), D03 (8), D04 (8) from Part F, `view: 'surface'`; D03 step 5 **live** "ηφ_zz = … vs φ_z = …".
- **Code:** `r = ch07.free_surface_residuals({{x}}, {{t}}, a, k, H)` · `print(r["kinematic_exact"], r["kinematic_linear"])` ·
  `s = ch07.free_surface_residual_scan(np.logspace(-3, -0.4, 12), kH)` · `print(s["slope_abs"], s["slope_rel"])`.
- **Walkthrough (5 steps):** the leaf question → Taylor transfer (derive D03 step 5) → what is dropped (terms) → slope 2 on
  log–log → your turn (find the ka where the error reaches 10 %).
- **Equations:** (7.13), (7.16), (7.18), (7.20), (7.21) written out. **Check yourself (3):** "Halve ka: by what factor does
  the absolute misfit fall?" — "4 (slope 2)" · "Why must both approximations be made together?" — "both dropped terms are
  of relative size ka (D03 step 6)" · "Which condition limits shallow water?" — "a/H ≪ 1".
- **Selftest parity rows:** `{name: 'kin exact', js: resid(0.05, 2).kinematic_exact, py: 'ch07.free_surface_residuals(0.3,
  0.0, 0.05, 1.0, 2.0)["kinematic_exact"]', rtol: 1e-6}` · `{name: 'slope', js: scanSlope(), py:
  'ch07.free_surface_residual_scan([0.001, 0.01, 0.1], 1.0)["slope_abs"]', rtol: 1e-3}` · invariant `{name: 'linear holds',
  js: resid(0.3, 1).kinematic_linear, expect: 0, atol: 1e-8}`.
- **Fit plan:** 360×640: `surface` + `bars`; `loglog` hidden. Desktop: rows [1, 1].

---

## Part D — runtime budget (full run < 5 min on a laptop / Colab CPU)

Chapter 7 is mostly closed forms (µs per point). The costs are: FFT evolutions (`linear_evolve` on 4096 points × ≤ 200
times ≈ 50 ms; `linear_evolve_2d` on 256² × 60 frames ≈ 3 s), the KdV run (IF-RK4 pseudo-spectral, N = 512, ≤ 10⁴
steps ≈ 2 s, cached), exact path lines (DOP853, rtol 1e-10, 20 periods × 6 depths ≈ 1 s), `dblquad` energy checks (3 kH ×
2 ≈ 1 s), ray tracing (≤ 20 rays × 10 ms), five sympy cells (`surface_wave_sympy` ≈ 1 s, `two_layer_sympy` ≈ 3 s,
`boussinesq_linear_sympy` ≈ 2 s, `kdv_residual_sympy` ≈ 1 s, `stokes_expansion_sympy` ≈ 2 s) and the D20 check (≈ 2 s),
six animations and eight plotly figures. The ch06 notebook ran in ≈ 170 s; this one has 16 CORE blocks, 37 derivations
and 9 explainers.

| Section | Heaviest cells | Full | FAST (`FLUIDPY_FAST=1`) |
|---|---|---|---|
| setup + imports | numpy/scipy/sympy/plotly/pint, chapter modules | 9 s | 9 s |
| §7.1 C01 | crest tracking on 120 001 points × 5 times, `linear_evolve` 2 modes, Fig. 7.1 contours (201²) | 3 s | 2 s |
| §7.2 R01–R06, C02 | residual table (2 × 7 stencils), from-scratch residual, residual scan (12 ka × 64 x-points × 7 residuals), two figures | 4 s | 3 s (6 ka) |
| §7.2 C03 | `surface_wave_sympy` (cached), brentq/Newton, V5 table, ω(k) figure | 3 s | 3 s |
| §7.2 C04 | regime table, pressure table, two-panel figure, IF1 slider (30 steps × 3 traces × 300 points), IF3 slider (25 steps), live cell (not run on the page) | 5 s | 3 s (12 steps) |
| §7.2 C05 | `particle_path` linear vs exact (401 points), RK4 by hand (400 steps), four-panel figure, ψ figure (161²), **A1 video 60 frames** (3 panels × 35 particles) | 12 s | 7 s (30 frames) |
| §7.2 C06 | `wave_energy` quad × 3, `energy_flux` quad, midpoint sum 400 × 200, two-panel figure | 4 s | 3 s |
| §7.3 C07 | IAPWS σ, minimum scan 4001 points + `minimize_scalar`, Fig. 7.10 remake, IF2 slider (26 steps) | 3 s | 2 s |
| §7.4 C08 | seiche tables, standing-wave parity, two figures (161² ψ), IF7 slider (6 steps), **A3 video 48 frames** (ψ contours 121² per frame) | 10 s | 6 s (24 frames, 81²) |
| §7.5 C09 | beats, Figs. 7.13–7.15, D20 sympy check (≈ 2 s), from-scratch envelope tracking (`linear_evolve` 4096 × 11 times + Hilbert), **A2 video 100 frames** (packet + pond), IF4 `animate_figure` (30 frames), `pond_ripples` (256 modes × 1000 x × 3 times), Fig. 7.16 | 16 s | 9 s (1024 points, 50 frames, 15 IF4 frames) |
| §7.5 C10 | chirp residuals, Snell ray, `ray_trace` × 9 (+ island 11 rays), RK4 by hand, three-panel figure, IF5 slider (13 steps × 7 rays) | 8 s | 5 s (7 steps) |
| §7.6 C11 | simple-wave frames, jump tables, KdV linear-speed figure, `kdv_residual_sympy`, **A5 frames 30** (+ `kdv_solve` N = 512 cached in `outputs/ch07/kdv_run.npz`: 2 s first run) | 10 s | 6 s (N = 256, 16 frames) |
| §7.6 C12 | `stokes_expansion_sympy` (≈ 2 s), `stokes_drift_numeric` 3 depths × 20 periods, from-scratch 20-period path line, drift figure (6 numeric depths), **A4 video 80 frames** (`dyed_line` 9 particles, exact path lines precomputed once) | 14 s | 8 s (10 periods, 40 frames) |
| §7.7 C13 | interface fields/residuals, `interface_energy` quad, np.linalg.solve, Fig. 7.24 quiver | 3 s | 2 s |
| §7.7 C14 | `two_layer_sympy` ★★★ (D29 and D30 checks share the cached result, ≈ 3 s), np.roots check, three figures, IF6 slider (25 steps × 2 branches × 300 k) | 8 s | 6 s (12 steps) |
| §7.8 R12–R20, C15 | `boussinesq_linear_sympy` ★★★ (≈ 2 s), FD residual of (7.134), two-panel figure, IF8 slider (19 steps) | 5 s | 4 s |
| §7.8 C16 | energy budget residuals, interface limit (4 ε × `quad`), Figs. 7.29/7.31, St Andrew's cross (400², FAST 200²), **A6 video 60 frames** (`linear_evolve_2d` 256², FAST 128²) | 14 s | 7 s |
| explainers (9 × `show_viz`) | read the HTML files | 1 s | 1 s |
| **Total** | | **≈ 132 s** | **≈ 86 s** |

**FAST plan.** Every size-dependent choice is written `a if not FAST else b` in the cell: FFT grids 4096 → 1024 points
(domain length kept, so the physics is unchanged; resolution 2 → 8 points per metre is still ≥ 10 points per carrier
wavelength for λ₀ = 2π m); `linear_evolve_2d` 256² → 128²; KdV N = 512 → 256 (dt from the CFL of the nonlinear term);
Stokes-drift integrations 20 → 10 periods (drift error still ≈ 1 %); videos 60–100 → 30–50 frames and the frames players
30 → 16; plotly sliders ≤ 30 → ≤ 12 steps (≤ 4 traces × ≤ 400 points, < 250 kB each); contour grids 161² → 81² in
animations. **Cached arrays / results:** `ch07.surface_wave_sympy`, `two_layer_sympy`, `boussinesq_linear_sympy`,
`kdv_residual_sympy` and `stokes_expansion_sympy` are `functools.lru_cache`d (D05/D06, D29/D30, D32/D33 read the same
objects); `kdv_solve` writes/reads `outputs/ch07/kdv_run.npz` keyed by its arguments; the C05 path lines and orbit ghosts
are computed once and reused by the figure and A1; the C09 FFT packet run is computed once and reused by the
envelope-tracking check, Fig. 7.17's x–t diagram (C10) and A2; the C12 exact path lines for A4 are precomputed at all
frame times in one `solve_ivp` call per particle (P137 pattern). **Outputs:** 5 videos (A1, A2, A3, A4, A6) + 1 frames
player (A5) at dpi 80 (each < 2.5 MB), 8 plotly figures, ≈ 30 static figures — the page stays well under 15 MB. Sympy
cells never `simplify` expressions with more than two generic functions at once; D33's check substitutes a plane wave
after the operator moves (as `boussinesq_linear_sympy` does) instead of simplifying fourth-order generic expressions.

---

## Part E — prerequisite ledger
Every concept, symbol, maths tool and Python function or idiom the notebook or its explainers use, with where it is
explained. "primer (in Cxx)" = a 📎 primer placed in that block before first use (the primer term is the Concept text —
the builder uses it verbatim in `nb.primer`); "knowledge/primers.md: <term> (chNN Pnn) — reminder" = a one-line
reminder naming the earlier primer; a CORE/RECAP id alone = taught there; "Cxx (Nnnn)" = the NOTE of this chapter placed
in that block; "Cxx (gloss …)" = one sentence where it is used. Earlier chapters' material that is neither a primer nor a
ch07 RECAP is named by section ("Ch. 4 §4.4"). New primers P165–P178 in first-use order: phase of a wave (C01), Taylor
transfer of a boundary condition (C02), separation of variables for a PDE and hyperbolic functions (C03), curvature of a
plane curve and scipy.optimize.minimize_scalar (C07), sum-to-product identities (C08), Fourier integral and a packet's
spectrum and envelope with scipy.signal.hilbert (C09), first-order wave equation and characteristics and Snell's law for
waves (C10), complex amplitudes (C13), operator elimination for linear PDEs (C15), mean of a product of real parts (C16)
— 14 primers.

| Concept | First used in | Explained by |
|---|---|---|
| wave, crest, trough, amplitude a | C01 | C01 |
| wavelength λ and wavenumber k = 2π/λ | C01 | C01 |
| period T, cyclic frequency ν, angular frequency ω | C01 | C01 |
| phase of a wave | C01 | primer (in C01) |
| ω as a frequency (was vorticity in Ch. 3–6) | C01 | C01 (⚠️ callout) |
| sinusoidal travelling wave (7.1), (7.2) | C01 | C01 (N03) |
| crest condition (7.3) | C01 | C01 (N04, D01) |
| phase speed c = ω/k = λν (7.4) | C01 | C01 (D01, N05) |
| solving a linear equation for x | C01 | C01 (D01 step 3) |
| np.argmax and np.polyfit (from-scratch crest tracking) | C01 | C01 (gloss in the from-scratch cell) |
| Fourier superposition of linear waves | C01 | C01 (N02) |
| Fourier modes and the FFT | C01 | knowledge/primers.md: Fourier modes and the FFT Poisson solver (ch05 P142) — reminder |
| plane wave (7.5), wavenumber vector K and K² (7.6) | C01 | C01 (N06, N07) |
| dot product K·x | C01 | C01 (N06; Ch. 2 §2.2) |
| level sets and the normal to a surface | C01 | knowledge/primers.md: level sets and the directional derivative (ch02 P75) — reminder |
| wavelength along K (7.7) and phase velocity vector (7.8) | C01 | C01 (N08, N09) |
| trace velocities | C01 | C01 (N10) |
| Doppler shift (7.9), intrinsic vs observed frequency | C01 | C01 (N11) |
| change of frame x′ = x − Ut | C01 | knowledge/primers.md: frames of reference and relative velocity (ch03 P96) — reminder |
| np.linspace and np.logspace | C01 | knowledge/primers.md: np.linspace and np.logspace (ch01 P06) — reminder |
| matplotlib figures | C01 | knowledge/primers.md: matplotlib figures (ch01 P01) — reminder |
| Python dictionaries (fluidpy results) | C01 | knowledge/primers.md: Python dictionaries (ch01 P23) — reminder |
| f-strings | C01 | knowledge/primers.md: f-strings (ch01 P04) — reminder |
| assert np.allclose | C01 | knowledge/primers.md: assert np.allclose (ch01 P15) — reminder |
| np.meshgrid and the grid layout | C01 | knowledge/primers.md: np.meshgrid and the project grid layout (ch02 P76) — reminder |
| contour and quiver plots | C01 | knowledge/primers.md: plt.contour, plt.quiver and plt.streamplot (ch02 P78) — reminder |
| numpy broadcasting | C01 | knowledge/primers.md: numpy broadcasting (ch02 P77) — reminder |
| ch07.wave_parameters, crest_positions, sinusoid, plane_wave, phase_velocity_vector, trace_velocities, doppler_frequency, linear_evolve | C01 | C01 (code explain) |
| velocity potential φ (7.10) | C02 | R01 |
| Laplace's equation (7.11) | C02 | R02 |
| potential.laplacian_residual | C02 | R02 (code comment) |
| functions as arguments and lambda | C02 | knowledge/primers.md: functions as arguments and lambda (ch01 P29) — reminder |
| no flow through the bottom (7.12) | C02 | R03 |
| kinematic condition at a moving surface (7.13) | C02 | R04 |
| η as elevation vs level-set function (notation change) | C02 | R04, C02 (N13, ⚠️ callout) |
| interfaces.kinematic_bc_residual, ch04.linear_wave_fields | C02 | R04 (code comment) |
| gauge pressure and the stress-free surface (7.19) | C02 | R05 |
| unsteady Bernoulli (4.83) and its linear form (7.20) | C02 | R06 |
| bernoulli.unsteady_bernoulli_pressure, ch07.linear_bernoulli_pressure | C02 | R06 (code explain) |
| small slope and small amplitude, ka ≪ 1, a/H ≪ 1 | C02 | C02 (N12) |
| Taylor transfer of a boundary condition | C02 | primer (in C02) |
| first-order Taylor expansion | C02 | knowledge/primers.md: first-order Taylor expansion (ch01 P26) — reminder |
| unit normal of the surface (7.14) | C02 | C02 (N13, D02) |
| moving level set and its normal speed | C02 | knowledge/primers.md: moving level set and its normal speed (ch04 P132) — reminder |
| chain rule along a path | C02 | knowledge/primers.md: multivariable chain rule along a path (ch03 P91) — reminder |
| surface velocity U_s (7.15) | C02 | C02 (N14) |
| exact kinematic condition (7.16) and D(z − η)/Dt = 0 | C02 | C02 (D02, N15; Ch. 4 §4.10) |
| limits and orders of smallness | C02 | knowledge/primers.md: limits and orders of smallness (ch02 P68) — reminder |
| order-of-magnitude scaling | C02 | knowledge/primers.md: order-of-magnitude scaling (ch04 P130) — reminder |
| linearised kinematic condition (7.17), (7.18) | C02 | C02 (D03, N16) |
| linearised dynamic condition (7.21) | C02 | C02 (D04, N17) |
| the linear free-surface problem (7.11), (7.12), (7.18), (7.21) | C02 | C02 |
| central differences | C02 | knowledge/primers.md: finite differences (ch01 P21) — reminder |
| power laws, log–log slopes, observed order | C02 | knowledge/primers.md: power laws and log–log plots (ch01 P13) — reminder |
| ch07.free_surface_residuals, free_surface_residual_scan, surface_normal | C02 | C02 (code explain) |
| dispersion relation (7.28) | C03 | C03 (D06) |
| separation of variables for a PDE | C03 | primer (in C03) |
| linear second-order ODE and the e^{λz} trial | C03 | knowledge/primers.md: linear second-order ODE (ch01 P44) — reminder |
| hyperbolic functions cosh, sinh, tanh | C03 | primer (in C03) |
| exponent rules | C03 | knowledge/primers.md: exponent rules (ch01 P43) — reminder |
| two linear equations in two unknowns | C03 | C03 (D05 steps 7–8) |
| matching coefficients of cos(kx − ωt) | C03 | C03 (gloss in D06 step 4) |
| potential (7.26) and velocities (7.27) | C03 | C03 (D05, N23, N24) |
| T–λ form of (7.28), coth | C03 | C03 (D06 step 7) |
| sympy | C03 | knowledge/primers.md: sympy (ch01 P40) — reminder |
| scipy.optimize.brentq | C03 | knowledge/primers.md: scipy.optimize.brentq (ch03 P108) — reminder |
| Newton's method | C03 | knowledge/primers.md: Newton's method for complex zeros (ch06 P152) — reminder |
| overflow of cosh/sinh and the stable ratio | C03 | C03 (gloss in the code explain); knowledge/primers.md: np.expm1 and cancellation near zero (ch03 P107) — reminder |
| explicit dispersion approximations (Fenton–McKee, Guo) | C03 | C03 (V5 code explain) |
| ch07.omega_gravity, period_from_wavelength, wavenumber_from_omega, surface_wave_sympy, wave_fields | C03 | C03 (code explain) |
| perturbation pressure p′ = p + ρgz (7.30) | C04 | R07 |
| phase speed (7.29) and dispersion | C04 | C04 (D07) |
| product rule | C04 | knowledge/primers.md: product rule for differentials (ch01 P38) — reminder |
| deep-water limit (7.45) | C04 | C04 (D08, N38) |
| shallow-water limit (7.49), nondispersive | C04 | C04 (D09, N43) |
| wave pressure (7.31) | C04 | C04 (N26) |
| deep-water pressure (7.48), bottom gauge as a low-pass filter | C04 | C04 (N42) |
| hydrostatic pressure of shallow waves (7.52) | C04 | C04 (N46) |
| limits of the hyperbolic functions (Fig. 7.7) | C04 | C04 (N37) |
| Taylor series of tanh and of a square root | C04 | C04 (D09 steps 1–3); knowledge/primers.md: sympy expand, series, removeO, collect and subs (ch04 P117) — reminder |
| bore speed of Ch. 4 Example 4.3 | C04 | C04 (N43; Ch. 4 §4.4) |
| ch07.phase_speed, depth_regime, wavelength_from_period, pressure_response, dispersion_state | C04 | C04 (code explain) |
| slider_figure | C04 | knowledge/primers.md: slider_figure (ch01 P17) — reminder |
| live widgets | C04 | knowledge/primers.md: live widgets (ch01 P47) — reminder |
| show_viz | C04 | knowledge/primers.md: show_viz (ch01 P18) — reminder |
| path lines (7.32) | C05 | R08 |
| ch03.example_3_1 | C05 | R08 (code comment) |
| multivariable Taylor expansion | C05 | knowledge/primers.md: multivariable first-order Taylor expansion (ch03 P98) — reminder |
| linearised path lines (7.34a, b), mean position | C05 | C05 (D10, N28) |
| antiderivatives of sin and cos | C05 | C05 (D10 steps 6–7) |
| particle excursions (7.35a, b) | C05 | C05 (N29) |
| ellipse, semi-axes, foci (7.36) | C05 | C05 (D11); knowledge/primers.md: linear map of a circle is an ellipse (ch03 P104) — reminder |
| sin² + cos² = 1 | C05 | C05 (D11 step 3) |
| clockwise orbits | C05 | C05 (D11 step 6) |
| deep circles (7.46), velocities (7.47) | C05 | C05 (N40, N41) |
| shallow ellipses (7.50), velocities (7.51) | C05 | C05 (N44, N45) |
| stream function of the wave (7.37), u = ∂ψ/∂z, w = −∂ψ/∂x | C05 | C05 (N30; Ch. 4 §4.3) |
| streamfunction.velocity_from_streamfunction_2d | C05 | C05 (N30 code comment) |
| scipy.integrate.solve_ivp | C05 | knowledge/primers.md: scipy.integrate.solve_ivp (ch01 P31) — reminder |
| solve_ivp options and DOP853 tolerances | C05 | C05 (gloss in the code explain); knowledge/primers.md: solve_ivp options: t_eval, dense_output, events, backward integration (ch03 P94) — reminder |
| RK4 by hand | C05 | knowledge/primers.md: RK4 by hand (ch03 P95) — reminder |
| animate and show_animation | C05 | knowledge/primers.md: animate and show_animation (ch01 P16) — reminder |
| ch07.orbit_semi_axes, orbit_linear, particle_path | C05 | C05 (code explain) |
| kinetic energy per unit area (7.38) | C06 | C06 (N31, D12) |
| wavelength average (overbar) vs period average ⟨⟩ | C06 | C06 (gloss in the idea cell) |
| integrals of sin² and cos² over a period | C06 | knowledge/primers.md: integrals of sines and cosines over a full period (ch06 P151) — reminder |
| definite integral | C06 | knowledge/primers.md: definite integral (ch01 P27) — reminder |
| integrals of cosh² and sinh² | C06 | C06 (D12 steps 5–6) |
| potential energy (7.40), column swap | C06 | C06 (N33, D13) |
| equipartition (7.41) | C06 | C06 (N34) |
| wave energy E = ½ρga² (7.42) | C06 | C06 (D13) |
| energy flux as pressure work (7.43) | C06 | C06 (N35, D14) |
| power of a force | C06 | knowledge/primers.md: power of a force and heat flux through a surface (ch04 P127) — reminder |
| flux = energy × group speed (7.44) | C06 | C06 (N36), C09 |
| scipy.integrate.quad and dblquad | C06 | knowledge/primers.md: scipy.integrate.quad and dblquad (ch03 P87) — reminder |
| midpoint double sums | C06 | knowledge/primers.md: volume and surface integrals as midpoint sums (ch02 P83) — reminder |
| energy per unit area vs per unit volume | C06 | C06 (⚠️ callout), C16 |
| ch07.wave_energy, energy_flux, wave_energy_density | C06 | C06 (code explain) |
| Laplace pressure jump (1.5) | C07 | R09 |
| interfaces.laplace_jump_from_balance | C07 | R09 (code comment) |
| capillary waves | C07 | C07 (N49) |
| curvature of a plane curve | C07 | primer (in C07) |
| surface-tension boundary condition (7.53)–(7.55) | C07 | C07 (D15, N50–N52) |
| capillary–gravity dispersion (7.56) and speed (7.57) | C07 | C07 (D15, N53) |
| slope zero at a minimum | C07 | C07 (gloss in D16 step 3); knowledge/primers.md: ordinary derivative as a slope (ch01 P19) — reminder |
| minimum phase speed (7.58), air–water numbers (7.59) | C07 | C07 (D16, N54, N55) |
| pure capillary waves (7.60) | C07 | C07 (N57) |
| scipy.optimize.minimize_scalar | C07 | primer (in C07) |
| ch01.surface_tension_water (IAPWS) | C07 | C07 (code comment; Ch. 1 §1.6) |
| ch07.omega_capillary_gravity, capillary_minimum, capillary_surface_pressure, curvature, capillary_state | C07 | C07 (code explain) |
| left-going wave (7.61) | C08 | C08 (N58) |
| sum-to-product identities | C08 | primer (in C08) |
| standing wave, nodes and antinodes | C08 | C08 (D17) |
| standing-wave stream function (7.62) and velocity (7.63) | C08 | C08 (D17, N59, N61) |
| seiche and the eigenvalue idea | C08 | C08 (N60); knowledge/primers.md: eigenvalues and eigenvectors (ch02 P80) — reminder |
| allowed wavelengths (7.64) and natural frequencies (7.65) | C08 | C08 (D18, N62, N63) |
| rectangular-basin modes | C08 | C08 (N63) |
| ch07.standing_wave_fields, seiche_modes, basin_modes, seiche_state | C08 | C08 (code explain) |
| beats (7.66) and the printed ½Δω x slip | C09 | C09 (D19, N65) |
| group velocity (7.67), chord vs tangent | C09 | C09 (D19) |
| derivative as the limit of a difference quotient | C09 | C09 (gloss in D19 step 7); knowledge/primers.md: ordinary derivative as a slope (ch01 P19) — reminder |
| Fourier integral and a packet's spectrum | C09 | primer (in C09) |
| substitution in an integral | C09 | knowledge/primers.md: substitution in an integral (ch03 P106) — reminder |
| wave packet and envelope (7.68) | C09 | C09 (D20, N66, N67) |
| spreading of a packet (the ω″ term) | C09 | C09 (D20 step 11) |
| group velocity of water waves (7.69), its limits (7.70), c_g = 3c/2 for ripples | C09 | C09 (D21, N68, N69) |
| energy flux F = E c_g (7.71) | C09 | C09 (D21, N70) |
| 3-D group velocity ∂ω/∂K_i | C09 | C09 (N71) |
| complex-step derivative | C09 | C09 (gloss in the code explain) |
| envelope with scipy.signal.hilbert | C09 | primer (in C09) |
| periodic FFT domain and wrap-around | C09 | C09 (gloss in the from-scratch cell); knowledge/primers.md: Fourier modes and the FFT Poisson solver (ch05 P142) — reminder |
| stone in a pond, minimum group velocity | C09 | C09 (N72) |
| viscous decay of waves | C09 | C09 (N72, named; Ch. 8) |
| a crest running through a group | C09 | C09 (N73) |
| animate_figure (plotly time slider) | C09 | C09 (code explain); knowledge/primers.md: slider_figure (ch01 P17) — reminder |
| ch07.group_velocity, group_velocity_numeric, group_velocity_vector, beat_wave, gaussian_packet, packet_state, envelope, pond_ripples, min_group_velocity, viscous_decay | C09 | C09 (code explain) |
| slowly varying wave train (7.72), local k and ω (7.73) | C10 | C10 (N74, N75) |
| mixed partial derivatives commute | C10 | knowledge/primers.md: Schwarz's theorem (ch04 P121) — reminder |
| crest conservation (7.74) | C10 | C10 (D22, N76) |
| first-order wave equation and characteristics | C10 | primer (in C10) |
| k carried at c_g (7.75), x–t diagram | C10 | C10 (D22, N77) |
| inhomogeneous medium ω(k, x) (7.76)–(7.78) | C10 | C10 (N78–N80, D23) |
| frequency constant along rays (7.79) | C10 | C10 (D23) |
| ray equations dx/dt = ∇_k ω, dk/dt = −∇_x ω | C10 | C10 (D23 step 8, D24 step 1) |
| Snell's law for waves | C10 | primer (in C10) |
| refraction on a beach and round an island | C10 | C10 (N47, N48) |
| np.arcsin, np.radians, np.degrees | C10 | C10 (gloss in the code explain); knowledge/primers.md: cosines of angles between unit vectors (ch02 P66) — reminder |
| ch07.local_wavenumber_frequency, crest_conservation_residual, ray_trace, snell_ray_plane_beach, refraction_state | C10 | C10 (code explain) |
| Froude number (4.104), super- and subcritical flow | C11 | R10 |
| similarity.froude_number | C11 | R10 (code comment) |
| nonlinear steepening c = c′ + u | C11 | C11 (N81) |
| simple wave and breaking time (our extension) | C11 | C11 (N81); primer "first-order wave equation and characteristics" (in C10) |
| hydraulic jump, stationary and moving | C11 | C11 (N82) |
| control-volume momentum balance (4.17) | C11 | C11 (D25 step 4; Ch. 4 §4.4); knowledge/primers.md: momentum flux through a surface (ch04 P114) — reminder |
| hydrostatic force on a vertical face ½ρgH² | C11 | C11 (D25 step 3) |
| Bélanger relation (7.80), (7.81) | C11 | C11 (D25, N83) |
| quadratic formula, np.roots, the physical root | C11 | C11 (gloss in D25 step 9); knowledge/primers.md: Vieta's formulas (ch02 P71) — reminder |
| jump energy loss and the second law | C11 | C11 (D26, N84) |
| factorising a cubic difference | C11 | C11 (D26 steps 5–6) |
| ch04.bore_speed | C11 | C11 (code comment; Ch. 4 §4.4) |
| dispersion balancing steepening | C11 | C11 (N85) |
| Korteweg–de Vries equation (7.87) | C11 | C11 (N93) |
| linearised KdV speed | C11 | C11 (N94) |
| Ursell number aλ²/H³ | C11 | C11 (N95) |
| cnoidal waves | C11 | C11 (N96) |
| solitary wave (7.88), sech² and its derivatives | C11 | C11 (N97, gloss); primer "hyperbolic functions cosh, sinh, tanh" (in C03) |
| pseudo-spectral KdV solver and its invariants | C11 | C11 (A5 code explain) |
| ch07.hydraulic_jump, jump_momentum_residual, jump_state, simple_wave_evolve, nonlinear_wavelet_speed, kdv_solve, kdv_invariants, kdv_linear_phase_speed, ursell_number, solitary_wave, cnoidal_wave, kdv_residual_sympy | C11 | C11 (code explain) |
| Stokes wave (7.82), amplitude-dependent speed (7.83), limiting steepness | C12 | C12 (N86, N87, N88) |
| perturbation series checked in sympy | C12 | C12 (N86 code explain); knowledge/primers.md: sympy expand, series, removeO, collect and subs (ch04 P117) — reminder |
| Stokes drift | C12 | C12 (D27) |
| Lagrangian vs Eulerian mean | C12 | C12 (gloss in N89; Ch. 3 §3.2) |
| Taylor-expanded path lines (7.84a, b) | C12 | C12 (N90) |
| Stokes drift at any depth (7.86) | C12 | C12 (N91) |
| zero Eulerian mean current | C12 | C12 (N92) |
| fundamental theorem of calculus | C12 | knowledge/primers.md: fundamental theorem of calculus (ch02 P84) — reminder |
| ch07.stokes_wave_profile, stokes_wave_speed, stokes_expansion_sympy, stokes_drift, stokes_drift_numeric, eulerian_mean_u, dyed_line | C12 | C12 (code explain) |
| Laplace in each fluid (7.90) | C13 | R11 |
| complex amplitudes | C13 | primer (in C13) |
| Euler's formula and i² = −1 | C13 | knowledge/primers.md: square root of a negative number (ch01 P45) — reminder |
| the complex plane in numpy | C13 | knowledge/primers.md: the complex plane in numpy (ch06 P153) — reminder |
| np.linalg.solve with complex numbers | C13 | knowledge/primers.md: np.linalg.solve (ch02 P57) — reminder |
| decaying potentials (7.91), (7.92) | C13 | C13 (D28, N101, N102) |
| two-sided interface conditions (7.93), (7.94) | C13 | C13 (N103, N104) |
| interfacial dispersion ε√(gk) (7.95) | C13 | C13 (D28) |
| interfacial energy (7.96) | C13 | C13 (N105, N106) |
| vortex sheet at the interface | C13 | C13 (N107; Ch. 5 §5.7) |
| ch05.vortex_sheet_strength | C13 | C13 (code comment) |
| heavier fluid on top (Rayleigh–Taylor, named) | C13 | C13 (code explain; Ch. 11) |
| dead water | C13 | C13 (N108) |
| ch07.interface_omega, interface_fields, interface_residuals, interface_energy, real_field | C13 | C13 (code explain) |
| two-layer problem (7.97)–(7.105) and the printed (7.105) slip | C14 | C14 (N110–N118) |
| complex 2 × 2 systems and factorising in sympy | C14 | C14 (gloss in the D29/D30 checks); knowledge/primers.md: sympy for complex series and residues (ch06 P158) — reminder |
| two-layer constants (7.106)–(7.109) | C14 | C14 (D29, N119–N122) |
| two-layer dispersion relation (7.110) | C14 | C14 (D30, N123) |
| barotropic mode (7.111), (7.112) | C14 | C14 (D31, N124, N125) |
| baroclinic mode (7.113), (7.114) | C14 | C14 (D31, N126) |
| long-wave limit (7.115) and reduced gravity (7.116), (7.117) | C14 | C14 (D31, N127, N128) |
| reduced gravity and its two conventions | C14 | C14 (⚠️ callout); knowledge/primers.md: reduced gravity and buoyancy (ch04 P131) — reminder |
| hydrostatic upper layer (7.118), (7.119) | C14 | C14 (N129, N130) |
| two thin layers in the Boussinesq limit, rigid lids | C14 | C14 (N131) |
| ch07.two_layer_free_surface_omega, two_layer_long_wave_speed, reduced_gravity_book, two_layer_modes, two_layer_state, two_layer_sympy, two_layer_rigid_lid_omega | C14 | C14 (code explain) |
| assumptions of §7.8 (Boussinesq, inviscid, linear, no rotation) | C15 | R12 |
| density equation (4.9) and continuity (4.10) | C15 | R13 |
| linear Boussinesq momentum (7.120)–(7.122) | C15 | R14 |
| why Dρ/Dt = 0 | C15 | R15 |
| hydrostatic resting state (7.123) | C15 | R16 |
| background plus perturbation (7.124) | C15 | R17 |
| buoyancy frequency N² (7.127) | C15 | R18 |
| stratification.brunt_vaisala_sq | C15 | R18 (code comment) |
| potential density | C15 | R19 |
| internal waves are rotational | C15 | R20 |
| linearisation about a base state | C15 | C15 (gloss before D32) |
| linearised density equation (7.125), (7.126) | C15 | C15 (D32, N132, N133) |
| perturbation equations (7.128)–(7.131) | C15 | C15 (D32, N134, N135) |
| operator elimination for linear PDEs | C15 | primer (in C15) |
| horizontal Laplacian ∇_H² | C15 | C15 (D33 step 4) |
| w-equation (7.132)–(7.134) | C15 | C15 (D33, N136–N138) |
| anisotropic dispersion ω(K) (7.135) | C15 | C15 (N139) |
| plane-wave trial (7.136), ∂ → ik and −iω | C15 | C15 (N140, D34); primer "complex amplitudes" (in C13) |
| internal-wave dispersion (7.137)–(7.139) | C15 | C15 (D34, N141, N142) |
| angle of a vector, np.arccos | C15 | C15 (gloss in N139) |
| limits ω = N and ω → 0 | C15 | C15 (N145) |
| stratification.parcel_displacement | C15 | C15 (N145 code comment) |
| steady layered flow (7.142) | C15 | C15 (N146) |
| blocking in strong stratification | C15 | C15 (N147, N196) |
| ch07.internal_wave_omega, beam_angle, boussinesq_linear_sympy, w_equation_residual, internal_wave_fields, layered_flow_check | C15 | C15 (code explain) |
| K·u = 0 (7.141), transverse waves | C16 | C16 (D35, N143, N144) |
| group velocity vector (7.143) | C16 | C16 (N148), C09 (N71) |
| gradient in wavenumber space | C16 | C16 (gloss in N148); knowledge/primers.md: partial derivative (ch01 P25) — reminder |
| quotient rule | C16 | C16 (D36 step 3) |
| phase and group velocity of internal waves (7.144), (7.145) | C16 | C16 (D36, N149) |
| c_g ⟂ c (7.146) | C16 | C16 (D36) |
| sign of k in (7.138) and (7.145) | C16 | C16 (N149, D36 step 7) |
| St Andrew's cross and beam angles | C16 | C16 (N150) |
| N(z) and WKB (named) | C16 | C16 (N151) |
| internal-wave energy equation (7.147) | C16 | C16 (N152; Ch. 4 §4.7) |
| available potential energy per volume (7.148)–(7.150) | C16 | C16 (N153, N154, N155) |
| density jump as a delta function (7.151), (7.152) | C16 | C16 (N156, N157); knowledge/primers.md: 2-D divergence theorem and the Dirac delta in the plane (ch06 P149) — reminder |
| mean of a product of real parts | C16 | primer (in C16) |
| complex conjugate | C16 | knowledge/primers.md: complex conjugate (ch02 P81) — reminder |
| polarization relations (7.153) | C16 | C16 (D37, N158, N159) |
| internal-wave energies and equipartition (7.154)–(7.157) | C16 | C16 (N160–N163) |
| energy flux F = c_g E (7.158), (7.159) | C16 | C16 (D37, N164, N165, N166) |
| ch07.internal_wave_velocities, internal_wave_energy, internal_energy_budget_residual, internal_pe_interface_limit, st_andrews_cross, internal_wave_state, linear_evolve_2d | C16 | C16 (code explain) |
| explainer tabs (Walkthrough, Explore, Explain, Derivation, Equations, Code, Check) | C04 | knowledge/primers.md: show_viz (ch01 P18) — reminder |

---

## Part F — derivation storyboards

Builders copy these word for word into `nb.derivation(key, title, goal=…, start=(tex, plain), plan=[…], uses=[…],
steps=[dict(did, tex, why, plain)], result=(tex, plain), interpret=…, check=…, check_src=…)` and into the explainer's
`derivations: [...]` (phones may shorten *why* to its first sentence; `live`, `set` and `watch` are the explainer's and
are listed in Part B). Every step is one move; *why* names the rule and says why we make the move (≤ 35 words, ≥ 6);
*did* ≤ 8 words. The book's own moves were read on the rendered pages (printed p255 for D01, p257–p258 for D02–D04,
p259–p260 for D05–D07, p265–p267 for D08–D09, p260–p262 for D10–D11, p263–p265 for D12–D14, p269–p270 for D15–D16,
p271–p273 for D17–D18, p273–p275 for D19–D21 (D20: the book cites Phillips only), p277–p278 for D22–D23, p267–p268 for
D24 (words only), p282–p283 for D25–D26, p284–p285 for D27, p287–p289 for D28, p290–p292 for D29–D31, p293–p295 for
D32–D33, p296–p297 for D34–D35, p299 for D36, p303–p304 for D37); the gaps listed in `analysis/ch07.md` §2b are filled
and the notebook says so ("the book skips this move"). Every equation named by number is written out. Colours: surface
blue, particles teal, phase orange, group/energy purple, pressure amber, tension rose, gravity blue. (No line in this part
starts with a table bar; absolute values are written with \lvert \rvert or in words.)

### D01 · Crest speed $c=\omega/k=\lambda\nu$ (7.4) from the crest condition (7.3) — ★, 5 steps, in C01 (notebook · `dispersion_relation`)
- **Goal.** Show that every crest of $\eta=a\cos(kx-\omega t)$ (7.2) moves at the same speed ω/k, and write that speed as
  wavelength × frequency.
- **Start.** $\eta(x,t)=a\cos(kx-\omega t)$ (7.2) — *in words:* a cosine surface whose phase changes in space and time.
- **Plan.** (1) Say where a crest is. (2) Turn that into an equation for its position. (3) Follow it for a time Δt.
- **Tools.** phase of a wave (P165) · cos θ = 1 only at θ = 2nπ · solving a linear equation.
- **Assumptions.** A single sinusoid of constant k and ω (step 2).
- **Steps.**
  1. *did:* Say where η is largest · *tex:* $\eta=a\iff\cos(kx-\omega t)=1$ · *why:* A crest is the highest point of the
     surface; since a cosine never exceeds 1, η reaches its maximum a exactly where the cosine equals 1. · *plain:* crests
     sit where the cosine peaks.
  2. *did:* Write the crest's phase · *tex:* $kx_{crest}-\omega t=2n\pi$ (7.3) · *why:* cos θ = 1 exactly when θ is a whole
     multiple of 2π; the integer n labels which crest we follow. · *plain:* each crest keeps its own fixed phase 2nπ.
  3. *did:* Solve for the crest position · *tex:* $x_{crest}=\frac{\omega}{k}t+\frac{2n\pi}{k}$ · *why:* A linear
     equation in x_crest: add ωt to both sides and divide by k, which is not zero for a wave. · *plain:* crest n starts at
     2nπ/k and moves steadily in time.
  4. *did:* Compare two times Δt apart · *tex:* $\Delta x_{crest}=\frac{\omega}{k}\Delta t$ · *why:* Subtract the positions
     at t and t + Δt: the label term 2nπ/k cancels, so every crest moves the same distance. · *plain:* in Δt every crest
     moves (ω/k)Δt. · *live:* "Δx = c × 1 s = … m".
  5. *did:* Divide by Δt and rewrite · *tex:* $c=\frac{\omega}{k}=\frac{2\pi\nu}{2\pi/\lambda}=\lambda\nu$ (7.4) · *why:* Speed
     is distance over time; ω = 2πν and k = 2π/λ by definition, so the 2π cancel. · *plain:* a crest moves one wavelength
     every period.
- **Result.** $c=\omega/k=\lambda\nu$ (7.4) — *in words:* the phase speed is the speed of any point of fixed phase.
- **Check.** Units: (rad/s)/(rad/m) = m/s ✓. Numbers: λ = 100 m, T = 8 s → ω/k = 0.785/0.0628 = 12.5 m/s = λ/T ✓ (the
  from-scratch crest tracking of C01 gives 12.500 m/s).
- **What it means.** c is a speed of pattern, not of water: the buoy under a crest does not travel. In a dispersive
  medium c depends on k, and the energy travels at a different speed (C09).
- **Traps.** Mixing ν [Hz] and ω [rad/s] loses a 2π (c = λν, but also c = λω/2π). The label n does not change the speed.

### D02 · The exact kinematic condition (7.16) from $(\mathbf n\cdot\mathbf u)_{z=\eta}=\mathbf n\cdot\mathbf U_s$ (7.13) — ★★, 8 steps, in C02 (notebook)
- **Goal.** Turn "the water on the surface stays on the surface" into an equation for φ and η.
- **Start.** $(\mathbf n\cdot\mathbf u)_{z=\eta}=\mathbf n\cdot\mathbf U_s$ (7.13) — *in words:* the water's velocity across the
  surface equals the surface's own normal speed.
- **Plan.** (1) Describe the surface by a level-set function. (2) Build its normal. (3) Choose a simple surface velocity.
  (4) Clear the square root and insert φ.
- **Tools.** level sets and the normal (P75) · moving level set (P132) · chain rule (P91) · R04 (Ch. 4 §4.10).
- **Assumptions.** The surface is a single-valued graph z = η(x, t) — no overturning (step 1); 2-D motion in x–z.
- **Steps.**
  1. *did:* Write the surface as a level set · *tex:* $f(x,z,t)=z-\eta(x,t)=0$ · *why:* The surface is where the height z
     equals η; f vanishes exactly there. Needs a single-valued η (no breaking). ⚠️ Ch. 4 called the level-set function η;
     here η is the height. · *plain:* surface points are those with f = 0.
  2. *did:* Take the gradient of f · *tex:* $\nabla f=-\frac{\partial\eta}{\partial x}\mathbf e_x+\mathbf e_z$ · *why:*
     ∂f/∂x = −∂η/∂x and ∂f/∂z = 1 since η does not depend on z; a level set's gradient is normal to it (P75). · *plain:* a
     vector perpendicular to the surface, pointing up out of the water.
  3. *did:* Normalise to a unit vector · *tex:* $\mathbf n=\frac{-(\partial\eta/\partial x)\mathbf e_x+\mathbf e_z}{\sqrt{
     (\partial\eta/\partial x)^2+1}}$ (7.14) · *why:* Dividing by $\lvert\nabla f\rvert$ gives length 1; the positive z
     part means it points out of the liquid. · *plain:* the unit normal leans back against the slope.
  4. *did:* Choose the surface velocity · *tex:* $\mathbf U_s=\frac{\partial\eta}{\partial t}\mathbf e_z$ (7.15) · *why:*
     Only n·U_s appears in (7.13); following the surface point straight above a fixed x is a valid choice with the right
     normal part. · *plain:* we let the surface point at fixed x go straight up and down.
  5. *did:* Multiply (7.13) by the length of ∇f · *tex:* $(\nabla f\cdot\mathbf u)_{z=\eta}=\nabla f\cdot\mathbf U_s$ · *why:*
     n = ∇f/|∇f| and |∇f| ≥ 1 is never zero, so multiplying both sides is allowed; it removes the square root. · *plain:*
     the same condition without normalising.
  6. *did:* Work out both dot products · *tex:* $\big(-u\frac{\partial\eta}{\partial x}+w\big)_{z=\eta}=\frac{\partial\eta}
     {\partial t}$ · *why:* ∇f·u with u = u e_x + w e_z gives −uη_x + w; ∇f·U_s = e_z·(η_t e_z) = η_t. · *plain:* upward
     velocity minus what the slope turns into the vertical equals the rise rate.
  7. *did:* Insert the velocity potential · *tex:* $\big(\frac{\partial\phi}{\partial z}\big)_{z=\eta}=\frac{\partial\eta}
     {\partial t}+\frac{\partial\eta}{\partial x}\big(\frac{\partial\phi}{\partial x}\big)_{z=\eta}$ (7.16) · *why:* (7.10)
     u = ∂φ/∂x, w = ∂φ/∂z; move the slope term to the right-hand side. · *plain:* the exact kinematic condition in φ.
  8. *did:* Recognise the material derivative · *tex:* $\frac{D(z-\eta)}{Dt}=w-\frac{\partial\eta}{\partial t}-u\frac{\partial
     \eta}{\partial x}=0$ on $z=\eta$ · *why:* D/Dt = ∂/∂t + u∂/∂x + w∂/∂z applied to f = z − η gives step 6: Ch. 4's (4.91)
     with the level-set function z − η. · *plain:* a surface particle keeps f = 0 — it stays on the surface.
- **Result.** $\big(\frac{\partial\phi}{\partial z}\big)_{z=\eta}=\frac{\partial\eta}{\partial t}+\frac{\partial\eta}{\partial x}
  \big(\frac{\partial\phi}{\partial x}\big)_{z=\eta}$ (7.16) — *in words:* exact, nonlinear, and applied on the unknown surface.
- **Check.** Units: every term m/s ✓. Still water (η = 0, φ = 0) satisfies it ✓. A surface rising uniformly (η_x = 0)
  gives w = η_t ✓. The linear wave of C03 misses it by O(a²kω) (C02 code).
- **What it means.** The difficulty of water waves in one line: the condition is nonlinear and lives on a moving,
  unknown boundary. D03 removes both difficulties for small waves.
- **Traps.** Using ch04's η (a level-set function) as if it were the height; forgetting that only the normal part of U_s
  matters; applying it to overturning (breaking) waves.

### D03 · Linearising and moving the kinematic condition: (7.16) → (7.17) → (7.18) — ★★, 8 steps, in C02 (notebook)
- **Goal.** Replace the exact condition (7.16) by a linear one applied on the flat level z = 0, and measure what that
  costs.
- **Start.** $\big(\frac{\partial\phi}{\partial z}\big)_{z=\eta}=\frac{\partial\eta}{\partial t}+\frac{\partial\eta}{\partial x}
  \big(\frac{\partial\phi}{\partial x}\big)_{z=\eta}$ (7.16).
- **Plan.** (1) Size every term. (2) Drop the product of small quantities. (3) Taylor-transfer the left side to z = 0.
  (4) Show both errors are of the same order ka.
- **Tools.** Taylor transfer of a boundary condition (P166) · orders of smallness (P68) · order-of-magnitude scaling
  (P130).
- **Assumptions.** ka ≪ 1 (steps 2, 5); a/H ≪ 1 (step 8); the wave field varies on the scale 1/k.
- **Steps.**
  1. *did:* Size the terms · *tex:* $\frac{\partial\eta}{\partial t}\sim a\omega,\quad\frac{\partial\eta}{\partial x}\sim ka,
     \quad\frac{\partial\phi}{\partial x}\sim a\omega$ · *why:* For η = a cos(kx − ωt) each derivative brings ω or k;
     velocities are aω (C03's (7.27)). Order-of-magnitude scaling (P130). · *plain:* rise rate and velocity are about aω;
     the slope is about ka.
  2. *did:* Size the product term · *tex:* $\frac{\partial\eta}{\partial x}\frac{\partial\phi}{\partial x}\sim ka\cdot a\omega$ ·
     *why:* A product of two small quantities; compared with the kept terms it is smaller by the factor ka ≪ 1 (P68). ·
     *plain:* the slope-times-velocity term is ka times too small to matter.
  3. *did:* Drop it · *tex:* $\big(\frac{\partial\phi}{\partial z}\big)_{z=\eta}\cong\frac{\partial\eta}{\partial t}$ (7.17) ·
     *why:* Linearisation keeps only terms of first order in the amplitude; the relative error is O(ka). · *plain:* the
     water's upward speed at the surface equals the surface's rise rate.
  4. *did:* Taylor-expand the left side about z = 0 · *tex:* $\big(\frac{\partial\phi}{\partial z}\big)_{z=\eta}=\big(\frac{
     \partial\phi}{\partial z}\big)_{z=0}+\eta\big(\frac{\partial^2\phi}{\partial z^2}\big)_{z=0}+\ldots$ · *why:* φ_z is smooth in
     z, so its value at the small height η follows from its value and slope at z = 0 (P166, P26). · *plain:* the velocity
     at the true surface is the velocity at the flat level plus a small correction.
  5. *did:* Size the correction · *tex:* $\eta\frac{\partial^2\phi}{\partial z^2}\sim a\cdot k\cdot a\omega=ka\cdot a\omega$ ·
     *why:* φ_z ∼ aω changes over a vertical distance about 1/k, so φ_zz ∼ k aω; η is about a. · *plain:* the correction is
     also ka times the kept term.
  6. *did:* Compare the two dropped terms · *tex:* $\eta\,\phi_{zz}\sim\eta_x\phi_x\sim ka\cdot a\omega$ · *why:* Consistency:
     step 3 already threw away a term of this size, so keeping this one would be false precision; both go together. ·
     *plain:* the two approximations stand or fall together.
  7. *did:* Keep the leading term · *tex:* $\big(\frac{\partial\phi}{\partial z}\big)_{z=0}\cong\frac{\partial\eta}{\partial t}$
     (7.18) · *why:* Drop the O(ka) correction of step 5; the condition now sits on the known level z = 0. · *plain:* the
     vertical velocity at the still level equals the surface's rise rate.
  8. *did:* State the price · *tex:* $\text{relative error}=O(ka)+O(a/H)$ · *why:* In shallow water φ varies over the depth
     H, not 1/k, so φ_zz ∼ aω/H and the correction is (a/H)·aω; ka = 2πa/λ ≪ 1 is the book's condition. · *plain:*
     for swell (ka ≈ 0.06) the error is a few percent. · *live (B1):* "ηφ_zz = … vs φ_z = …".
- **Result.** $\big(\frac{\partial\phi}{\partial z}\big)_{z=0}\cong\frac{\partial\eta}{\partial t}$ (7.18) — *in words:* linear and
  applied on a fixed level.
- **Check.** C02's residual scan: the exact condition evaluated with the linear solution misses by O((ka)²) in absolute
  terms (slope 2) and O(ka) relative to the kept terms (slope 1) ✓.
- **What it means.** The method of every wave problem in this chapter, of Kelvin–Helmholtz (Ch. 11) and of the linear
  shallow-water equations (Ch. 13): linearise, then transfer the condition to the undisturbed boundary.
- **Traps.** Dropping one O(ka) term and keeping the other; forgetting a/H ≪ 1 in shallow water; thinking the condition
  is still applied on z = η.

### D04 · The linearised dynamic condition: (7.19) → (7.20) → (7.21) — ★★, 8 steps, in C02 (notebook)
- **Goal.** Turn "the surface feels only atmospheric pressure" into a linear condition on φ at z = 0.
- **Start.** $(p)_{z=\eta}=0$ (7.19) and $\frac{\partial\phi}{\partial t}+\tfrac12\lvert\nabla\phi\rvert^2+\frac p\rho+gz=C(t)$
  (4.83) — *in words:* gauge pressure zero on the surface; unsteady Bernoulli everywhere in the water.
- **Plan.** (1) Fix the Bernoulli constant. (2) Drop the square. (3) Evaluate on the surface. (4) Transfer to z = 0.
- **Tools.** unsteady Bernoulli (R06) · gauge pressure (R05) · Taylor transfer (P166) · orders of smallness (P68).
- **Assumptions.** Irrotational, constant density, gravity only (step 1); no surface tension (step 6); ka ≪ 1 (steps 3,
  7).
- **Steps.**
  1. *did:* Write unsteady Bernoulli · *tex:* $\frac{\partial\phi}{\partial t}+\tfrac12\lvert\nabla\phi\rvert^2+\frac p\rho+gz=C(t)$ ·
     *why:* The flow is irrotational with constant density and a conservative body force, so (4.83) holds at every point of
     the water (R06). · *plain:* four terms add up to a function of time only.
  2. *did:* Fix the constant far away · *tex:* $C(t)=0$ · *why:* Far from the waves the water is at rest (φ_t = 0, ∇φ = 0) on
     its still surface z = 0 where the gauge pressure is 0; any leftover C(t) could be absorbed into φ. · *plain:* with gauge
     pressure and z = 0 at the still level, the constant vanishes.
  3. *did:* Size the square term · *tex:* $\tfrac12\lvert\nabla\phi\rvert^2\sim(a\omega)^2,\quad\frac{\partial\phi}{\partial t}
     \sim\frac{a\omega^2}{k}$ · *why:* velocities are aω and φ ∼ aω/k (7.26), so φ_t ∼ aω²/k; their ratio is ka. · *plain:*
     the kinetic term is ka times smaller.
  4. *did:* Drop it · *tex:* $\frac{\partial\phi}{\partial t}+\frac p\rho+gz\cong0$ (7.20) · *why:* A product of two
     first-order quantities is second order; linearisation drops it (P68). · *plain:* linear Bernoulli.
  5. *did:* Evaluate on the surface · *tex:* $\big(\frac{\partial\phi}{\partial t}\big)_{z=\eta}+\frac{(p)_{z=\eta}}{\rho}+g\eta
     \cong0$ · *why:* (7.20) holds everywhere in the water, including just below the surface, where z = η, so gz becomes gη
     (not zero). · *plain:* at the surface, pressure, height and φ_t balance.
  6. *did:* Apply the free-surface condition · *tex:* $\big(\frac{\partial\phi}{\partial t}\big)_{z=\eta}\cong-g\eta$ · *why:*
     (7.19): with no surface tension the gauge pressure just below the surface is zero, so the p term drops. · *plain:* the
     surface height sets the rate of change of φ there.
  7. *did:* Taylor-transfer to z = 0 · *tex:* $\big(\frac{\partial\phi}{\partial t}\big)_{z=\eta}=\big(\frac{\partial\phi}{
     \partial t}\big)_{z=0}+\eta\big(\frac{\partial^2\phi}{\partial z\partial t}\big)_{z=0}+\ldots$ · *why:* Same move as D03 step 4;
     the correction ηφ_tz ∼ a·aω² is ka times φ_t ∼ aω²/k. · *plain:* moving to the flat level costs O(ka).
  8. *did:* Keep the leading term · *tex:* $\big(\frac{\partial\phi}{\partial t}\big)_{z=0}\cong-g\eta$ (7.21) · *why:* Drop
     the O(ka) correction, consistent with (7.18): both surface conditions now live on z = 0. · *plain:* the linear dynamic
     condition.
- **Result.** $\big(\frac{\partial\phi}{\partial t}\big)_{z=\eta}\cong\big(\frac{\partial\phi}{\partial t}\big)_{z=0}\cong-g\eta$
  (7.21) — *in words:* where the surface is high, φ is decreasing.
- **Check.** Units: m²/s² on both sides ✓. For the wave of C03, φ_t(0) = −(aω²/k)coth kH cos(…) and −gη agree exactly when
  (7.28) holds (D06) ✓.
- **What it means.** Together with (7.18) this closes the linear problem: two conditions on z = 0 for one unknown surface
  — one will give φ, the other ω(k).
- **Traps.** Setting gz = 0 at the surface (it is gη); the sign (φ_t = −gη); forgetting that surface tension would add a
  pressure (C07).

### D05 · The velocity potential (7.26) from the separable trial (7.22) — ★★, 11 steps, in C03 (notebook · `dispersion_relation`)
- **Goal.** Find the flow under a small sinusoidal surface wave on water of depth H.
- **Start.** $\nabla^2\phi=0$ (7.11), $\frac{\partial\phi}{\partial z}=0$ at $z=-H$ (7.12), $\big(\frac{\partial\phi}{\partial
  z}\big)_{z=0}=\frac{\partial\eta}{\partial t}$ (7.18), with $\eta=a\cos(kx-\omega t)$ (7.2) — *in words:* Laplace in a strip,
  no flow through the bottom, the surface moves with the water.
- **Plan.** (1) Try a product of a depth profile and a travelling sine. (2) Laplace gives an ODE for the profile. (3) The
  bottom fixes its shape, the surface its size. (4) Rewrite with cosh.
- **Tools.** separation of variables for a PDE (P167) · e^{λz} trial (P44) · hyperbolic functions (P168) · two linear
  equations.
- **Assumptions.** Linear problem of C02; flat bottom; a single sinusoid. The dynamic condition (7.21) is **not** used yet
  (it gives ω(k) in D06).
- **Steps.**
  1. *did:* Try a product form · *tex:* $\phi=f(z)\sin(kx-\omega t)$ (7.22) · *why:* (7.18) needs φ_z ∝ η_t ∝ sin and (7.21)
     needs φ_t ∝ η ∝ cos: a sine of the phase does both. Separation of variables (P167). · *plain:* every depth moves with
     the same travelling shape, with its own strength f(z).
  2. *did:* Put it into Laplace's equation · *tex:* $f''(z)\sin(kx-\omega t)-k^2f(z)\sin(kx-\omega t)=0$ · *why:*
     ∂²/∂x² of sin(kx − ωt) is −k² sin(kx − ωt); ∂²/∂z² acts on f only. · *plain:* the curvature in x must cancel the
     curvature in z.
  3. *did:* Divide by the sine · *tex:* $\frac{d^2f}{dz^2}-k^2f=0$ (7.23) · *why:* The equation holds for all x and t, in
     particular where the sine is not zero; dividing there leaves an ODE in z alone. · *plain:* the depth profile obeys a
     simple ODE.
  4. *did:* Solve the ODE · *tex:* $f(z)=Ae^{kz}+Be^{-kz}$ (7.23) · *why:* The trial e^{λz} (P44) gives λ² = k², λ = ±k —
     real because of the + sign: growth and decay in z, not oscillation. · *plain:* the strength changes exponentially with
     depth.
  5. *did:* Apply the bottom condition · *tex:* $k(Ae^{-kH}-Be^{kH})=0\ \Rightarrow\ B=Ae^{-2kH}$ (7.24) · *why:* (7.12):
     f′(−H) = 0 for all x and t; multiply by e^{kH}/k. · *plain:* the bottom fixes the mix of the two exponentials.
  6. *did:* Apply the kinematic condition · *tex:* $k(A-B)\sin(kx-\omega t)=\omega a\sin(kx-\omega t)$ (7.25) · *why:* (7.18):
     φ_z at z = 0 is k(A − B) sin(…), and η_t = aω sin(kx − ωt) from (7.2). · *plain:* the surface rises as fast as the
     water under it.
  7. *did:* Cancel the sine and insert B · *tex:* $A(1-e^{-2kH})=\frac{\omega a}{k}$ · *why:* Matching coefficients (true for
     all x, t); B = Ae^{−2kH} from step 5. · *plain:* one equation for A.
  8. *did:* Solve for A and B · *tex:* $A=\frac{a\omega}{k(1-e^{-2kH})},\quad B=\frac{a\omega e^{-2kH}}{k(1-e^{-2kH})}$ · *why:*
     1 − e^{−2kH} > 0 for kH > 0, so we may divide; B = Ae^{−2kH}. · *plain:* both constants from a, ω, k, H.
  9. *did:* Regroup f as a cosh · *tex:* $f=Ae^{-kH}\big(e^{k(z+H)}+e^{-k(z+H)}\big)=2Ae^{-kH}\cosh k(z+H)$ · *why:*
     Ae^{kz} + Ae^{−2kH}e^{−kz}: factor Ae^{−kH}; cosh u = (eᵘ + e⁻ᵘ)/2 (P168). · *plain:* the depth profile is a cosh
     measured up from the bottom. · *set (E1):* the amber profile lights up.
  10. *did:* Simplify the prefactor · *tex:* $2Ae^{-kH}=\frac{a\omega}{k}\,\frac{2e^{-kH}}{1-e^{-2kH}}=\frac{a\omega}{k\sinh kH}$ ·
      *why:* Multiply top and bottom by e^{kH}: 2/(e^{kH} − e^{−kH}) = 1/sinh kH. · *plain:* the strength at the bottom
      of the cosh.
  11. *did:* Assemble φ · *tex:* $\phi=\frac{a\omega}{k}\frac{\cosh k(z+H)}{\sinh kH}\sin(kx-\omega t)$ (7.26) · *why:* Put
      steps 9–10 into (7.22). ω is still free: the dynamic condition will fix it (D06). · *plain:* the flow under a linear
      wave. · *live (E1):* "φ amplitude aω/(k sinh kH) = … m²/s".
- **Result.** $\phi=\frac{a\omega}{k}\frac{\cosh k(z+H)}{\sinh kH}\sin(kx-\omega t)$ (7.26), and by one derivative each
  $u=a\omega\frac{\cosh k(z+H)}{\sinh kH}\cos(kx-\omega t)$, $w=a\omega\frac{\sinh k(z+H)}{\sinh kH}\sin(kx-\omega t)$ (7.27).
- **Check.** Units: aω/k [m²/s] ✓. ∇²φ = (k² − k²)φ = 0 ✓; φ_z(−H) ∝ sinh 0 = 0 ✓; φ_z(0) = aω sin(…) = η_t ✓ (the
  `surface_wave_sympy` residuals are all 0; parity with Ch. 4's test field to 1e-12).
- **What it means.** The motion decays from the surface on the scale 1/k and is squeezed by the bottom when kH is small;
  this cosh ratio reappears in the pressure (7.31), the orbits (7.35) and the energy (7.38).
- **Traps.** Choosing cos instead of sin in (7.22); expecting sin/cos in z (the + sign gives exponentials); using the
  dynamic condition too early; forgetting that sinh kH → 0 in very shallow water (the ratio stays finite).

### D06 · The dispersion relation $\omega=\sqrt{gk\tanh kH}$ (7.28) — ★★, 7 steps, in C03 (notebook · `dispersion_relation`)
- **Goal.** Use the last condition — the free surface feels no pressure — to find the one frequency that goes with each
  wavenumber.
- **Start.** $\phi=\frac{a\omega}{k}\frac{\cosh k(z+H)}{\sinh kH}\sin(kx-\omega t)$ (7.26) and $\big(\frac{\partial\phi}{\partial
  t}\big)_{z=0}\cong-g\eta$ (7.21).
- **Plan.** (1) Compute φ_t at z = 0. (2) Insert into (7.21). (3) Match coefficients and solve for ω.
- **Tools.** hyperbolic functions and coth (P168) · matching coefficients (gloss) · D05.
- **Assumptions.** As D05; ω > 0 for a right-going wave (step 6).
- **Steps.**
  1. *did:* Differentiate φ in time · *tex:* $\frac{\partial\phi}{\partial t}=-\frac{a\omega^2}{k}\frac{\cosh k(z+H)}{\sinh kH}
     \cos(kx-\omega t)$ · *why:* d/dt sin(kx − ωt) = −ω cos(kx − ωt); the depth factor does not depend on t. · *plain:* φ_t
     is a cosine, in step with η.
  2. *did:* Evaluate at z = 0 · *tex:* $\big(\frac{\partial\phi}{\partial t}\big)_{z=0}=-\frac{a\omega^2}{k}\frac{\cosh kH}{\sinh kH}
     \cos(kx-\omega t)$ · *why:* (7.21) is applied on the still level; cosh k(0 + H) = cosh kH. · *plain:* the value at the
     surface.
  3. *did:* Insert into the dynamic condition · *tex:* $-\frac{a\omega^2}{k}\coth kH\cos(kx-\omega t)=-ag\cos(kx-\omega t)$ ·
     *why:* (7.21) with η = a cos(kx − ωt) (7.2); coth = cosh/sinh. This is the book's unnumbered line. · *plain:* the
     surface's inertia must balance gravity everywhere.
  4. *did:* Match the coefficients · *tex:* $\frac{\omega^2}{k}\coth kH=g$ · *why:* The identity holds for all x and t, so
     the coefficients of cos(kx − ωt) must be equal; then divide by −a ≠ 0. · *plain:* one equation linking ω and k. ·
     *set (E1):* the dot on ω(k) lights up.
  5. *did:* Solve for ω² · *tex:* $\omega^2=gk\tanh kH$ · *why:* Multiply by k and by tanh kH = 1/coth kH. · *plain:* the
     frequency squared.
  6. *did:* Take the positive root · *tex:* $\omega=\sqrt{gk\tanh kH}$ (7.28) · *why:* −ω gives the same wave moving the
     other way (7.61); for our right-going wave ω > 0. · *plain:* the dispersion relation. · *live (E1):* "ω² = 9.81 × k ×
     tanh kH = …".
  7. *did:* Rewrite with T and λ · *tex:* $T=\sqrt{\frac{2\pi\lambda}{g}\coth\frac{2\pi H}{\lambda}}$ (7.28) · *why:* ω = 2π/T, k =
     2π/λ: T² = 4π²/[(2πg/λ) tanh(2πH/λ)] = (2πλ/g) coth(2πH/λ). · *plain:* the period that goes with each wavelength.
- **Result.** $\omega=\sqrt{gk\tanh kH}$, or $T=\sqrt{\frac{2\pi\lambda}{g}\coth\frac{2\pi H}{\lambda}}$ (7.28).
- **Check.** Units: g k [1/s²] ✓. Deep limit tanh → 1: ω² = gk ✓; shallow tanh kH ≈ kH: ω = k√(gH) ✓. Number: H = 1 m,
  k = 1 rad/m → ω = 2.73 rad/s (C03 tiny example).
- **What it means.** The water will not oscillate at any frequency for a given wavelength: gravity (the restoring force)
  and the inertia of a layer of depth ~min(1/k, H) fix it. Everything in §7.2–§7.5 is read off this relation.
- **Traps.** Cancelling cos(kx − ωt) as if it were a number that could be zero (it is legal because the identity holds for
  all x, t); flipping coth and tanh; dropping one of the two 2π in the T–λ form.

### D07 · Phase speed (7.29) and why longer waves are faster — ★, 5 steps, in C04 (notebook · `dispersion_relation`)
- **Goal.** Get the speed of crests from (7.28) and show it always grows with wavelength.
- **Start.** $\omega=\sqrt{gk\tanh kH}$ (7.28) and $c=\omega/k$ (7.4).
- **Plan.** (1) Divide. (2) Rewrite in λ. (3) Differentiate c² with respect to λ.
- **Tools.** product and chain rules (P38, P49) · hyperbolic functions (P168).
- **Assumptions.** k > 0, H > 0.
- **Steps.**
  1. *did:* Divide ω by k · *tex:* $c=\frac{\omega}{k}=\frac{\sqrt{gk\tanh kH}}{k}$ · *why:* This is the definition of the phase speed, (7.4) c = ω/k; we apply it to the dispersion relation. · *plain:*
     crest speed from the dispersion relation.
  2. *did:* Bring k under the root · *tex:* $c=\sqrt{\frac gk\tanh kH}$ (7.29) · *why:* k = √(k²) for k > 0, and gk/k² = g/k.
     · *plain:* the speed depends on g/k and on kH.
  3. *did:* Rewrite in wavelength · *tex:* $c=\sqrt{\frac{g\lambda}{2\pi}\tanh\frac{2\pi H}{\lambda}}$ (7.29) · *why:* k = 2π/λ,
     both in g/k and inside the tanh. · *plain:* the same speed in terms of λ.
  4. *did:* Differentiate c² with respect to λ · *tex:* $\frac{d}{d\lambda}\Big[\lambda\tanh\frac{2\pi H}{\lambda}\Big]=
     \frac{\tfrac12\sinh2x-x}{\cosh^2x},\ x=\frac{2\pi H}{\lambda}$ · *why:* Product and chain rules: tanh x + λ sech²x·(−x/λ)
     = (sinh x cosh x − x)/cosh²x; sinh 2x = 2 sinh x cosh x. · *plain:* the slope of c² against λ.
  5. *did:* Read its sign · *tex:* $\tfrac12\sinh 2x>x\ \Rightarrow\ \frac{dc}{d\lambda}>0$ · *why:* sinh y > y for every y > 0
     (its series has only positive terms), with y = 2x. So c grows with λ at any depth: dispersive. · *plain:* longer
     waves outrun shorter ones. · *watch (E1):* "the black curve rises to the right".
- **Result.** $c=\sqrt{\frac gk\tanh kH}=\sqrt{\frac{g\lambda}{2\pi}\tanh\frac{2\pi H}{\lambda}}$ (7.29), increasing with λ.
- **Check.** Units m/s ✓. H = 10 m: λ = 50 m → 8.15 m/s, λ = 500 m → 9.88 m/s — faster ✓, and never above √(gH) = 9.90 m/s.
- **What it means.** An impulse spreads into a train with the longest waves in front (the stone in a pond, C09). For
  long waves the growth stops at √(gH): no dispersion (D09).
- **Traps.** √(g tanh kH/k²)·k is not the same as √((g/k) tanh kH)·… — keep one k under the root; change k in the tanh
  too when writing in λ.

### D08 · Deep-water limits (7.45)–(7.48) — ★, 7 steps, in C04 (notebook · `dispersion_relation`)
- **Goal.** Simplify speed, orbits, velocity and pressure when the water is much deeper than the wave is long, and say
  how deep is deep.
- **Start.** (7.29) $c=\sqrt{\frac gk\tanh kH}$, the depth ratios $\frac{\cosh k(z+H)}{\sinh kH}$, $\frac{\sinh k(z+H)}{\sinh kH}$
  of (7.27), (7.35) and $\frac{\cosh k(z+H)}{\cosh kH}$ of (7.31).
- **Plan.** (1) Let kH be large in the speed. (2) Size the error at kH = 2. (3) Rewrite the depth ratios with exponentials
  and take the limit.
- **Tools.** hyperbolic functions and their large-x limits (P168, N37).
- **Assumptions.** kH > 2 (steps 1–3); points not too close to the bottom, z + H ≫ 1/k (step 5).
- **Steps.**
  1. *did:* Let kH grow · *tex:* $\tanh kH=\frac{1-e^{-2kH}}{1+e^{-2kH}}\to1$ · *why:* Divide sinh and cosh by e^{kH}/2; the
     small exponentials vanish as kH → ∞ (P168). · *plain:* the bottom stops mattering.
  2. *did:* Simplify the speed · *tex:* $c=\sqrt{g/k}=\sqrt{g\lambda/2\pi}$ (7.45) · *why:* Substitute tanh kH = 1 in (7.29).
     · *plain:* deep-water speed depends on wavelength only.
  3. *did:* Size the error at kH = 2 · *tex:* $1-\sqrt{\tanh2}=1-\sqrt{0.964}=0.018$ · *why:* c/c_deep = √(tanh kH); kH = 2 ⇔
     H = 2λ/2π = 0.32λ. · *plain:* deeper than a third of a wavelength is 'deep' to 2 %. · *live (E1):* "tanh kH = …: c within
     … % of √(g/k)".
  4. *did:* Rewrite the depth ratio exactly · *tex:* $\frac{\cosh k(z+H)}{\sinh kH}=e^{kz}\,\frac{1+e^{-2k(z+H)}}{1-e^{-2kH}}$ ·
     *why:* Multiply top and bottom by 2e^{−kH}: 2e^{−kH}cosh k(z + H) = e^{kz} + e^{−kz−2kH}, 2e^{−kH}sinh kH = 1 −
     e^{−2kH}. · *plain:* an exact form that shows e^{kz}.
  5. *did:* Take the deep limit · *tex:* $\frac{\cosh k(z+H)}{\sinh kH}\approx\frac{\sinh k(z+H)}{\sinh kH}\approx e^{kz}$ · *why:* For
     kH ≫ 1 and z + H ≫ 1/k both small exponentials vanish; the sinh ratio has the same limit. · *plain:* motion decays like
     e^{kz} below the surface.
  6. *did:* Apply to the orbits, velocity and pressure · *tex:* $\xi\cong-ae^{kz_0}\sin(kx_0-\omega t),\ \zeta\cong ae^{kz_0}\cos(kx_0-
     \omega t)$ (7.46); $p'=\rho gae^{kz}\cos(kx-\omega t)$ (7.48) · *why:* The same ratios sit in (7.35), (7.27) and (7.31)
     (cosh k(z + H)/cosh kH → e^{kz} too); equal semi-axes make the orbits circles, and (7.47) u, w = aωe^{kz}(cos, sin). ·
     *plain:* circles, and velocity and pressure fading like e^{kz}.
  7. *did:* Evaluate half a wavelength down · *tex:* $e^{-k\lambda/2}=e^{-\pi}=0.043$ · *why:* z = −λ/2 means kz = −π. ·
     *plain:* at depth λ/2 only 4 % of the surface pressure and orbit size remain.
- **Result.** $c=\sqrt{g/k}$ (7.45); orbits (7.46) circles of radius ae^{kz₀}; $u,w=a\omega e^{kz}(\cos,\sin)(kx-\omega t)$ (7.47);
  $p'=\rho gae^{kz}\cos(kx-\omega t)$ (7.48).
- **Check.** kH = 3: error 0.25 % ✓ (E1 preset). e^{kz} at z = −λ/2 = 4.3 % ✓. Units unchanged.
- **What it means.** Wind waves in the ocean (λ ~ 100 m, H ~ 4 km) never feel the sea floor until they reach the shelf;
  bottom pressure gauges filter them out.
- **Traps.** kH > 2 is 1.8 % error, not zero; the e^{kz} form fails near the bottom when H is finite.

### D09 · Shallow-water limits (7.49)–(7.52) — ★, 7 steps, in C04 (notebook · `dispersion_relation`)
- **Goal.** Simplify speed, orbits, velocity and pressure when the wave is much longer than the depth, and size the
  error.
- **Start.** (7.29), (7.35), (7.27), (7.31) as in D08.
- **Plan.** (1) Expand tanh for small kH. (2) Keep the leading term for the speed and the next for the error. (3) Apply
  the small-argument forms to orbits, velocities and pressure.
- **Tools.** small-argument forms of cosh, sinh, tanh (P168, N37) · Taylor series (P117) · √(1 − ε) ≈ 1 − ε/2 (P26).
- **Assumptions.** kH ≪ 1 (all steps); H < 0.07λ for 3 % (step 4).
- **Steps.**
  1. *did:* Expand tanh for small kH · *tex:* $\tanh kH\approx kH-\tfrac13(kH)^3$ · *why:* Taylor series of tanh about 0
     (P117). · *plain:* for long waves tanh kH is almost kH.
  2. *did:* Keep the leading term · *tex:* $c=\sqrt{\frac gk\,kH}=\sqrt{gH}$ (7.49) · *why:* Keeping only the first term of step 1 is allowed for kH ≪ 1; then (g/k)·kH = gH and k cancels. · *plain:* the
     speed depends only on the depth — every long wave travels together.
  3. *did:* Keep the next term for the error · *tex:* $\frac{c}{\sqrt{gH}}=\sqrt{1-\tfrac13(kH)^2}\approx1-\tfrac16(kH)^2$ · *why:*
     (g/k)(kH − (kH)³/3) = gH(1 − (kH)²/3), and √(1 − ε) ≈ 1 − ε/2. · *plain:* the relative error is (kH)²/6. · *live (E1):*
     "1 − (kH)²/6 = …".
  4. *did:* Evaluate at H = 0.07λ · *tex:* $kH=2\pi\times0.07=0.44,\quad\sqrt{\tanh kH/kH}=0.969$ · *why:* The exact ratio
     √(0.4134/0.4398) gives 3 % error (the two-term estimate 0.968). · *plain:* waves 14 times longer than the depth are
     'shallow' to 3 %.
  5. *did:* Simplify the orbits · *tex:* $\xi\cong-\frac{a}{kH}\sin(kx_0-\omega t),\ \zeta\cong a\big(1+\frac{z_0}H\big)\cos(kx_0-\omega t)$
     (7.50) · *why:* In (7.35) put cosh k(z₀ + H) ≈ 1, sinh k(z₀ + H) ≈ k(z₀ + H), sinh kH ≈ kH. · *plain:* flat ellipses: the
     same width a/kH at every depth, height falling to 0 at the bottom.
  6. *did:* Simplify the velocities · *tex:* $u=\frac{a\omega}{kH}\cos(kx-\omega t),\ w=a\omega\big(1+\frac zH\big)\sin(kx-\omega t)$
     (7.51) · *why:* The same substitutions in (7.27); w/u ≤ kH ≪ 1. · *plain:* the motion is almost horizontal and uniform
     with depth.
  7. *did:* Simplify the pressure · *tex:* $p'=\rho ga\cos(kx-\omega t)=\rho g\eta$ (7.52) · *why:* cosh k(z + H)/cosh kH ≈ 1/1 in
     (7.31); a cos(kx − ωt) = η by (7.2). · *plain:* the pressure is just the weight of the extra water above: hydrostatic.
- **Result.** $c=\sqrt{gH}$ (7.49), flat-ellipse orbits (7.50), depth-uniform u (7.51), hydrostatic $p'=\rho g\eta$ (7.52).
- **Check.** H = 4 km: √(gH) = 198 m/s ✓ (a 200 km tsunami has 197.6 m/s from (7.29)). Matches Ch. 4's bore speed for a small
  step ✓.
- **What it means.** Shallow water is hydrostatic and nondispersive: the model of tides, tsunamis and storm surges (Ch. 13).
- **Traps.** The semi-major axis a/kH is large, not small; c/√(gH) − 1 is quadratic in kH, so 'shallow' can be quite
  generous; p′ = ρgη is independent of z.

### D10 · Linearised particle orbits: (7.33) → (7.34a, b) → (7.35a, b) — ★★, 9 steps, in C05 (notebook · `particle_orbits`)
- **Goal.** Find the path of one water parcel under the wave, to first order in the amplitude.
- **Start.** $\frac{dx_p}{dt}=u(x_p,z_p,t)$, $\frac{dz_p}{dt}=w(x_p,z_p,t)$ (7.32) with u, w from (7.27).
- **Plan.** (1) Insert the velocities. (2) Split the position into mean + small excursion. (3) Freeze the right sides at the
  mean position. (4) Integrate in time.
- **Tools.** path lines (R08) · multivariable Taylor (P98) · antiderivatives of sin and cos.
- **Assumptions.** a ≪ λ (steps 3–4); the mean position is the time average of the path (step 8).
- **Steps.**
  1. *did:* Insert the velocities · *tex:* $\frac{dx_p}{dt}=a\omega\frac{\cosh k(z_p+H)}{\sinh kH}\cos(kx_p-\omega t)$ and the same
     for z_p with sinh, sin (7.33) · *why:* (7.32) says a parcel moves with the local velocity; (7.27) gives it. · *plain:*
     exact but nonlinear equations for the path.
  2. *did:* Split the position · *tex:* $x_p=x_0+\xi(t),\quad z_p=z_0+\zeta(t)$ · *why:* We expect the parcel to stay near a
     fixed mean position (x₀, z₀), with excursions ξ, ζ of order a, small compared with λ. · *plain:* the parcel wanders
     slightly round a fixed point.
  3. *did:* Size the velocity change over an excursion · *tex:* $u(x_0+\xi,\ldots)-u(x_0,\ldots)\approx\xi\frac{\partial u}{\partial
     x}\sim a\cdot k\cdot a\omega$ · *why:* First-order Taylor (P98): the field changes on the scale 1/k, so the difference
     is ka times u itself. · *plain:* moving by ξ changes the velocity felt by only ka of it.
  4. *did:* Freeze the argument · *tex:* $\frac{d\xi}{dt}\cong a\omega\frac{\cosh k(z_0+H)}{\sinh kH}\cos(kx_0-\omega t)$ (7.34a) ·
     *why:* Drop the O(ka) correction of step 3 (linearisation); dx_p/dt = dξ/dt because x₀ is fixed. · *plain:* the parcel
     feels the velocity at its mean position. · *watch (E2):* "the black arrow is the neglected ξ ∂u/∂x".
  5. *did:* The same for the vertical · *tex:* $\frac{d\zeta}{dt}\cong a\omega\frac{\sinh k(z_0+H)}{\sinh kH}\sin(kx_0-\omega t)$ (7.34b)
     · *why:* Same move on the z-equation of (7.33). · *plain:* the vertical velocity at the mean position.
  6. *did:* Integrate the horizontal equation · *tex:* $\xi=-a\frac{\cosh k(z_0+H)}{\sinh kH}\sin(kx_0-\omega t)+\xi_c$ · *why:*
     ∫cos(kx₀ − ωt)dt = −sin(kx₀ − ωt)/ω, since d/dt sin(kx₀ − ωt) = −ω cos(…); aω/ω = a. · *plain:* the horizontal
     excursion oscillates.
  7. *did:* Integrate the vertical equation · *tex:* $\zeta=a\frac{\sinh k(z_0+H)}{\sinh kH}\cos(kx_0-\omega t)+\zeta_c$ · *why:*
     ∫sin(kx₀ − ωt)dt = cos(kx₀ − ωt)/ω. · *plain:* the vertical excursion oscillates a quarter period out of step. ·
     *live (E2):* "ξ amplitude = … m".
  8. *did:* Set the constants to zero · *tex:* $\xi_c=\zeta_c=0$ ⇒ (7.35a, 7.35b) · *why:* (x₀, z₀) is by definition the
     mean position; the oscillating parts already average to zero, so any constant would shift the mean. · *plain:* the
     parcel circles its own mean position.
  9. *did:* Check nothing grows · *tex:* $\lvert\xi\rvert\le a\coth kH,\ \lvert\zeta\rvert\le a$ · *why:* No term grows with t (no
     secular term), so ξ, ζ stay small and the freezing of step 4 stays valid forever. The dropped ξu_x returns at second
     order as Stokes drift (C12). · *plain:* the approximation does not wear out.
- **Result.** $\xi\cong-a\frac{\cosh k(z_0+H)}{\sinh kH}\sin(kx_0-\omega t)$, $\zeta\cong a\frac{\sinh k(z_0+H)}{\sinh kH}\cos(kx_0-
  \omega t)$ (7.35a, 7.35b).
- **Check.** Units: m ✓. dξ/dt of the result returns (7.34a) ✓. RK4 on (7.34) (C05 from scratch) matches to 1e-6 ✓.
- **What it means.** Every parcel oscillates about a fixed point; parcels of one vertical column are in phase.
- **Traps.** The minus sign in ξ from ∫cos; putting the constants of integration anywhere but zero; thinking the dropped
  term is always negligible (its average is the drift, C12).

### D11 · The orbits are ellipses (7.36): axes, foci, sense — ★, 6 steps, in C05 (notebook · `particle_orbits`)
- **Goal.** Show that each parcel traces an ellipse, find its axes and foci, and its direction of travel.
- **Start.** (7.35a, 7.35b) from D10.
- **Plan.** (1) Name the amplitudes. (2) Normalise, square, add. (3) Read off the shape and the sense.
- **Tools.** sin² + cos² = 1 · ellipse geometry (P104) · cosh² − sinh² = 1 (P168).
- **Assumptions.** Linear orbits (D10); z₀ above the bottom.
- **Steps.**
  1. *did:* Name the semi-axes · *tex:* $A=a\frac{\cosh k(z_0+H)}{\sinh kH},\quad B=a\frac{\sinh k(z_0+H)}{\sinh kH}$ · *why:* The
     amplitudes of ξ and ζ in (7.35); they depend on the depth only. · *plain:* half-width and half-height of the loop.
  2. *did:* Divide by the semi-axes · *tex:* $\frac\xi A=-\sin(kx_0-\omega t),\quad\frac\zeta B=\cos(kx_0-\omega t)$ · *why:* A, B > 0
     above the bottom, so dividing is allowed. · *plain:* two normalised oscillations.
  3. *did:* Square and add · *tex:* $\frac{\xi^2}{A^2}+\frac{\zeta^2}{B^2}=\sin^2+\cos^2=1$ (7.36) · *why:* sin² + cos² = 1
     removes the phase — the time. · *plain:* every parcel traces an ellipse.
  4. *did:* Compare the axes · *tex:* $\frac BA=\tanh k(z_0+H)\le1$ · *why:* The ratio of step 1's expressions; tanh < 1, → 1
     in deep water (circles), → 0 at the bottom (a flat line). · *plain:* the loops flatten with depth. · *live (E2):* "B/A
     = …".
  5. *did:* Find the focal distance · *tex:* $\sqrt{A^2-B^2}=\frac{a}{\sinh kH}\sqrt{\cosh^2-\sinh^2}=\frac{a}{\sinh kH}$ · *why:*
     An ellipse with A ≥ B has foci √(A² − B²) from its centre; cosh² − sinh² = 1 (P168). · *plain:* the foci are the same
     distance apart at every depth.
  6. *did:* Find the sense of rotation · *tex:* $\zeta=+B\Rightarrow\frac{d\xi}{dt}=+a\omega\frac{\cosh k(z_0+H)}{\sinh kH}>0$ · *why:*
     ζ = B when cos(kx₀ − ωt) = 1; then (7.34a) gives a forward velocity. Forward at the top, backward at the bottom:
     clockwise with x right and z up. · *plain:* parcels go round clockwise for a wave travelling to the right. · *watch
     (E2):* "the selected parcel runs clockwise".
- **Result.** $\xi^2/A^2+\zeta^2/B^2=1$ (7.36) with the axes of step 1; foci a/sinh kH from the centre; clockwise.
- **Check.** kH = 3 at the surface: A = 1.005a, B = a ✓ (C05 code); kH = 0.3: A = 3.43a ✓.
- **What it means.** The picture of what water does under a wave; the constant focal distance and the in-phase columns
  are visible in the A1 animation.
- **Traps.** Squaring before dividing by the axes; calling the rotation anticlockwise (it depends on the direction of
  travel); forgetting that B = 0 at the bottom (the bottom parcel only slides).

### D12 · The kinetic energy: (7.38) → (7.39) $E_k=\tfrac12\rho g\overline{\eta^2}$ — ★★, 10 steps, in C06 (notebook)
- **Goal.** Find the kinetic energy of the orbiting water per square metre of sea surface.
- **Start.** $E_k=\frac{\rho}{2\lambda}\int_0^\lambda\int_{-H}^0(u^2+w^2)\,dz\,dx$ — *in words:* ½ρ|u|² summed over the depth, averaged
  over a wavelength.
- **Plan.** (1) Insert (7.27). (2) Do the x-averages. (3) Do the z-integrals with hyperbolic identities. (4) Use (7.28).
- **Tools.** mean of cos² and sin² (P151) · hyperbolic identities (P168) · definite integral (P27) · (7.28).
- **Assumptions.** Linear wave; integrate to z = 0, not η (step 2).
- **Steps.**
  1. *did:* Write the energy per unit area · *tex:* $E_k=\frac{\rho}{2\lambda}\int_0^\lambda\int_{-H}^0(u^2+w^2)\,dz\,dx$ · *why:*
     ½ρ(u² + w²) is the kinetic energy per unit volume; integrate over the depth and average over one wavelength. ·
     *plain:* kinetic energy under one square metre.
  2. *did:* Stop the depth integral at z = 0 · *tex:* $\int_{-H}^{\eta}\approx\int_{-H}^{0}$ · *why:* The slab between 0 and η is
     O(a) thick with u² ∼ (aω)²: an O(a³) change, beyond the O(a²) accuracy of a linear energy. · *plain:* consistent
     with linear theory.
  3. *did:* Insert the velocities · *tex:* $E_k=\frac{\rho\omega^2}{2\sinh^2kH}\Big[\frac1\lambda\int_0^\lambda a^2\cos^2dx\int_{-H}^0
     \cosh^2k(z+H)dz+\frac1\lambda\int_0^\lambda a^2\sin^2dx\int_{-H}^0\sinh^2k(z+H)dz\Big]$ (7.38) · *why:* u² and w² from (7.27);
     each splits into an x-part and a z-part. · *plain:* horizontal and vertical motions contribute separately.
  4. *did:* Do the x-averages · *tex:* $\frac1\lambda\int_0^\lambda a^2\cos^2dx=\frac1\lambda\int_0^\lambda a^2\sin^2dx=\frac{a^2}2=
     \overline{\eta^2}$ · *why:* Over a whole wavelength the mean of cos² and of sin² is ½ (P151); η² = a² cos² has the same
     mean. · *plain:* both averages are the mean-square surface height.
  5. *did:* Integrate cosh² over the depth · *tex:* $\int_{-H}^0\cosh^2k(z+H)dz=\frac H2+\frac{\sinh2kH}{4k}$ · *why:* cosh²u =
     (1 + cosh 2u)/2 (P168); substitute s = z + H from 0 to H. · *plain:* the depth weight of the horizontal motion.
  6. *did:* Integrate sinh² over the depth · *tex:* $\int_{-H}^0\sinh^2k(z+H)dz=-\frac H2+\frac{\sinh2kH}{4k}$ · *why:* sinh²u =
     (cosh 2u − 1)/2, same substitution. · *plain:* the depth weight of the vertical motion.
  7. *did:* Add with equal weights · *tex:* $E_k=\frac{\rho\omega^2\overline{\eta^2}}{2\sinh^2kH}\cdot\frac{\sinh2kH}{2k}$ · *why:* Both
     x-averages equal the mean square, so the ±H/2 cancel and the two sinh 2kH/4k add. · *plain:* the depth-integrated kinetic
     energy.
  8. *did:* Use the double-angle identity · *tex:* $E_k=\frac{\rho\omega^2\overline{\eta^2}\cosh kH}{2k\sinh kH}$ · *why:* sinh 2kH =
     2 sinh kH cosh kH (P168); one sinh cancels. · *plain:* a simpler form.
  9. *did:* Insert the dispersion relation · *tex:* $E_k=\frac{\rho\,gk\tanh kH\,\overline{\eta^2}}{2k\tanh kH}$ · *why:* (7.28)
     ω² = gk tanh kH; cosh/sinh = 1/tanh. · *plain:* gravity enters.
  10. *did:* Cancel · *tex:* $E_k=\tfrac12\rho g\overline{\eta^2}$ (7.39) · *why:* k and tanh kH cancel exactly. · *plain:*
      kinetic energy depends only on the mean-square surface height.
- **Result.** $E_k=\tfrac12\rho g\overline{\eta^2}$ (7.39) — *in words:* independent of depth and wavelength.
- **Check.** Units: ρg·m² = J/m² ✓. `wave_energy(method="quad")` (dblquad of the first form) agrees to 1e-8 for kH = 0.5,
  1, 3 ✓; midpoint sum from scratch to 1e-4 ✓.
- **What it means.** Gravity waves carry exactly as much kinetic energy as D13 will find potential energy: equipartition.
- **Traps.** Forgetting the H/2 terms cancel only because both x-averages are equal; dropping the ½ of the double-angle;
  integrating to η.

### D13 · Potential energy (7.40), equipartition (7.41) and $E=\tfrac12\rho ga^2$ (7.42) — ★, 7 steps, in C06 (notebook)
- **Goal.** Find the energy stored in the lifted and lowered water, compare it with D12, and add them.
- **Start.** Potential energy per unit volume ρgz — *in words:* water higher up stores more energy.
- **Plan.** (1) Define E_p as wavy minus flat. (2) Cancel the common part of the integrals. (3) Integrate. (4) Compare and
  add.
- **Tools.** splitting an integral · definite integral (P27) · mean of cos² (P151) · D12.
- **Assumptions.** Linear wave; sinusoidal η for the last step.
- **Steps.**
  1. *did:* Define E_p as a difference · *tex:* $E_p=\frac{\rho g}{\lambda}\int_0^\lambda\!\int_{-H}^{\eta}z\,dz\,dx-\frac{\rho g}
     {\lambda}\int_0^\lambda\!\int_{-H}^{0}z\,dz\,dx$ · *why:* Potential energy is measured from the flat, undisturbed state:
     the work needed to deform the surface (ρgz per unit volume, per unit length along the crest). · *plain:* wave PE =
     wavy water minus flat water.
  2. *did:* Split the first z-integral · *tex:* $\int_{-H}^{\eta}z\,dz=\int_{-H}^{0}z\,dz+\int_0^{\eta}z\,dz$ · *why:* Integrals are
     additive over adjacent intervals; the first piece cancels the flat-state integral. · *plain:* only the water between
     the still level and the surface counts.
  3. *did:* Integrate · *tex:* $\int_0^{\eta}z\,dz=\tfrac12\eta^2$ · *why:* Antiderivative z²/2 from 0 to η; the H² terms have
     cancelled; for η < 0 it is also +η²/2 (removing water below the level also stores energy). · *plain:* each column
     stores ρgη²/2 per unit area.
  4. *did:* Write E_p · *tex:* $E_p=\frac{\rho g}{2\lambda}\int_0^\lambda\eta^2dx$ (7.40) · *why:* Put step 3 into step 1. The
     column swap of Fig. 7.6 agrees: over half a wavelength, lifting a column of mass ρη dx by η costs ρgη² dx. · *plain:*
     the potential energy per unit area.
  5. *did:* Use the mean square · *tex:* $E_p=\tfrac12\rho g\overline{\eta^2}$ (7.41) · *why:* The mean square over a wavelength is by definition $\overline{\eta^2}=\frac1\lambda\int_0^\lambda\eta^2dx$; we use it to compare with (7.39). · *plain:* PE in terms of the mean-square height.
  6. *did:* Compare with D12 · *tex:* $E_p=E_k$ · *why:* (7.39) and (7.41) are the same expression: equipartition, as for
     any small oscillation of a conservative system without rotation. · *plain:* motion and height share the energy equally.
  7. *did:* Add and use a sinusoid · *tex:* $E=E_k+E_p=\rho g\overline{\eta^2}=\tfrac12\rho ga^2$ (7.42) · *why:* For η = a cos(kx −
     ωt), $\overline{\eta^2}$ = a²/2 (⟨cos²⟩ = ½, P151). · *plain:* total wave energy per square metre.
- **Result.** $E_p=\tfrac12\rho g\overline{\eta^2}=E_k$ (7.41), $E=\tfrac12\rho ga^2$ (7.42).
- **Check.** a = 1 m → E = 4905 J/m², half each (C06 code with `quad`) ✓. Units J/m² ✓.
- **What it means.** Wave energy grows with the square of the height: a 2 m sea holds four times the energy of a 1 m sea.
  Equipartition fails once the Earth's rotation matters (Ch. 13).
- **Traps.** Forgetting that the flat-state energy must be subtracted; $\overline{\eta^2}=a^2/2$ only for a sinusoid.

### D14 · The energy flux: (7.43) → (7.44) — ★★, 10 steps, in C06 (notebook)
- **Goal.** Find how much energy per second crosses a vertical line under the wave, and show it is energy × a speed that
  is not c.
- **Start.** $F=\big\langle\int_{-H}^0pu\,dz\big\rangle$ — *in words:* the water on one side pushes on the other with pressure p
  while moving with u.
- **Plan.** (1) Split the pressure. (2) Kill the background part. (3) Insert p′ and u. (4) Average and integrate. (5)
  Regroup as energy × speed.
- **Tools.** power of a force (P127) · mean of cos² (P151) · ∫cosh² from D12 · hyperbolic identities (P168) · (7.28).
- **Assumptions.** Linear wave; integrate to z = 0 (as D12).
- **Steps.**
  1. *did:* Write the pressure work across x = 0 · *tex:* $F=\big\langle\int_{-H}^0pu\,dz\big\rangle$ · *why:* Power of a force is
     force × velocity (P127): pressure p over the strip dz times u; ⟨⟩ averages over a period. · *plain:* energy per second
     through a vertical line, per metre of crest.
  2. *did:* Split the pressure · *tex:* $p=p'-\rho gz$ · *why:* (7.30) p′ ≡ p + ρgz: wave part plus still-water part. ·
     *plain:* two pieces.
  3. *did:* Separate the two parts · *tex:* $F=\big\langle\int p'u\,dz\big\rangle-\rho g\int_{-H}^0z\langle u\rangle dz$ · *why:*
     Linearity of the integral; z does not depend on t, so the time average moves inside onto u. · *plain:* the background
     pushes on a velocity that averages to zero.
  4. *did:* Use ⟨u⟩ = 0 · *tex:* $F=\big\langle\int_{-H}^0p'u\,dz\big\rangle$ (7.43) · *why:* u ∝ cos(kx − ωt) averages to zero over
     a period at every depth — which is why the book may pull ⟨u⟩ out of the integral. · *plain:* only the wave pressure
     does net work.
  5. *did:* Insert p′ and u · *tex:* $p'u=\frac{\rho a^2\omega^3}{k}\frac{\cosh^2k(z+H)}{\sinh^2kH}\cos^2(kx-\omega t)$ · *why:* First
     form of (7.31) p′ = ρ(aω²/k)(cosh/sinh)cos, and u from (7.27) (the book's 'u from (7.28)' means (7.27)). · *plain:*
     pressure and velocity are in step, so the product is never negative.
  6. *did:* Average in time · *tex:* $F=\frac{\rho a^2\omega^3}{2k\sinh^2kH}\int_{-H}^0\cosh^2k(z+H)dz$ · *why:* The time average of cos² over a period is ½ (P151); nothing else in the product depends on time. · *plain:* half the peak product.
  7. *did:* Integrate over the depth · *tex:* $F=\frac{\rho a^2\omega^3}{2k\sinh^2kH}\Big(\frac H2+\frac{\sinh2kH}{4k}\Big)$ · *why:* The
     cosh² integral of D12 step 5. · *plain:* the depth integral done.
  8. *did:* Factor out sinh 2kH/4k · *tex:* $F=\frac{\rho a^2\omega^3\sinh2kH}{8k^2\sinh^2kH}\Big(1+\frac{2kH}{\sinh2kH}\Big)$ · *why:*
     Write H/2 = (sinh 2kH/4k)(2kH/sinh 2kH), an identity, so that the depth term becomes a correction to 1. · *plain:* a depth factor appears.
  9. *did:* Simplify with (7.28) · *tex:* $\frac{\sinh2kH}{\sinh^2kH}=\frac{2}{\tanh kH},\ \omega^3=\omega\,gk\tanh kH\ \Rightarrow\
     F=\frac{\rho a^2g\omega}{4k}\Big(1+\frac{2kH}{\sinh2kH}\Big)$ · *why:* Double angle (P168); one ω² replaced by (7.28); tanh kH
     cancels. · *plain:* the flux in terms of g.
  10. *did:* Regroup as energy × speed · *tex:* $F=\big[\tfrac12\rho ga^2\big]\Big[\frac c2\Big(1+\frac{2kH}{\sinh2kH}\Big)\Big]$ (7.44) ·
      *why:* ρa²gω/(4k) = ½ρga² · ω/(2k) and ω/k = c. · *plain:* flux = energy density × a speed that is not c.
- **Result.** $F=\big[\tfrac12\rho ga^2\big]\big[\frac c2\big(1+\frac{2kH}{\sinh2kH}\big)\big]$ (7.44).
- **Check.** Units: J/m² × m/s = W/m ✓. Deep water: F = E c/2; shallow: F = E c. T = 10 s, a = 1 m deep: 38.3 kW/m, equal to
  `energy_flux(method="quad")` ✓.
- **What it means.** Energy moves at (c/2)(1 + 2kH/sinh 2kH), not at c: the group velocity of C09, (7.69). Only in
  shallow water do crests and energy travel together.
- **Traps.** Using the second form of (7.31) (it works but hides the ω³); forgetting ⟨u⟩ = 0 holds at every depth (so
  pulling it out is legal); losing the factor ½ of ⟨cos²⟩.

### D15 · The surface-tension condition (7.53)–(7.55) and capillary–gravity dispersion (7.56)–(7.57) — ★★, 10 steps, in C07 (notebook · `capillary_gravity_waves`)
- **Goal.** Add surface tension to the dynamic condition and find the new dispersion relation.
- **Start.** $\Delta p=\sigma(1/R_1+1/R_2)$ (1.5) — *in words:* a curved surface pushes harder on the side of its centre of
  curvature.
- **Plan.** (1) Keep one radius (straight crests). (2) Use the curvature of a graph. (3) Put the surface pressure into
  linear Bernoulli. (4) Redo D06 with the extra term.
- **Tools.** Laplace jump (R09) · curvature of a plane curve (P169) · gauge pressure (R05) · D04, D05, D06.
- **Assumptions.** Straight crests (∂/∂y = 0); small slope; air pressure uniform.
- **Steps.**
  1. *did:* Keep one radius of curvature · *tex:* $\frac1{R_2}=0,\quad\frac1{R_1}=\frac1R$ · *why:* Nothing varies along the
     crest (y), so the surface is straight in that direction: that radius is infinite. · *plain:* only the curvature in the
     x–z plane counts.
  2. *did:* Write the pressure jump · *tex:* $p_a-(p)_{z=\eta}=\sigma\frac1R$ · *why:* (1.5) between the air (p_a) and the water
     just below; 1/R > 0 when the surface bends up (η_xx > 0, centre in the air, higher pressure in the air). · *plain:*
     the air–water pressure difference is σ × curvature.
  3. *did:* Check the sign under a crest · *tex:* $\eta_{xx}<0\ \Rightarrow\ (p)_{z=\eta}>p_a$ · *why:* At a crest the centre
     of curvature is below, in the water, and (1.5) puts the higher pressure on that side. · *plain:* tension squeezes the
     crest downward. · *watch (E3):* "the rose arrow points into the water".
  4. *did:* Insert the curvature of a graph · *tex:* $\frac1R=\frac{\partial^2\eta/\partial x^2}{[1+(\partial\eta/\partial x)^2]^{3/2}}
     \cong\frac{\partial^2\eta}{\partial x^2}$ (7.53) · *why:* Curvature of z = η(x) (P169); for small slopes the denominator is
     1 + O((ka)²). · *plain:* for gentle waves curvature is η_xx.
  5. *did:* Use gauge pressure · *tex:* $(p)_{z=\eta}=-\sigma\frac{\partial^2\eta}{\partial x^2}$ (7.54) · *why:* Set p_a = 0 in
     step 2 with step 4. · *plain:* the water pressure just below the surface.
  6. *did:* Put it into linear Bernoulli · *tex:* $\big(\frac{\partial\phi}{\partial t}\big)_{z=0}=\frac\sigma\rho\frac{\partial^2\eta}
     {\partial x^2}-g\eta$ (7.55) · *why:* In D04 step 5 the pressure term is now −(σ/ρ)η_xx instead of zero; move it right
     and transfer to z = 0 as before. · *plain:* tension adds a curvature term to the dynamic condition.
  7. *did:* Keep the same potential · *tex:* $\phi=\frac{a\omega}{k}\frac{\cosh k(z+H)}{\sinh kH}\sin(kx-\omega t)$ (7.26) · *why:*
     Laplace, the bottom and the kinematic condition are unchanged, so D05 still holds; only D06 changes. · *plain:* the
     flow pattern is the same.
  8. *did:* Differentiate the sinusoid twice · *tex:* $\frac{\partial^2\eta}{\partial x^2}=-k^2a\cos(kx-\omega t)=-k^2\eta$ · *why:*
     η = a cos(kx − ωt) (7.2). · *plain:* for a sinusoid, curvature is −k² × height.
  9. *did:* Insert into (7.55) · *tex:* $\big(\frac{\partial\phi}{\partial t}\big)_{z=0}=-\Big(g+\frac{\sigma k^2}\rho\Big)\eta$ · *why:*
     Substitute η_xx = −k²η from step 8: (σ/ρ)(−k²η) − gη = −(g + σk²/ρ)η; both terms are now proportional to η. · *plain:* tension acts like extra gravity σk²/ρ.
  10. *did:* Redo D06 with g → g + σk²/ρ · *tex:* $\omega=\sqrt{k\Big(g+\frac{\sigma k^2}\rho\Big)\tanh kH}$ (7.56), $c=\sqrt{\Big(
      \frac gk+\frac{\sigma k}\rho\Big)\tanh kH}$ (7.57) · *why:* D06 used g only through (7.21); replacing it gives (7.56);
      dividing by k gives (7.57) (λ form with k = 2π/λ). · *plain:* capillary–gravity dispersion and speed. · *live (E3):*
      "ω = √(k(g + σk²/ρ)) = …".
- **Result.** (7.53)–(7.57).
- **Check.** Units: σk²/ρ [N/m × 1/m² × m³/kg = m/s²] ✓ — an acceleration, like g. σ = 0 returns (7.28) ✓. Number: a = 1 mm,
  λ = 1 cm: p under the crest = σak² = +28.7 Pa ✓.
- **What it means.** Short waves are restored mainly by tension, long ones by gravity; the crossover is near
  √(σ/ρg) ≈ 2.7 mm × 2π.
- **Traps.** The sign of the pressure under a crest (higher, not lower); using both radii for straight crests;
  forgetting that η_xx = −k²η only for a sinusoid.

### D16 · The minimum phase speed (7.58) — ★★, 7 steps, in C07 (notebook · `capillary_gravity_waves`)
- **Goal.** Find the slowest capillary–gravity wave and its wavelength.
- **Start.** $c=\sqrt{\big(\frac gk+\frac{\sigma k}\rho\big)\tanh kH}$ (7.57).
- **Plan.** (1) Deep water. (2) Minimise c² instead of c. (3) Slope zero, solve, check it is a minimum.
- **Tools.** slope zero at a minimum (gloss; P19) · deep limit (D08).
- **Assumptions.** Deep water, tanh kH ≈ 1 (ponds are deep for centimetre ripples).
- **Steps.**
  1. *did:* Take deep water · *tex:* $c^2=\frac gk+\frac{\sigma k}\rho$ · *why:* tanh kH ≈ 1 when H > 0.32λ (D08); here λ is
     centimetres. · *plain:* one term falls with k, the other grows.
  2. *did:* Minimise c² instead of c · *tex:* $\min c\iff\min c^2$ · *why:* c > 0 and squaring is increasing for positive
     numbers, so both have their minimum at the same k. · *plain:* easier algebra, same answer.
  3. *did:* Set the slope to zero · *tex:* $\frac{d(c^2)}{dk}=-\frac g{k^2}+\frac\sigma\rho=0$ · *why:* At a smooth minimum the
     derivative vanishes (P19). · *plain:* where gravity's decrease matches tension's increase.
  4. *did:* Solve for k · *tex:* $k_m=\sqrt{\frac{\rho g}\sigma},\qquad\frac g{k_m}=\frac{\sigma k_m}\rho$ · *why:* Multiply by k²:
     k² = ρg/σ; at this k the two terms of c² are equal. · *plain:* at the minimum both restoring effects are equal. · *set
     (E3):* λ = λ_m.
  5. *did:* Convert to wavelength · *tex:* $\lambda_m=\frac{2\pi}{k_m}=2\pi\sqrt{\frac\sigma{\rho g}}$ (7.58) · *why:* Wavelength and wavenumber are related by λ = 2π/k (7.1), so the slowest wavenumber gives the slowest wavelength. ·
     *plain:* the wavelength of the slowest wave.
  6. *did:* Evaluate c² there · *tex:* $c_{min}^2=\frac{2g}{k_m}=2\sqrt{\frac{g\sigma}\rho}$ · *why:* Both terms equal g/k_m, and
     g/k_m = g√(σ/ρg) = √(gσ/ρ). · *plain:* the minimum speed squared. · *live (E3):* "c_min² = 2√(gσ/ρ) = …".
  7. *did:* Take the root, check the minimum · *tex:* $c_{min}=\Big[\frac{4g\sigma}\rho\Big]^{1/4}$ (7.58), $\frac{d^2(c^2)}{dk^2}=
     \frac{2g}{k^3}>0$ · *why:* √(2√(gσ/ρ)) = (4gσ/ρ)^{1/4}; a positive second derivative means a minimum. · *plain:* no
     water wave travels slower than c_min.
- **Result.** $c_{min}=[4g\sigma/\rho]^{1/4}$ at $\lambda_m=2\pi\sqrt{\sigma/\rho g}$ (7.58).
- **Check.** Clean water at 20 °C (σ = 72.74 mN/m, ρ = 998.2): 23.1 cm/s at 1.71 cm (7.59) ✓; the log-grid scan +
  `minimize_scalar` agree to 1e-6 ✓. Units: (m/s² × N/m ÷ kg/m³)^{1/4} = m/s ✓.
- **What it means.** At λ_m the tangent to ω(k) passes through the origin, so c_g = c there; the group velocity has its own,
  lower minimum (≈ 17.8 cm/s) — the calm ring of C09.
- **Traps.** Minimising in λ instead of k gives the same point (both are monotone maps), but not the same derivative;
  forgetting deep water was assumed.

### D17 · The standing wave, its ψ (7.62) and u (7.63) — ★★, 8 steps, in C08 (notebook · `seiche_standing_waves`)
- **Goal.** Add a right-going and a left-going wave and find the surface, the streamlines and the velocity of the result.
- **Start.** $\eta=a\cos(kx-\omega t)$ (7.2), $\eta=a\cos[kx+\omega t]$ (7.61), $\psi=\frac{a\omega}k\frac{\sinh k(z+H)}{\sinh kH}
  \cos(kx-\omega t)$ (7.37).
- **Plan.** (1) Add the surfaces with sum-to-product. (2) Find the nodes. (3) Add the stream functions (with the right
  sign). (4) Differentiate for u.
- **Tools.** sum-to-product identities (P171) · superposition (N02) · ψ convention u = ∂ψ/∂z (N30).
- **Assumptions.** Linear problem; equal amplitudes and wavenumbers.
- **Steps.**
  1. *did:* Take the left-going wave · *tex:* $\eta_2=a\cos[kx+\omega t]$ (7.61) · *why:* Replacing ω by −ω in (7.2) reverses
     the direction; the dispersion relation contains only ω², so it is also a solution. · *plain:* an identical wave
     moving to −x. · *set (E4):* comp = left-going.
  2. *did:* Add the two surfaces · *tex:* $\eta=a\cos(kx-\omega t)+a\cos(kx+\omega t)$ · *why:* The problem is linear, so a sum of
     solutions is a solution (superposition). · *plain:* both waves at once.
  3. *did:* Use sum-to-product · *tex:* $\eta=2a\cos kx\cos\omega t$ · *why:* cos A + cos B = 2 cos((A + B)/2) cos((A − B)/2) with
     A = kx − ωt, B = kx + ωt; (A − B)/2 = −ωt and cos is even (P171). · *plain:* a fixed shape breathing in time. · *set
     (E4):* comp = both.
  4. *did:* Find the nodes · *tex:* $\eta=0\ \text{for all }t\ \text{at}\ kx=\pm\tfrac\pi2,\pm\tfrac{3\pi}2,\ldots$ · *why:* cos kx = 0
     there whatever cos ωt is. · *plain:* some points never move.
  5. *did:* Write the left-going stream function · *tex:* $\psi_2=-\frac{a\omega}k\frac{\sinh k(z+H)}{\sinh kH}\cos(kx+\omega t)$ ·
     *why:* (7.37) with ω → −ω: the prefactor aω/k flips sign — the reversed wave's water moves the other way. · *plain:*
     the reversed wave's ψ carries a minus sign.
  6. *did:* Add the stream functions · *tex:* $\psi=\frac{a\omega}k\frac{\sinh k(z+H)}{\sinh kH}[\cos(kx-\omega t)-\cos(kx+\omega t)]$
     · *why:* ψ is linear in the velocity, so it superposes too. · *plain:* the total ψ.
  7. *did:* Use the difference identity · *tex:* $\psi=\frac{2a\omega}k\frac{\sinh k(z+H)}{\sinh kH}\sin kx\sin\omega t$ (7.62) · *why:*
     cos A − cos B = −2 sin((A + B)/2) sin((A − B)/2) = −2 sin kx sin(−ωt) = 2 sin kx sin ωt (P171). · *plain:* the streamlines
     of a standing wave. · *live (E4):* "u amplitude 2aω cosh kH/sinh kH = …".
  8. *did:* Differentiate for u · *tex:* $u=\frac{\partial\psi}{\partial z}=2a\omega\frac{\cosh k(z+H)}{\sinh kH}\sin kx\sin\omega t$
     (7.63) · *why:* u = ∂ψ/∂z (our convention); d/dz sinh k(z + H) = k cosh k(z + H). · *plain:* horizontal velocity:
     largest at the nodes of η, zero at its crests and troughs.
- **Result.** $\eta=2a\cos kx\cos\omega t$, (7.62), (7.63).
- **Check.** `standing_wave_fields` equals the sum of two `wave_fields` (C08 code) ✓. Units: ψ [m²/s], u [m/s] ✓.
- **What it means.** Energy sloshes between potential (ωt = 0, η largest, u = 0) and kinetic (ωt = π/2, η = 0, u largest);
  vertical lines under crests and troughs are streamlines — walls can stand there (D18).
- **Traps.** Forgetting the minus sign of ψ₂ (it would give cos kx instead of sin kx); nodes of η are antinodes of u.

### D18 · Seiche wavelengths (7.64) and frequencies (7.65) — ★, 5 steps, in C08 (notebook · `seiche_standing_waves`)
- **Goal.** Find which standing waves fit a basin with vertical walls at x = 0 and x = L, and their frequencies.
- **Start.** $u=2a\omega\frac{\cosh k(z+H)}{\sinh kH}\sin kx\sin\omega t$ (7.63).
- **Plan.** (1) Satisfy the wall at x = 0. (2) Impose the wall at x = L. (3) Put the allowed k into (7.28).
- **Tools.** zeros of sine · dispersion relation (7.28).
- **Assumptions.** Vertical walls, uniform depth, no variation across the basin.
- **Steps.**
  1. *did:* Put a wall at x = 0 · *tex:* $u(0,z,t)\propto\sin0=0$ · *why:* A wall needs u = 0 at all depths and times; sin kx
     vanishes at x = 0 automatically — the book's choice of origin (an antinode of η at the wall). · *plain:* the first
     wall costs nothing.
  2. *did:* Put the second wall at x = L · *tex:* $\sin kL=0$ · *why:* u = 0 at x = L for all z and t requires the only
     x-factor, sin kx, to vanish there. · *plain:* the second wall selects k.
  3. *did:* Solve for the allowed k · *tex:* $kL=(n+1)\pi,\ n=0,1,2,\ldots\ \Rightarrow\ \lambda=\frac{2L}{n+1}$ (7.64) · *why:* The
     zeros of sine are multiples of π; k = 0 is no wave, so the count starts at π (n = 0). · *plain:* only 2L, L, 2L/3, …
     fit. · *set (E4):* n = 1.
  4. *did:* Insert into the dispersion relation · *tex:* $\omega^2=\frac{g(n+1)\pi}L\tanh\frac{(n+1)\pi H}L$ · *why:* Each
     travelling component obeys (7.28) ω² = gk tanh kH; put k = (n + 1)π/L. · *plain:* each allowed wavelength has its
     frequency.
  5. *did:* Take the root and the shallow form · *tex:* $\omega=\sqrt{\frac{\pi g(n+1)}L\tanh\frac{(n+1)\pi H}L}$ (7.65);
     $T\approx\frac{2L}{(n+1)\sqrt{gH}}$ · *why:* Positive root; for kH ≪ 1, tanh x ≈ x gives ω = (n + 1)π√(gH)/L. · *plain:*
     the natural frequencies of the basin. · *live (E4):* "T = …".
- **Result.** $\lambda=\frac{2L}{n+1}$ (7.64), $\omega=\sqrt{\frac{\pi g(n+1)}L\tanh\big[\frac{(n+1)\pi H}L\big]}$ (7.65).
- **Check.** Bathtub L = 1.5 m, H = 0.2 m: T₀ = 2.20 s; lake L = 50 km, H = 100 m: 53.2 min = 2L/√(gH) ✓.
- **What it means.** Boundaries pick discrete modes — an eigenvalue problem; seiches, harbour resonance and tidal basins
  (Ch. 13) are this.
- **Traps.** Starting n at 1 (the book counts from 0); using λ = L for the gravest mode (it is 2L); forgetting tanh for
  deep, short basins.

### D19 · Beats (7.66) and the group velocity $c_g=d\omega/dk$ (7.67) — ★, 7 steps, in C09 (notebook · `group_velocity_packets`)
- **Goal.** Show that two waves of nearly equal length make groups, and that the groups move at dω/dk.
- **Start.** $\eta=a\cos(k_1x-\omega_1t)+a\cos(k_2x-\omega_2t)$, $\omega_i=\omega(k_i)$.
- **Plan.** (1) Sum-to-product. (2) Read the speed of the fast factor and of the slow factor. (3) Let the waves merge.
- **Tools.** sum-to-product identities (P171) · D01's phase argument · derivative as a limit (P19).
- **Assumptions.** Equal amplitudes; Δk small (step 7).
- **Steps.**
  1. *did:* Add two nearby waves · *tex:* $\eta=a\cos(k_1x-\omega_1t)+a\cos(k_2x-\omega_2t)$ · *why:* Superposition; the medium
     fixes each ω_i through its dispersion relation. · *plain:* two waves of almost the same wavelength.
  2. *did:* Define means and differences · *tex:* $k=\tfrac{k_1+k_2}2,\ \omega=\tfrac{\omega_1+\omega_2}2,\ \Delta k=k_2-k_1,\
     \Delta\omega=\omega_2-\omega_1$ · *why:* These are the half-sums and half-differences that the sum-to-product identity (P171) needs; they also name the mean wave. · *plain:* a central wave and a
     small mismatch.
  3. *did:* Apply sum-to-product · *tex:* $\eta=2a\cos\big(\tfrac12\Delta k\,x-\tfrac12\Delta\omega\,t\big)\cos(kx-\omega t)$ (7.66) ·
     *why:* With A = k₁x − ω₁t, B = k₂x − ω₂t: (A + B)/2 = kx − ωt and (A − B)/2 = −(½Δk x − ½Δω t); cos is even. · *plain:* a
     fast carrier times a slow envelope.
  4. *did:* Correct the printed slip · *tex:* $\text{printed: }\cos\big(\tfrac12\Delta k\,x-\tfrac12\Delta\omega\,x\big)\ \to\ \text{read}\ \tfrac12\Delta\omega\,t$ ·
     *why:* The identity gives ½Δω t; the book's next line and its figure use t. With x the envelope would not depend on
     time and could never move. · *plain:* the envelope must move in time. · *set (E5):* show the printed ghost.
  5. *did:* Read the carrier's speed · *tex:* $kx-\omega t=\text{const}\Rightarrow\frac{dx}{dt}=\frac\omega k=c$ · *why:* D01's argument
     applied to the fast factor. · *plain:* crests move at the phase speed.
  6. *did:* Read the envelope's speed · *tex:* $\tfrac12\Delta k\,x-\tfrac12\Delta\omega\,t=\text{const}\Rightarrow\frac{dx}{dt}=\frac{\Delta
     \omega}{\Delta k}$ · *why:* The same argument for the slow factor; its wavelength is 4π/Δk and its period 4π/Δω. ·
     *plain:* nodes of the envelope move at Δω/Δk. · *live (E5):* "Δω/Δk = … vs dω/dk = …".
  7. *did:* Let the waves merge · *tex:* $c_g=\lim_{\Delta k\to0}\frac{\Delta\omega}{\Delta k}=\frac{d\omega}{dk}$ (7.67) · *why:* The
     difference quotient tends to the derivative (P19): the chord slope becomes the tangent slope of ω(k). · *plain:* groups
     move at the slope of the dispersion curve.
- **Result.** (7.66) and $c_g=\Delta\omega/\Delta k=d\omega/dk$ (7.67).
- **Check.** Deep water, k₁ = 0.9, k₂ = 1.1 rad/m: Δω/Δk = 1.568 m/s vs ½√(g/k) = 1.566 m/s at k = 1 ✓.
- **What it means.** No energy crosses a node of the envelope, so energy travels at the node speed — c_g (made precise
  for packets in D20 and for the flux in D21).
- **Traps.** Taking the printed Δω x literally; mixing the envelope's wavelength (4π/Δk) with the beat length between
  nodes (2π/Δk).

### D20 · A packet's envelope moves at c_g: $\eta=a(x-c_gt)\cos(kx-\omega t)$ (7.68) — ★★★, 12 steps, in C09 (notebook · `group_velocity_packets`)
- **Goal.** Show that a narrow packet keeps its envelope shape and carries it at c_g = dω/dk, and see what eventually
  spreads it. (The book cites Phillips (1977) without the steps; we write them out.)
- **Start.** $\eta(x,t)=\mathrm{Re}\int A(k)\,e^{i(kx-\omega(k)t)}dk$ with A(k) concentrated near k₀ — *in words:* a group is a
  continuous sum of sinusoids, each moving at its own speed.
- **Plan.** (1) Measure k from the carrier. (2) Taylor-expand ω(k). (3) Split the phase into carrier + envelope + small
  rest. (4) Drop the rest for short times and recognise the shifted envelope.
- **Tools.** Fourier integral and a packet's spectrum (P172) · Taylor expansion (P26) · substitution in an integral (P106) ·
  complex amplitudes (Re taken at the end).
- **Assumptions.** Narrow spectrum δk ≪ k₀ (step 2); short times t ≪ 1/(ω″δk²) (step 7); A symmetric about k₀ for a real
  envelope (step 10).
- **Steps.**
  1. *did:* Write the packet as a Fourier integral · *tex:* $\eta=\mathrm{Re}\int A(k)e^{i(kx-\omega(k)t)}dk$ · *why:* Linear waves
     superpose, and each Fourier component moves with its own ω(k) (N02, P172). · *plain:* a group is a sum of sinusoids.
  2. *did:* Say what 'narrow' means · *tex:* $A(k)\approx0\ \text{unless}\ \lvert k-k_0\rvert\lesssim\delta k,\ \delta k\ll k_0$ · *why:*
     A long, nearly sinusoidal group has a narrow spectrum (length ~ 1/δk, N66); this is the key assumption. · *plain:* only
     wavenumbers near k₀ are present.
  3. *did:* Measure k from the carrier · *tex:* $k=k_0+\kappa,\quad dk=d\kappa$ · *why:* A change of variable (P106); κ is small
     wherever A is not zero. · *plain:* count wavenumbers from the centre.
  4. *did:* Look at the initial shape · *tex:* $\eta(x,0)=\mathrm{Re}\big[e^{ik_0x}a(x)\big],\quad a(x)\equiv\int A(k_0+\kappa)e^{i\kappa x}d\kappa$ ·
     *why:* At t = 0 factor e^{ik₀x} out of the integral; what remains contains only small κ, so it varies slowly in x. ·
     *plain:* a carrier times a slowly varying envelope.
  5. *did:* Taylor-expand the dispersion relation · *tex:* $\omega(k_0+\kappa)=\omega_0+c_g\kappa+\tfrac12\omega''\kappa^2+\ldots$ · *why:*
     ω(k) is smooth and κ is small (P26); c_g = ω′(k₀), ω″ = ω″(k₀). · *plain:* near k₀ the curve is almost its tangent. ·
     *watch (E5):* "the purple tangent is the linear term".
  6. *did:* Substitute into the phase · *tex:* $kx-\omega t=(k_0x-\omega_0t)+\kappa(x-c_gt)-\tfrac12\omega''\kappa^2t$ · *why:* kx =
     k₀x + κx; subtract t × the Taylor series; group the terms by powers of κ. · *plain:* carrier phase + envelope phase +
     a small correction.
  7. *did:* Drop the quadratic term for short times · *tex:* $\tfrac12\lvert\omega''\rvert\kappa^2t\ll1\ \text{while}\ t\ll\frac{1}{\lvert
     \omega''\rvert\delta k^2}$ · *why:* κ is at most about δk, so this term changes the phase by much less than a radian for
     such times. · *plain:* for a while the curvature of ω(k) does not matter.
  8. *did:* Factor out the carrier · *tex:* $\eta\approx\mathrm{Re}\Big[e^{i(k_0x-\omega_0t)}\int A(k_0+\kappa)e^{i\kappa(x-c_gt)}d\kappa\Big]$ ·
     *why:* The carrier phase does not depend on κ, so it leaves the integral. · *plain:* the carrier moves out in front.
  9. *did:* Recognise the envelope · *tex:* $\int A(k_0+\kappa)e^{i\kappa(x-c_gt)}d\kappa=a(x-c_gt)$ · *why:* It is step 4's integral
     with x replaced by x − c_g t. · *plain:* the envelope is the initial envelope, shifted by c_g t. · *set (E5):* order = 1
     (no spreading).
  10. *did:* Take the real part · *tex:* $\eta=a(x-c_gt)\cos(k_0x-\omega_0t)$ (7.68) · *why:* For A symmetric about k₀ the
      envelope a is real, and Re[a e^{iθ}] = a cos θ. · *plain:* envelope at c_g, carrier at c = ω₀/k₀.
  11. *did:* Look at the dropped term · *tex:* $\text{later: width}\sim\lvert\omega''\rvert\,\delta k\,t$ · *why:* Components with
      different κ pick up different phases −½ω″κ²t and drift apart: dispersion spreads the packet; zero only if ω″ = 0
      (shallow water). · *plain:* the dropped term is what makes packets spread. · *set (E5):* order = 2, narrow packet.
  12. *did:* Identify the energy speed · *tex:* $E(x,t)=\tfrac12\rho g\,a(x-c_gt)^2$ · *why:* The local energy density (7.42) is
      carried by the envelope, so it moves at c_g — the same speed D21 finds for the flux, (7.71) F = E c_g. · *plain:* the
      energy travels with the group. · *live (E5):* "envelope peak at c_g t = … m".
- **Result.** $\eta=a(x-c_gt)\cos(kx-\omega t)$ (7.68), $c_g=d\omega/dk$, valid for $t\ll1/(\lvert\omega''\rvert\delta k^2)$.
- **Check.** Units: c_g t [m] ✓. Shallow water: ω″ = 0, the result is exact for all times ✓. Numbers: the FFT packet of C09
  (k₀ = 1 rad/m, deep) moves its envelope peak at 1.566 m/s within 1 % ✓.
- **sympy check intent (★★★).** Re-run steps 4–10 for a Gaussian spectrum and an exactly quadratic ω: (a) with ω″ = 0 the
  envelope integral equals a(x − c_g t) identically; (b) with ω″ ≠ 0 the squared modulus of the envelope still peaks at
  ξ = x − c_g t = 0 for every t (the group moves at c_g) while its width grows. `check_src` sketch:
  ```python
  import sympy as sp                                    # symbolic algebra
  kap, xi, x = sp.symbols('kappa xi x', real=True)       # wavenumber offset, moving coordinate x - c_g t, position
  t, s, w2 = sp.symbols('t s omega2', positive=True)     # time, spectral width δk, curvature ω″ (> 0 for the check)
  A = sp.exp(-kap**2/(2*s**2))                           # a narrow Gaussian spectrum centred on k0 (step 2)
  a0 = sp.integrate(A*sp.exp(sp.I*kap*x), (kap, -sp.oo, sp.oo))    # the initial envelope a(x), step 4
  env1 = sp.integrate(A*sp.exp(sp.I*kap*xi), (kap, -sp.oo, sp.oo)) # step 9 with the quadratic term dropped
  print(sp.simplify(env1 - a0.subs(x, xi)))              # 0: the envelope is a(x - c_g t) exactly (steps 9-10)
  env2 = sp.integrate(A*sp.exp(sp.I*kap*xi - sp.I*w2*kap**2*t/2), (kap, -sp.oo, sp.oo), conds='none')  # keep ω″
  m2 = sp.simplify(sp.expand_complex(env2*sp.conjugate(env2)).replace(sp.exp_polar, sp.exp))   # |envelope|², real
  print(sp.simplify(sp.diff(m2, xi).subs(xi, 0)))        # 0: the peak stays at xi = 0, i.e. at x = c_g t (step 12)
  print(sp.simplify(m2.subs(xi, 0)/m2.subs({xi: 0, t: 0})))   # 1/sqrt(1 + (omega2 s² t)²): spreading sets in at t ~ 1/(ω″δk²), step 7
  ```
- **What it means.** Groups, energy and any signal travel at dω/dk; crests only at ω/k. Swell forecasting, tsunami
  arrival and every "where does the energy go" question of Ch. 13 use this.
- **Traps.** Thinking the carrier moves at c_g (it moves at c); dropping the quadratic term forever (packets do spread);
  forgetting that a narrow spectrum is essential (a single spike has no group).

### D21 · $c_g=\frac c2\big[1+\frac{2kH}{\sinh2kH}\big]$ (7.69), its limits (7.70) and $F=Ec_g$ (7.71) — ★★, 11 steps, in C09 (notebook · `group_velocity_packets`)
- **Goal.** Differentiate the dispersion relation of water waves, read off its limits, and recognise the second factor of
  the energy flux (7.44).
- **Start.** $\omega^2=gk\tanh kH$ (7.28).
- **Plan.** (1) Differentiate ω² implicitly. (2) Divide by c and simplify. (3) Take the limits. (4) Compare with (7.44).
- **Tools.** chain and product rules (P49, P38) · hyperbolic identities (P168) · D14.
- **Assumptions.** Linear gravity waves, constant depth.
- **Steps.**
  1. *did:* Differentiate ω² · *tex:* $2\omega\frac{d\omega}{dk}=\frac{d}{dk}(gk\tanh kH)$ · *why:* Differentiating ω² with the
     chain rule (P49) avoids differentiating a square root. · *plain:* the slope of ω² in terms of the slope of ω.
  2. *did:* Use the product rule · *tex:* $\frac{d}{dk}(gk\tanh kH)=g\tanh kH+gkH\,\mathrm{sech}^2kH$ · *why:* Product of k and
     tanh kH (P38); d/dk tanh kH = H sech² kH (chain rule, tanh′ = sech²). · *plain:* two contributions.
  3. *did:* Solve for dω/dk · *tex:* $c_g=\frac{g\tanh kH+gkH\,\mathrm{sech}^2kH}{2\omega}$ · *why:* Divide both sides by 2ω, which is not zero for a wave; the left side is then dω/dk alone. · *plain:*
     the group velocity.
  4. *did:* Divide by c = ω/k · *tex:* $\frac{c_g}c=\frac{k\,(g\tanh kH+gkH\,\mathrm{sech}^2kH)}{2\omega^2}$ · *why:* Since c = ω/k, dividing by c is multiplying by k/ω; the ratio c_g/c is easier to simplify.
     · *plain:* the ratio to the phase speed.
  5. *did:* Replace ω² by (7.28) · *tex:* $\frac{c_g}c=\frac{gk\tanh kH+gk\cdot kH\,\mathrm{sech}^2kH}{2gk\tanh kH}$ · *why:* Replace ω² in the denominator by the dispersion relation (7.28) ω² = gk tanh kH, so g and k can cancel. · *plain:* g and k are about to cancel.
  6. *did:* Divide through · *tex:* $\frac{c_g}c=\frac12\Big(1+\frac{kH\,\mathrm{sech}^2kH}{\tanh kH}\Big)$ · *why:* Each term divided
     by gk tanh kH. · *plain:* one half, plus a depth correction.
  7. *did:* Simplify sech²/tanh · *tex:* $c_g=\frac c2\Big[1+\frac{2kH}{\sinh(2kH)}\Big]$ (7.69) · *why:* sech²x/tanh x =
     1/(sinh x cosh x) = 2/sinh 2x (P168). · *plain:* the group velocity of water waves. · *live (E5):* "c_g/c = …".
  8. *did:* Take deep water · *tex:* $kH\to\infty:\ \frac{2kH}{\sinh2kH}\to0\ \Rightarrow\ c_g=\frac c2$ (7.70) · *why:* sinh grows
     exponentially, faster than its argument. · *plain:* energy moves at half the crest speed.
  9. *did:* Take shallow water · *tex:* $kH\to0:\ \frac{2kH}{\sinh2kH}\to1\ \Rightarrow\ c_g=c$ (7.70) · *why:* sinh y ≈ y for small
     y (P168). · *plain:* groups and crests together — no dispersion.
  10. *did:* Compare with the energy flux · *tex:* $F=\big[\tfrac12\rho ga^2\big]\Big[\frac c2\Big(1+\frac{2kH}{\sinh2kH}\Big)\Big]=Ec_g$
      (7.71) · *why:* (7.44)'s first factor is E (7.42) and its second is exactly (7.69). · *plain:* energy flux = energy ×
      group velocity. · *live (E5):* "F = E c_g = … W/m".
  11. *did:* Read it as an energy speed · *tex:* $c_E\equiv\frac FE=c_g$ · *why:* A flux divided by the density of what flows is
      the speed at which it flows (as mass flux ρu over ρ is u). · *plain:* wave energy travels at the group velocity.
- **Result.** $c_g=\frac c2\big[1+\frac{2kH}{\sinh(2kH)}\big]$ (7.69), $c_g=c/2$ (deep), $c_g=c$ (shallow) (7.70),
  $F=Ec_g$ (7.71).
- **Check.** H = 10 m, λ = 50 m: c = 8.15 m/s, c_g = 5.74 m/s (ratio 0.705) ✓ (numeric central difference agrees to 1e-8).
  Pure capillary ω² = σk³/ρ gives c_g = 3c/2 (Exercise 7.9, `group_velocity(g=0)`) ✓.
- **What it means.** C06's cliff-hanger is solved: the speed in (7.44) is the group velocity. Swell energy from a distant
  storm arrives at c_g.
- **Traps.** Forgetting the product rule (tanh kH depends on k); deep-water c_g = c/2, not 2c; confusing sinh(2kH) with
  2 sinh kH.

### D22 · Crest conservation (7.74) and k carried at c_g (7.75) — ★★, 8 steps, in C10 (notebook · `wave_rays_refraction`)
- **Goal.** Show that crests are conserved and that the local wavenumber is carried at the group velocity.
- **Start.** $\eta=a(x,t)\cos[\theta(x,t)]$ (7.72) — *in words:* a wave train whose amplitude and wavelength change slowly.
- **Plan.** (1) Define local k and ω from the phase. (2) Cross-differentiate. (3) Use the local dispersion relation. (4)
  Read the result along characteristics.
- **Tools.** mixed partials commute (P121) · chain rule (P91) · first-order wave equation and characteristics (P174).
- **Assumptions.** Slowly varying train (each piece obeys the uniform-medium ω(k)) from step 5; homogeneous medium.
- **Steps.**
  1. *did:* Write a slowly varying train · *tex:* $\eta=a(x,t)\cos[\theta(x,t)]$ (7.72) · *why:* Amplitude and wavelength change
     little over one wavelength and one period; a uniform train has θ = kx − ωt. · *plain:* a wave whose properties drift
     slowly.
  2. *did:* Define local k and ω · *tex:* $k\equiv\frac{\partial\theta}{\partial x},\qquad\omega\equiv-\frac{\partial\theta}{\partial t}$ (7.73)
     · *why:* For θ = kx − ωt these return the usual k and ω; in general they are the local rates of change of phase. ·
     *plain:* crests per metre and crests per second, measured locally.
  3. *did:* Cross-differentiate · *tex:* $\frac{\partial k}{\partial t}=\frac{\partial^2\theta}{\partial t\,\partial x}=\frac{\partial^2\theta}{\partial x\,
     \partial t}=-\frac{\partial\omega}{\partial x}$ · *why:* θ is smooth, so mixed partial derivatives commute (Schwarz, P121). ·
     *plain:* the time change of k is minus the space change of ω.
  4. *did:* Write it as a conservation law · *tex:* $\frac{\partial k}{\partial t}+\frac{\partial\omega}{\partial x}=0$ (7.74) · *why:*
     It has the form of a continuity equation with density k (crests per metre) and flux ω (crests per second). · *plain:*
     crests are neither created nor destroyed.
  5. *did:* Use the local dispersion relation · *tex:* $\omega(x,t)=\omega\big(k(x,t)\big)$ · *why:* In a homogeneous medium each
     short stretch of a slowly varying train obeys the same ω(k) as a uniform wave. · *plain:* frequency follows the local
     wavenumber. · *watch (E6):* "crest lines end at the group's front".
  6. *did:* Apply the chain rule · *tex:* $\frac{\partial\omega}{\partial x}=\frac{d\omega}{dk}\frac{\partial k}{\partial x}=c_g\frac{\partial k}{\partial x}$
     · *why:* ω depends on x only through k (P91); dω/dk = c_g (7.67). · *plain:* ω changes in space because k does.
  7. *did:* Substitute · *tex:* $\frac{\partial k}{\partial t}+c_g\frac{\partial k}{\partial x}=0$ (7.75) · *why:* Substitute step 6 into crest conservation (7.74) ∂k/∂t + ∂ω/∂x = 0; only k is left as unknown.
     · *plain:* k obeys a first-order wave equation with speed c_g.
  8. *did:* Read along characteristics · *tex:* $\frac{dk}{dt}=0\ \text{along}\ \frac{dx}{dt}=c_g$ · *why:* P174: ∂q/∂t + c ∂q/∂x = 0
     keeps q constant along dx/dt = c; here q = k and c = c_g, so ω(k) is constant too. · *plain:* an observer moving at c_g
     always sees the same wavelength. · *live (E6):* "c_g = … m/s is the slope of the thick lines".
- **Result.** $\partial k/\partial t+\partial\omega/\partial x=0$ (7.74), $\partial k/\partial t+c_g\,\partial k/\partial x=0$ (7.75).
- **Check.** A chirped train θ = 0.5x + 0.002x² − 2t satisfies (7.74) to 1e-9 (C10 code) ✓. The FFT packet's x–t diagram
  shows k constant along lines of slope c_g ✓.
- **What it means.** The group velocity is the speed at which wavenumbers (and, as D23 shows, frequencies) are advected —
  the kinematic basis of ray tracing.
- **Traps.** k is constant along dx/dt = c_g, not along dx/dt = c; the local dispersion relation is an assumption about
  slow variation.

### D23 · Frequency is constant along rays: (7.76) → (7.79) — ★★, 8 steps, in C10 (notebook · `wave_rays_refraction`)
- **Goal.** Let the depth change slowly with x and show that ω — not k — is carried unchanged at c_g.
- **Start.** $\frac{\partial k}{\partial t}+\frac{\partial\omega}{\partial x}=0$ (7.74) and $\omega=\sqrt{gk\tanh[kH(x)]}$.
- **Plan.** (1) Write ω(k, x). (2) Differentiate ω in time at fixed x. (3) Multiply (7.74) by c_g and substitute. (4) Add
  the companion equation for k.
- **Tools.** chain rule with a variable held fixed (P39, P91) · characteristics (P174).
- **Assumptions.** Steady medium (H does not change in time); slowly varying depth.
- **Steps.**
  1. *did:* Let the depth vary with x · *tex:* $\omega=\omega(k,x)$ (7.76) · *why:* With H = H(x), the local dispersion relation
     ω = √(gk tanh kH(x)) depends on x directly as well as through k. · *plain:* frequency depends on the wavenumber and on
     where you are.
  2. *did:* Define the local group velocity · *tex:* $\frac{\partial\omega(k,x)}{\partial k}=c_g$ (7.77) · *why:* The slope of ω
     against k at a fixed position (P39: derivative with x held fixed). · *plain:* the group velocity of the local medium.
  3. *did:* Differentiate ω in time at fixed x · *tex:* $\frac{\partial\omega}{\partial t}=\Big(\frac{\partial\omega}{\partial k}\Big)_x
     \frac{\partial k}{\partial t}$ · *why:* Chain rule; at fixed x the medium does not change in time, so ω changes only
     through k. · *plain:* the local frequency changes only if the local wavenumber does.
  4. *did:* Recognise c_g · *tex:* $c_g\frac{\partial k}{\partial t}=\frac{\partial\omega}{\partial t}$ (7.78) · *why:* At fixed x the slope (∂ω/∂k)_x is the local group velocity c_g of (7.77), so step 3 reads c_g ∂k/∂t. · *plain:* the book's (7.78).
  5. *did:* Multiply crest conservation by c_g · *tex:* $c_g\frac{\partial k}{\partial t}+c_g\frac{\partial\omega}{\partial x}=0$ · *why:*
     (7.74) times c_g — not c — to create the combination of step 4. · *plain:* prepare to substitute.
  6. *did:* Substitute (7.78) · *tex:* $\frac{\partial\omega}{\partial t}+c_g\frac{\partial\omega}{\partial x}=0$ (7.79) · *why:* Use (7.78) c_g ∂k/∂t = ∂ω/∂t to replace the first term; the equation then contains ω only. · *plain:* ω obeys a first-order wave equation with speed c_g. · *watch (E6):* "the black ω line
     stays flat".
  7. *did:* Read along rays · *tex:* $\frac{d\omega}{dt}=0\ \text{along}\ \frac{dx}{dt}=c_g;\qquad\frac{\partial\omega}{\partial t}+\mathbf c_g\cdot
     \nabla\omega=0$ · *why:* Characteristics (P174); the book states the 3-D form. · *plain:* frequency is conserved along
     rays while k, c and c_g change.
  8. *did:* Add the companion equation (ours) · *tex:* $\frac{dk}{dt}=-\Big(\frac{\partial\omega}{\partial x}\Big)_k\ \text{along}\ \frac{dx}{dt}=c_g$
     · *why:* In (7.74) write ∂ω/∂x = c_g∂k/∂x + (∂ω/∂x)_k; the first part joins ∂k/∂t. k changes where the medium changes —
     this bends rays (D24). · *plain:* the wavenumber changes along a ray as the depth changes. · *live (E6):* "ω drift …".
- **Result.** $\frac{\partial\omega}{\partial t}+c_g\frac{\partial\omega}{\partial x}=0$ (7.79), with the ray equations $\frac{dx}{dt}=
  \frac{\partial\omega}{\partial k}$, $\frac{dk}{dt}=-\frac{\partial\omega}{\partial x}$.
- **Check.** `ray_trace` over a 1:50 beach conserves ω to < 1e-8 relative while k grows from 0.071 to 0.181 rad/m ✓.
- **What it means.** In a steady medium the period of a wave never changes on its way; wavelength and direction do. This
  is the ray (WKB) theory Ch. 13 uses for internal waves in N(z) and for Rossby waves.
- **Traps.** Multiplying (7.74) by c instead of c_g; assuming k is conserved in an inhomogeneous medium; forgetting the
  medium must be steady (a changing tide breaks ω-conservation).

### D24 · Refraction on a plane beach: $k\sin\alpha=\text{const}$ (Snell) — ★★, 9 steps, in C10 (notebook · `wave_rays_refraction`)
- **Goal.** Show why crests arriving at an angle turn parallel to the shore. (The book gives the argument in words only;
  we write it out.)
- **Start.** The ray equations of D23 in two dimensions, with depth H(x) on a plane beach.
- **Plan.** (1) Write 2-D rays. (2) Use straight contours: k_y is conserved. (3) ω is conserved: |k| grows as H falls. (4)
  Combine into Snell's law.
- **Tools.** ray equations (D23 step 8) · Snell's law for waves (P175) · dispersion relation (7.28).
- **Assumptions.** Straight, parallel depth contours (H = H(x)); slowly varying depth; steady medium.
- **Steps.**
  1. *did:* Write the 2-D ray equations · *tex:* $\frac{d\mathbf x}{dt}=\nabla_{\mathbf k}\omega,\qquad\frac{d\mathbf k}{dt}=-\nabla_{\mathbf x}\omega$ ·
     *why:* D23 steps 7–8 in two dimensions with ω = ω(|k|, H(x, y)) (our extension of (7.79): Hamilton's equations for rays).
     · *plain:* a ray moves at c_g, and its wavenumber changes where the medium changes.
  2. *did:* Use straight depth contours · *tex:* $H=H(x)\ \Rightarrow\ \frac{\partial\omega}{\partial y}=0$ · *why:* On a plane beach the
     depth changes only with distance x from the shoreline. · *plain:* nothing changes along the shore.
  3. *did:* Conclude k_y is conserved · *tex:* $\frac{dk_y}{dt}=0\ \Rightarrow\ k_y=\text{const}$ · *why:* The y-component of the
     second ray equation with ∂ω/∂y = 0. · *plain:* the along-shore wavenumber never changes.
  4. *did:* Recall ω is conserved · *tex:* $\omega=\text{const along the ray}$ · *why:* (7.79): the medium is steady (D23). ·
     *plain:* the period stays the same all the way in.
  5. *did:* Get |k| from the local depth · *tex:* $\omega^2=g\lvert\mathbf k\rvert\tanh(\lvert\mathbf k\rvert H(x))$ · *why:* (7.28) with ω
     fixed: as H falls, tanh(|k|H) falls, so |k| must rise to keep ω (C03's `wavenumber_from_omega`). · *plain:* the wave
     shortens as the water shallows. · *watch (E6):* "the teal α curve falls".
  6. *did:* Define the angle to the shore normal · *tex:* $k_y=\lvert\mathbf k\rvert\sin\alpha$ · *why:* α is measured between k and
     the x-axis, the normal to the shoreline and the contours. · *plain:* α = 0 means crests parallel to the shore.
  7. *did:* Combine · *tex:* $\lvert\mathbf k\rvert\sin\alpha=\lvert\mathbf k_0\rvert\sin\alpha_0$ · *why:* Steps 3 and 6: the product is
     the conserved k_y (P175). · *plain:* Snell's law for water waves.
  8. *did:* Solve for α · *tex:* $\sin\alpha=\frac{\lvert\mathbf k_0\rvert\sin\alpha_0}{\lvert\mathbf k(x)\rvert}$ · *why:* Divide by |k| >
     0; |k| grows toward the shore (step 5), so sin α shrinks. · *plain:* the ray turns toward the shore normal. · *live
     (E6):* "sin α = … × 0.5/… = …".
  9. *did:* Take the shallow limit · *tex:* $\lvert\mathbf k\rvert\approx\frac{\omega}{\sqrt{gH}}\to\infty\ \text{as}\ H\to0\ \Rightarrow\ \alpha\to0$ ·
     *why:* (7.49): c = ω/|k| = √(gH) → 0 at the shoreline. · *plain:* crests end up parallel to the depth contours.
- **Result.** $\lvert\mathbf k\rvert\sin\alpha=\text{const}$ along a ray on a plane beach; α → 0 at the shore.
- **Check.** T = 8 s, 30° at 20 m → 11.3° at 2 m ✓ (worked example); the RK4 ray and the closed form agree to 1e-6 ✓.
- **What it means.** Refraction is why breakers line up with beaches and why energy focuses on headlands and over
  ridges; the same law bends light (optics) and sound.
- **Traps.** Conserving k instead of ω along the ray; measuring α from the shoreline instead of its normal; forgetting
  that the law needs straight contours (an island needs the full ray equations).

### D25 · The Bélanger relation: (7.80) → (7.81) — ★★, 10 steps, in C11 (notebook · `hydraulic_jump`)
- **Goal.** Find the downstream depth of a hydraulic jump from the upstream depth and speed, using only mass and momentum.
- **Start.** Mass $Q=u_1H_1=u_2H_2$ and the control-volume momentum balance (4.17) — *in words:* what flows in flows out, and
  momentum out − momentum in = net force.
- **Plan.** (1) Choose the CV and its face forces. (2) Balance momentum. (3) Eliminate the velocities and cancel the
  common factor. (4) Solve the quadratic and keep the physical root.
- **Tools.** CV momentum (Ch. 4 §4.4; momentum flux P114) · hydrostatic pressure (Ch. 1 §1.7) · Froude number (R10) ·
  quadratic formula and Vieta (P71).
- **Assumptions.** Steady jump; horizontal, frictionless bed; uniform velocity and hydrostatic pressure on the two faces
  (far from the roller); atmospheric pressure cancels (gauge).
- **Steps.**
  1. *did:* Draw the control volume · *tex:* $\text{face 1: }(H_1,u_1),\quad\text{face 2: }(H_2,u_2)$ · *why:* The jump is steady in
     this frame; vertical faces where the flow is uniform and hydrostatic keep the bookkeeping exact; the roller stays
     inside. · *plain:* a box round the jump. · *set (E7):* the CV faces light up.
  2. *did:* Write mass conservation · *tex:* $Q=u_1H_1=u_2H_2$ · *why:* Steady and incompressible: the volume flow per unit width
     entering must leave. · *plain:* the same flow rate crosses both faces.
  3. *did:* Find the face pressure forces · *tex:* $\int_0^{H}\rho g(H-z)\,dz=\tfrac12\rho gH^2$ · *why:* Hydrostatic gauge pressure
     ρg(H − z) on a vertical face of height H, integrated over the depth (per unit width). · *plain:* each face is pushed by
     the weight of its column. · *watch (E7):* "½ρgH² on each face".
  4. *did:* Balance momentum · *tex:* $\rho Q(u_2-u_1)=\tfrac12\rho gH_1^2-\tfrac12\rho gH_2^2$ · *why:* (4.17) with d/dt = 0 and no
     x-body force: momentum outflow − inflow = net pressure force; the frictionless bed gives no x-force. · *plain:* the
     momentum lost equals the difference of the pushes.
  5. *did:* Eliminate the velocities · *tex:* $Q^2\Big(\frac1{H_2}-\frac1{H_1}\Big)=\tfrac12g(H_1^2-H_2^2)$ (7.80) · *why:* u₁ = Q/H₁,
     u₂ = Q/H₂; divide by ρ. · *plain:* the jump relation in depths.
  6. *did:* Factor both sides · *tex:* $\frac{Q^2(H_1-H_2)}{H_1H_2}=\tfrac12g(H_1-H_2)(H_1+H_2)$ · *why:* 1/H₂ − 1/H₁ = (H₁ − H₂)/(H₁H₂)
     and a² − b² = (a − b)(a + b). · *plain:* both sides share H₁ − H₂.
  7. *did:* Cancel the common factor · *tex:* $Q^2=\tfrac12gH_1H_2(H_1+H_2)$ · *why:* H₁ ≠ H₂ for a jump (H₁ = H₂ is the trivial
     'no jump' solution). · *plain:* the flow rate links both depths.
  8. *did:* Make it dimensionless · *tex:* $2\mathrm{Fr}_1^2=r(1+r),\quad r=\frac{H_2}{H_1},\quad\mathrm{Fr}_1^2=\frac{Q^2}{gH_1^3}=\frac{u_1^2}{gH_1}$
     · *why:* Divide by ½gH₁³; (4.104) Fr = u/√(gH). · *plain:* one quadratic in the depth ratio.
  9. *did:* Use the quadratic formula · *tex:* $r^2+r-2\mathrm{Fr}_1^2=0\ \Rightarrow\ r=\tfrac12\big(-1\pm\sqrt{1+8\mathrm{Fr}_1^2}\big)$ ·
     *why:* Roots of r² + r + c = 0; their product is −2Fr₁² < 0 (Vieta, P71), so exactly one root is positive. · *plain:* two
     candidates, one negative. · *live (E7):* "r = ½(−1 + √(1 + 8 × …)) = …".
  10. *did:* Keep the positive root · *tex:* $\frac{H_2}{H_1}=\tfrac12\big(-1+\sqrt{1+8\mathrm{Fr}_1^2}\big)$ (7.81) · *why:* A depth ratio
      must be positive. · *plain:* the conjugate depth of the jump.
- **Result.** $\frac{H_2}{H_1}=\frac12\big(-1+\sqrt{1+8\mathrm{Fr}_1^2}\big)$ (7.81).
- **Check.** Fr₁ = 1 → r = 1 (no jump) ✓. H₁ = 0.1 m, Fr₁ = 3 → H₂ = 0.377 m; momentum terms 882.9 − 234.1 = 697.9 − 49.05 N/m
  ✓. Matches the published Bélanger form (V1 cross-check) and `np.roots` from scratch ✓.
- **What it means.** Momentum alone fixes the jump height; energy is not conserved and was never used — it decides the
  direction (D26).
- **Traps.** Using Bernoulli across the jump (energy is lost); forgetting the pressure forces are ½ρgH², not ρgH²;
  keeping the negative root.

### D26 · The energy loss of a jump and the second law — ★★, 8 steps, in C11 (notebook · `hydraulic_jump`)
- **Goal.** Compute the change of mechanical energy across a jump and show why only upward jumps exist.
- **Start.** $E=\frac{u^2}2+gH=\frac{Q^2}{2H^2}+gH$ and $Q^2=\tfrac12gH_1H_2(H_1+H_2)$ (D25 step 7).
- **Plan.** (1) Subtract the two faces. (2) Insert Q². (3) Put everything over one denominator and factor. (4) Read the
  sign.
- **Tools.** difference of squares · factorising · D25.
- **Assumptions.** As D25; E is per unit mass of a surface particle.
- **Steps.**
  1. *did:* Write the energy per unit mass · *tex:* $E=\frac{Q^2}{2H^2}+gH$ · *why:* Bernoulli's sum for a surface particle: gauge
     pressure 0, height H above the bed, speed u = Q/H. · *plain:* kinetic plus potential energy of surface water.
  2. *did:* Subtract the two faces · *tex:* $E_2-E_1=\frac{Q^2}2\Big(\frac1{H_2^2}-\frac1{H_1^2}\Big)+g(H_2-H_1)$ · *why:* Evaluate step 1
     at both faces. · *plain:* the change across the jump.
  3. *did:* Combine the fractions · *tex:* $\frac1{H_2^2}-\frac1{H_1^2}=\frac{(H_1-H_2)(H_1+H_2)}{H_1^2H_2^2}$ · *why:* Common
     denominator and a difference of squares. · *plain:* ready for Q².
  4. *did:* Insert Q² from D25 · *tex:* $\frac{Q^2}2\cdot\frac{(H_1-H_2)(H_1+H_2)}{H_1^2H_2^2}=\frac{g(H_1+H_2)^2(H_1-H_2)}{4H_1H_2}$ ·
     *why:* Q² = ½gH₁H₂(H₁ + H₂); one factor H₁H₂ cancels. · *plain:* the kinetic change in depths only.
  5. *did:* Put both terms over 4H₁H₂ · *tex:* $E_2-E_1=\frac{g(H_2-H_1)}{4H_1H_2}\big[4H_1H_2-(H_1+H_2)^2\big]$ · *why:* Write g(H₂ − H₁)
     = g(H₂ − H₁)·4H₁H₂/(4H₁H₂) and H₁ − H₂ = −(H₂ − H₁). · *plain:* factor out the common pieces.
  6. *did:* Simplify the bracket · *tex:* $4H_1H_2-(H_1+H_2)^2=-(H_2-H_1)^2$ · *why:* Expand: 4H₁H₂ − H₁² − 2H₁H₂ − H₂² = −(H₁² −
     2H₁H₂ + H₂²). · *plain:* a negative square.
  7. *did:* Collect · *tex:* $E_2-E_1=-\frac{g(H_2-H_1)^3}{4H_1H_2}$ · *why:* Multiply step 5 by step 6. (The extracted text drops the
     minus sign; the page has it.) · *plain:* the energy change is a cube of the depth change. · *set (E7):* Fr₁ = 0.7.
  8. *did:* Read the sign · *tex:* $H_2>H_1\Rightarrow E_2<E_1;\qquad H_2<H_1\Rightarrow E_2>E_1\ (\text{forbidden})$ · *why:* The
     turbulent roller can only dissipate mechanical energy (second law); a 'jump down' would create energy, so only Fr₁ > 1
     (H₂ > H₁) happens. Head loss = (E₁ − E₂)/g. · *plain:* water can jump up, never down. · *live (E7):* "E₂ − E₁ = … J/kg".
- **Result.** $E_2-E_1=-(H_2-H_1)\frac{g(H_2-H_1)^2}{4H_1H_2}$, negative for every physical jump.
- **Check.** H₁ = 0.1 m, Fr₁ = 3: E₂ − E₁ = −1.385 J/kg, head loss 0.141 m (the published Bélanger loss) ✓. Units J/kg ✓.
- **What it means.** The jump is a shock: like a gas-dynamic shock (Ch. 15) it takes supercritical to subcritical flow and
  loses mechanical energy; engineers build stilling basins to make jumps dissipate a spillway's energy safely.
- **Traps.** Taking E per unit mass of a particle at the bed (the pressure term changes); losing the minus sign; thinking
  the loss is small (it grows like the cube of the depth change).

### D27 · Deep-water Stokes drift: (7.84a) → (7.85) — ★★, 9 steps, in C12 (notebook · `particle_orbits`)
- **Goal.** Find the mean forward speed of a water parcel under a deep-water wave, to second order in the amplitude.
- **Start.** $\frac{dx_p}{dt}=u(x_0,z_0,t)+\xi\big(\frac{\partial u}{\partial x}\big)_{x_0,z_0}+\zeta\big(\frac{\partial u}{\partial z}\big)_{x_0,z_0}+\ldots$
  (7.84a) — *in words:* the velocity a parcel feels, corrected for where it actually is.
- **Plan.** (1) Average the Taylor-expanded velocity over a period. (2) Insert the linear orbits and velocities. (3) Watch
  sin² + cos² make the average constant.
- **Tools.** multivariable Taylor (P98) · mean of sin² and cos² (P151) · Lagrangian vs Eulerian mean (gloss) · D10, (7.46),
  (7.47).
- **Assumptions.** Deep water; second order in ka; the first-order orbit (7.46) is used in the correction terms.
- **Steps.**
  1. *did:* Expand u about the mean position · *tex:* $\frac{dx_p}{dt}=u(x_0,z_0,t)+\xi\frac{\partial u}{\partial x}+\zeta\frac{\partial u}{\partial z}
     +\ldots$ (7.84a) · *why:* First-order multivariable Taylor (P98) in the small excursions ξ, ζ — this time we keep the terms
     D10 step 4 dropped. · *plain:* the velocity the parcel really feels.
  2. *did:* Define the drift · *tex:* $\bar u_L=\Big\langle\frac{dx_p}{dt}\Big\rangle$ · *why:* The mean velocity following one parcel
     over a period (the Lagrangian mean). · *plain:* how fast the parcel advances on average.
  3. *did:* Average the first term · *tex:* $\langle u(x_0,z_0,t)\rangle=0$ · *why:* At a fixed point u ∝ cos(kx₀ − ωt), which
     averages to zero over a period. · *plain:* the linear velocity alone gives no drift.
  4. *did:* Insert the deep-water fields · *tex:* $\xi=-ae^{kz_0}\sin\varphi,\ \zeta=ae^{kz_0}\cos\varphi,\ u=a\omega e^{kz}\cos\varphi,\
     \varphi=kx_0-\omega t$ · *why:* (7.46) and (7.47) at the mean position. · *plain:* the known orbit and velocity.
  5. *did:* Compute ξ ∂u/∂x · *tex:* $\frac{\partial u}{\partial x}=-a\omega ke^{kz_0}\sin\varphi\ \Rightarrow\ \xi\frac{\partial u}{\partial x}=a^2\omega
     ke^{2kz_0}\sin^2\varphi$ · *why:* d/dx cos(kx − ωt) = −k sin(kx − ωt); multiplying two negative factors gives a positive
     product. · *plain:* ahead of its mean position the parcel meets faster forward water.
  6. *did:* Compute ζ ∂u/∂z · *tex:* $\frac{\partial u}{\partial z}=a\omega ke^{kz_0}\cos\varphi\ \Rightarrow\ \zeta\frac{\partial u}{\partial z}=a^2\omega
     ke^{2kz_0}\cos^2\varphi$ · *why:* Differentiate (7.47) in z: d/dz e^{kz} = k e^{kz}; then multiply by ζ from (7.46), both at the mean position. · *plain:* the parcel is highest, in faster water, when u is
     forward. · *live (E2):* "⟨ξu_x + ζu_z⟩ = …".
  7. *did:* Add the two terms · *tex:* $\xi\frac{\partial u}{\partial x}+\zeta\frac{\partial u}{\partial z}=a^2\omega ke^{2kz_0}\big(\sin^2\varphi+\cos^2\varphi\big)
     =a^2\omega ke^{2kz_0}$ · *why:* sin² + cos² = 1: the sum does not depend on time at all. · *plain:* the correction always
     points forward. · *watch (E2):* "the purple track creeps forward at this slope".
  8. *did:* Average · *tex:* $\bar u_L=a^2\omega ke^{2kz_0}$ (7.85) · *why:* The average of a constant is the constant — no factor ½
     appears. · *plain:* the deep-water Stokes drift.
  9. *did:* Compare the decay rates · *tex:* $\bar u_L\propto e^{2kz_0}\quad\text{vs}\quad\text{orbit radius}\propto e^{kz_0}$ · *why:* The drift
     is a product of two first-order fields, each ∝ e^{kz₀}. · *plain:* the drift dies twice as fast with depth; at the
     surface it is ka × aω.
- **Result.** $\bar u_L=a^2\omega ke^{2kz_0}$ (7.85); at any depth $\bar u_L=a^2\omega k\frac{\cosh2k(z_0+H)}{2\sinh^2kH}$ (7.86) by the
  same moves (N91).
- **Check.** Units: a²ωk [m²·1/s·1/m] = m/s ✓. a = 1 m, T = 8 s: 4.94 cm/s ✓. Exact path lines (C12 from scratch) drift at
  (7.85) within 2 % at ka = 0.05 ✓. H → ∞ in (7.86) gives (7.85) ✓.
- **What it means.** Waves carry water — a mass transport confined to the top ~λ/4, while a fixed meter reads no mean
  current (N92). In the ocean it drives the Stokes–Coriolis force and Langmuir cells (Ch. 13).
- **Traps.** Expecting a ½ from averaging (sin² + cos² is constant); averaging u at a fixed point (zero); using e^{kz₀} for
  the decay.

### D28 · Interfacial waves: (7.89)–(7.94) → $\omega=\varepsilon\sqrt{gk}$ (7.95) — ★★, 10 steps, in C13 (notebook · `two_layer_modes`)
- **Goal.** Find the frequency of a wave on the interface between two deep fluids of densities ρ₁ (above) < ρ₂ (below).
- **Start.** $\zeta=ae^{i(kx-\omega t)}$ (7.89), Laplace in each fluid (7.90), $\phi_1\to0$ as $z\to\infty$ (7.91), $\phi_2\to0$ as
  $z\to-\infty$ (7.92), $\frac{\partial\phi_1}{\partial z}=\frac{\partial\phi_2}{\partial z}=\frac{\partial\zeta}{\partial t}$ (7.93), $\rho_1\frac{\partial
  \phi_1}{\partial t}+\rho_1g\zeta=\rho_2\frac{\partial\phi_2}{\partial t}+\rho_2g\zeta$ (7.94), both at z = 0.
- **Plan.** (1) Choose decaying potentials. (2) Two kinematic conditions give A and B. (3) The pressure condition gives ω.
- **Tools.** complex amplitudes (P176) · e^{±kz} from Laplace (as D05) · matching coefficients.
- **Assumptions.** Both fluids infinitely deep; linear; no interfacial tension; ρ₂ > ρ₁ (stable).
- **Steps.**
  1. *did:* Write the interface in complex form · *tex:* $\zeta=ae^{i(kx-\omega t)}$ (7.89) · *why:* Complex notation (P176): Re is
     taken at the end; derivatives become multiplications. · *plain:* the interface displacement.
  2. *did:* Choose decaying potentials · *tex:* $\phi_1=Ae^{-kz}e^{i(kx-\omega t)},\quad\phi_2=Be^{kz}e^{i(kx-\omega t)}$ · *why:* Laplace
     gives e^{±kz} (as D05); (7.91) forbids e^{kz} above (it grows as z → ∞), (7.92) forbids e^{−kz} below. · *plain:* each
     fluid's motion dies away from the interface.
  3. *did:* Differentiate at z = 0 · *tex:* $\frac{\partial\phi_1}{\partial z}=-kA\,e^{i\theta},\ \frac{\partial\phi_2}{\partial z}=kB\,e^{i\theta},\
     \frac{\partial\zeta}{\partial t}=-i\omega a\,e^{i\theta}$ · *why:* d/dz e^{∓kz} = ∓k at z = 0; ∂/∂t brings −iω (P176); θ =
     kx − ωt. · *plain:* both vertical velocities and the interface's rise rate.
  4. *did:* Apply the upper kinematic condition · *tex:* $-kA=-i\omega a\ \Rightarrow\ A=\frac{i\omega a}k$ · *why:* (7.93) for fluid 1;
     cancel e^{iθ} (it holds for all x, t). · *plain:* the upper fluid follows the interface.
  5. *did:* Apply the lower kinematic condition · *tex:* $kB=-i\omega a\ \Rightarrow\ B=-\frac{i\omega a}k$ · *why:* The lower fluid also moves with the interface, (7.93) ∂φ₂/∂z = ∂ζ/∂t at z = 0; cancel e^{iθ} again.
     · *plain:* A = −B: the two potentials are opposite. · *live (E8):* "A = iωa/k = i × …".
  6. *did:* Differentiate the potentials in time · *tex:* $\frac{\partial\phi_1}{\partial t}=\frac{\omega^2a}k e^{i\theta},\quad\frac{\partial\phi_2}{\partial t}=
     -\frac{\omega^2a}k e^{i\theta}$ · *why:* ∂/∂t → −iω: (−iω)(iωa/k) = ω²a/k and (−iω)(−iωa/k) = −ω²a/k. · *plain:* the
     pressure-producing parts, opposite in the two fluids.
  7. *did:* Insert into the pressure condition · *tex:* $\rho_1\frac{\omega^2a}k+\rho_1ga=-\rho_2\frac{\omega^2a}k+\rho_2ga$ · *why:* (7.94) with
     ζ = ae^{iθ}; cancel e^{iθ}. · *plain:* pressure matches across the interface.
  8. *did:* Collect the ω² terms · *tex:* $\frac{\omega^2}k(\rho_1+\rho_2)=g(\rho_2-\rho_1)$ · *why:* Move ω² terms left, g terms right;
     divide by a ≠ 0. · *plain:* the inertia of both fluids against the weight difference.
  9. *did:* Solve for ω · *tex:* $\omega=\sqrt{gk\frac{\rho_2-\rho_1}{\rho_2+\rho_1}}=\varepsilon\sqrt{gk}$ (7.95), $\varepsilon^2\equiv\frac{\rho_2-\rho_1}
     {\rho_2+\rho_1}$ · *why:* Multiply by k/(ρ₁ + ρ₂); positive root. · *plain:* a deep-water wave slowed by ε. · *live (E8):*
     "ω² = gk(ρ₂ − ρ₁)/(ρ₂ + ρ₁) = …".
  10. *did:* Check the air–water limit · *tex:* $\frac{\rho_1}{\rho_2}\to0\ \Rightarrow\ \varepsilon\to1,\ \omega=\sqrt{gk}$ · *why:* With a
      negligible upper density the interface is a free surface: (7.45). · *plain:* the surface wave is recovered.
- **Result.** $\omega=\varepsilon\sqrt{gk}$, $\varepsilon^2=\frac{\rho_2-\rho_1}{\rho_2+\rho_1}$ (7.95); u₁ = −ωae^{−kz}e^{iθ}, u₂ = ωae^{kz}e^{iθ}
  (N107).
- **Check.** Units ✓. Δρ = 2 on 1000 kg/m³, λ = 100 m: T = 253 s vs 8.0 s at the surface ✓. `np.linalg.solve` from
  scratch agrees ✓. ρ₁ > ρ₂ gives ω² < 0: growth (Rayleigh–Taylor, Ch. 11).
- **What it means.** Weak stratification makes slow, tall waves; the velocity jumps across the interface — a vortex sheet,
  the seed of Kelvin–Helmholtz instability.
- **Traps.** Letting φ₁ grow upward; forgetting the minus from −iω; cancelling ρ₁ terms with the wrong sign.

### D29 · The two-layer constants (7.104)–(7.109) — ★★★, 12 steps, in C14 (notebook · `two_layer_modes`)
- **Goal.** For a layer of thickness H with a free surface over a deep fluid, express the flow and the interface amplitude
  through the surface amplitude a and ω.
- **Start.** $\phi_2\to0$ at $z\to-\infty$ (7.97); $\frac{\partial\phi_1}{\partial z}=\frac{\partial\eta}{\partial t}$ at z = 0 (7.98);
  $\frac{\partial\phi_1}{\partial t}+g\eta=0$ at z = 0 (7.99); $\frac{\partial\phi_1}{\partial z}=\frac{\partial\phi_2}{\partial z}=\frac{\partial\zeta}{\partial t}$ at
  z = −H (7.100); $\eta=ae^{i(kx-\omega t)}$ (7.102), $\zeta=be^{i(kx-\omega t)}$ (7.103).
- **Plan.** (1) Write the potentials. (2) Surface conditions give A, B. (3) Interface kinematics gives C and b.
- **Tools.** complex amplitudes (P176) · two linear equations (as D05) · sympy for complex systems (P158).
- **Assumptions.** Linear; lower layer infinitely deep; a real; b complex (a phase difference is allowed).
- **Steps.**
  1. *did:* Write the complex surfaces · *tex:* $\eta=ae^{i\theta},\quad\zeta=be^{i\theta},\quad\theta=kx-\omega t$ (7.102), (7.103) · *why:*
     Both boundaries oscillate with the same k and ω; a real sets the phase origin; b may be complex. · *plain:* surface
     and interface waves.
  2. *did:* Write the layer potentials · *tex:* $\phi_1=(Ae^{kz}+Be^{-kz})e^{i\theta}$ (7.104), $\phi_2=Ce^{kz}e^{i\theta}$ (7.105) · *why:* Laplace
     gives e^{±kz}; the finite upper layer keeps both; (7.97) keeps only e^{kz} below. The book prints (7.105) with
     e^{i(kz−ωt)}: read e^{i(kx−ωt)}. · *plain:* four unknowns A, B, C, b.
  3. *did:* Apply the surface kinematic condition · *tex:* $k(A-B)=-i\omega a$ · *why:* (7.98) at z = 0: ∂φ₁/∂z = k(A − B)e^{iθ},
     ∂η/∂t = −iωae^{iθ}; cancel e^{iθ}. · *plain:* the upper layer follows the free surface.
  4. *did:* Apply the surface dynamic condition · *tex:* $-i\omega(A+B)+ga=0\ \Rightarrow\ A+B=-\frac{iga}{\omega}$ · *why:* (7.99) at z = 0:
     ∂φ₁/∂t = −iω(A + B)e^{iθ}; 1/i = −i. · *plain:* the free surface stays at atmospheric pressure.
  5. *did:* Solve for A · *tex:* $A=-\frac{ia}2\Big(\frac\omega k+\frac g\omega\Big)$ (7.106) · *why:* A = ½[(A − B) + (A + B)] with A − B =
     −iωa/k. · *plain:* first constant.
  6. *did:* Solve for B · *tex:* $B=\frac{ia}2\Big(\frac\omega k-\frac g\omega\Big)$ (7.107) · *why:* B = ½[(A + B) − (A − B)]. · *plain:*
     second constant.
  7. *did:* Match vertical velocities at the interface · *tex:* $k(Ae^{-kH}-Be^{kH})=kCe^{-kH}$ · *why:* (7.100), first equality at
     z = −H: ∂φ₁/∂z = k(Ae^{kz} − Be^{−kz}) and ∂φ₂/∂z = kCe^{kz}, both at z = −H. · *plain:* the layers move together
     vertically at the interface. · *live (E8):* "A = −(ia/2)(ω/k + g/ω) = …".
  8. *did:* Solve for C · *tex:* $C=A-Be^{2kH}=-\frac{ia}2\Big(\frac\omega k+\frac g\omega\Big)-\frac{ia}2\Big(\frac\omega k-\frac g\omega\Big)e^{2kH}$
     (7.108) · *why:* Divide by ke^{−kH}; insert A and B. · *plain:* third constant.
  9. *did:* Match with the interface motion · *tex:* $k(Ae^{-kH}-Be^{kH})=-i\omega b$ · *why:* (7.100), second equality: ∂φ₁/∂z =
     ∂ζ/∂t = −iωbe^{iθ} at z = −H. · *plain:* the interface moves with the water.
  10. *did:* Solve for b · *tex:* $b=\frac{ik}\omega\big(Ae^{-kH}-Be^{kH}\big)$ · *why:* Divide by −iω; 1/(−i) = i. · *plain:* the
      interface amplitude from A and B.
  11. *did:* Insert A and B · *tex:* $\frac{ik}\omega Ae^{-kH}=\frac a2\Big(1+\frac{gk}{\omega^2}\Big)e^{-kH},\quad-\frac{ik}\omega Be^{kH}=\frac a2\Big(1-
      \frac{gk}{\omega^2}\Big)e^{kH}$ · *why:* i·(−i) = 1 and −i·i = 1: (k/ω)(a/2)(ω/k + g/ω) = (a/2)(1 + gk/ω²), and likewise the
      second term. · *plain:* two pieces of b.
  12. *did:* Add them · *tex:* $b=\frac a2\Big(1+\frac{gk}{\omega^2}\Big)e^{-kH}+\frac a2\Big(1-\frac{gk}{\omega^2}\Big)e^{kH}$ (7.109) · *why:* Sum of
      step 11; b is real here (the modes will be exactly in or out of phase). ω is still unknown: (7.101) fixes it (D30). ·
      *plain:* the interface amplitude. · *live (E8):* "b/a = …".
- **Result.** (7.106)–(7.109).
- **Check.** Units: A, B, C [m²/s]; b [m] ✓. H → 0 (no upper layer): b → a — the interface is the surface ✓.
- **sympy check intent (★★★).** Re-do steps 3–12 symbolically and test every condition; show the printed (7.105) fails.
  `check_src` sketch:
  ```python
  import sympy as sp                                         # symbolic algebra with complex numbers
  x, z, t = sp.symbols('x z t', real=True)                   # coordinates and time
  k, w, g, H, a = sp.symbols('k omega g H a', positive=True) # wavenumber, frequency, gravity, layer depth, surface amplitude
  A, B, C, b = sp.symbols('A B C b')                         # unknown (complex) constants
  E = sp.exp(sp.I*(k*x - w*t))                               # the common wave factor e^{iθ}
  phi1 = (A*sp.exp(k*z) + B*sp.exp(-k*z))*E                  # (7.104) upper-layer potential
  phi2 = C*sp.exp(k*z)*E                                     # (7.105) with e^{i(kx-ωt)} (the printed kz is a slip)
  eta, zeta = a*E, b*E                                       # (7.102), (7.103)
  eqs = [sp.diff(phi1, z).subs(z, 0) - sp.diff(eta, t),      # (7.98) surface kinematics
         sp.diff(phi1, t).subs(z, 0) + g*eta,                # (7.99) surface dynamics
         (sp.diff(phi1, z) - sp.diff(phi2, z)).subs(z, -H),  # (7.100) first equality
         sp.diff(phi1, z).subs(z, -H) - sp.diff(zeta, t)]    # (7.100) second equality
  sol = sp.solve([sp.simplify(e/E) for e in eqs], [A, B, C, b], dict=True)[0]   # steps 5-12
  book = {A: -sp.I*a/2*(w/k + g/w), B: sp.I*a/2*(w/k - g/w),
          C: -sp.I*a/2*(w/k + g/w) - sp.I*a/2*(w/k - g/w)*sp.exp(2*k*H),
          b: a/2*(1 + g*k/w**2)*sp.exp(-k*H) + a/2*(1 - g*k/w**2)*sp.exp(k*H)}   # (7.106)-(7.109)
  print([sp.simplify(sol[s] - book[s]) for s in (A, B, C, b)])                 # [0, 0, 0, 0]
  phi2_printed = C*sp.exp(k*z)*sp.exp(sp.I*(k*z - w*t))                          # the printed (7.105)
  r = (sp.diff(phi1, z) - sp.diff(phi2_printed, z)).subs(z, -H).subs(sol)
  print(sp.simplify(sp.diff(r, x)) == 0)                     # False: it cannot hold for every x
  ```
  (`ch07.two_layer_sympy()` does the same and is what the notebook cell prints; the cell shows the construction.)
- **What it means.** Everything is fixed by a and ω; the one remaining condition (pressure at the interface) will allow
  only two values of ω — two modes.
- **Traps.** Using the printed (7.105); forgetting that b may be complex until the algebra shows it is real; signs of i.

### D30 · The two-layer dispersion relation (7.110) — ★★★, 12 steps, in C14 (notebook · `two_layer_modes`)
- **Goal.** Use pressure continuity at the interface to find which ω are possible, and show the result factors into two
  modes. (The book: "after some algebraic manipulations", Exercise 7.19; written out here.)
- **Start.** $\rho_1\frac{\partial\phi_1}{\partial t}+\rho_1g\zeta=\rho_2\frac{\partial\phi_2}{\partial t}+\rho_2g\zeta$ at $z=-H$ (7.101) with
  (7.106)–(7.109).
- **Plan.** (1) Evaluate (7.101) with the constants. (2) Eliminate C, then A and B. (3) Introduce s = ω²/gk. (4) Factor out
  (s − 1) and turn exponentials into sinh and cosh.
- **Tools.** complex amplitudes (P176) · hyperbolic functions from exponentials (P168) · factorising (P158 in sympy) · D29.
- **Assumptions.** As D29; no interfacial tension.
- **Steps.**
  1. *did:* Write the interface pressure condition · *tex:* $\rho_1\frac{\partial\phi_1}{\partial t}+\rho_1g\zeta=\rho_2\frac{\partial\phi_2}{\partial t}+
     \rho_2g\zeta$ at $z=-H$ (7.101) · *why:* Pressure is continuous across the interface (no tension); linear Bernoulli in each
     layer, applied at the mean interface. · *plain:* the last condition — it fixes ω.
  2. *did:* Evaluate the time derivatives · *tex:* $\frac{\partial\phi_1}{\partial t}=-i\omega(Ae^{-kH}+Be^{kH})e^{i\theta},\quad\frac{\partial\phi_2}{\partial t}=
     -i\omega Ce^{-kH}e^{i\theta}$ · *why:* ∂/∂t → −iω; set z = −H in (7.104), (7.105). · *plain:* the pressure-producing parts.
  3. *did:* Collect the terms · *tex:* $-i\omega\big[\rho_1(Ae^{-kH}+Be^{kH})-\rho_2Ce^{-kH}\big]=(\rho_2-\rho_1)gb$ · *why:* Cancel e^{iθ}; move
     ρ₁gζ to the right. · *plain:* inertia difference = weight difference.
  4. *did:* Replace Ce^{−kH} · *tex:* $Ce^{-kH}=Ae^{-kH}-Be^{kH}$ · *why:* D29 step 7 (vertical velocities match). · *plain:*
     eliminate C.
  5. *did:* Regroup by A and B · *tex:* $-i\omega\big[(\rho_1-\rho_2)Ae^{-kH}+(\rho_1+\rho_2)Be^{kH}\big]=(\rho_2-\rho_1)gb$ · *why:* Expand ρ₁Ae^{−kH}
     + ρ₁Be^{kH} − ρ₂Ae^{−kH} + ρ₂Be^{kH}. · *plain:* two terms left.
  6. *did:* Insert −iωA and −iωB · *tex:* $-i\omega A=-\frac a2\Big(\frac{\omega^2}k+g\Big),\quad-i\omega B=\frac a2\Big(\frac{\omega^2}k-g\Big)$ · *why:*
     From (7.106), (7.107): (−iω)(−ia/2) = −ωa/2 and (−iω)(ia/2) = ωa/2. · *plain:* the constants times −iω.
  7. *did:* Introduce s = ω²/gk · *tex:* $\frac{\omega^2}k=gs,\qquad\frac{gk}{\omega^2}=\frac1s$ · *why:* One dimensionless unknown shortens
     the algebra; s = 1 is the deep-water surface wave. · *plain:* measure ω² in units of gk.
  8. *did:* Write the equation in s · *tex:* $(\rho_2-\rho_1)(s+1)e^{-kH}+(\rho_1+\rho_2)(s-1)e^{kH}=(\rho_2-\rho_1)\big[(1+\tfrac1s)e^{-kH}+(1-\tfrac1s)e^{kH}\big]$ ·
     *why:* Steps 5–6 with b from (7.109); multiply by 2/(ag); (ρ₁ − ρ₂)(−(s + 1)) = (ρ₂ − ρ₁)(s + 1). · *plain:* one equation for s.
  9. *did:* Move all to the left and factor · *tex:* $\frac{s-1}{s}\Big[(\rho_2-\rho_1)(s+1)e^{-kH}+\big((\rho_1+\rho_2)s-(\rho_2-\rho_1)\big)e^{kH}\Big]=0$ ·
     *why:* (s + 1) − (1 + 1/s) = (s − 1)(s + 1)/s, and (ρ₁ + ρ₂)(s − 1) − (ρ₂ − ρ₁)(1 − 1/s) = ((s − 1)/s)[(ρ₁ + ρ₂)s − (ρ₂ − ρ₁)].
     · *plain:* a common factor (s − 1)/s appears.
  10. *did:* Keep the factor s − 1 · *tex:* $s-1=\frac{\omega^2}{gk}-1$ · *why:* Do not cancel it as 'trivial': s = 1, ω² = gk, is a
      genuine solution — the surface wave (D31). 1/s ≠ 0 may be dropped. · *plain:* one factor per mode. · *set (E8):* the
      barotropic branch lights up.
  11. *did:* Turn the bracket into sinh and cosh · *tex:* $s\big[(\rho_2-\rho_1)e^{-kH}+(\rho_1+\rho_2)e^{kH}\big]+(\rho_2-\rho_1)(e^{-kH}-e^{kH})
      =2s(\rho_1\sinh kH+\rho_2\cosh kH)-2(\rho_2-\rho_1)\sinh kH$ · *why:* Collect the s-terms; ρ₂(e^{kH} + e^{−kH}) + ρ₁(e^{kH} −
      e^{−kH}) = 2ρ₂ cosh kH + 2ρ₁ sinh kH; e^{−kH} − e^{kH} = −2 sinh kH (P168). · *plain:* exponentials become hyperbolic
      functions.
  12. *did:* Assemble · *tex:* $\Big(\frac{\omega^2}{gk}-1\Big)\Big\{\frac{\omega^2}{gk}[\rho_1\sinh kH+\rho_2\cosh kH]-(\rho_2-\rho_1)\sinh kH\Big\}=0$ (7.110) ·
      *why:* Divide by the nonzero factor 2/s. The full residual of (7.101) is (ag²k/ω²) × this left side. · *plain:* the
      two-layer dispersion relation. · *set (E8):* the baroclinic branch lights up.
- **Result.** (7.110).
- **Check.** kH → ∞: the brace gives ω² = gk(ρ₂ − ρ₁)/(ρ₂ + ρ₁) (7.95) ✓. ρ₁ = ρ₂: the brace gives ω²(ρ sinh + ρ cosh) = 0 → no
  internal mode ✓. `np.roots` of the expanded quadratic in s agrees with `two_layer_free_surface_omega` to 1e-10 ✓.
- **sympy check intent (★★★).** Build the residual of (7.101) with the constants of D29, factor it, and compare with the
  book's product. `check_src` sketch:
  ```python
  import sympy as sp                                              # symbolic algebra
  k, w, g, H, a, r1, r2 = sp.symbols('k omega g H a rho1 rho2', positive=True)
  A = -sp.I*a/2*(w/k + g/w); B = sp.I*a/2*(w/k - g/w)             # (7.106), (7.107)
  C = A - B*sp.exp(2*k*H)                                         # (7.108), D29 step 8
  b = sp.I*k/w*(A*sp.exp(-k*H) - B*sp.exp(k*H))                   # D29 step 10
  lhs = -sp.I*w*(r1*(A*sp.exp(-k*H) + B*sp.exp(k*H)) - r2*C*sp.exp(-k*H))   # step 3, left side
  res = lhs - (r2 - r1)*g*b                                       # residual of (7.101) after cancelling e^{iθ}
  s = w**2/(g*k)                                                  # step 7
  book = (s - 1)*(s*(r1*sp.sinh(k*H) + r2*sp.cosh(k*H)) - (r2 - r1)*sp.sinh(k*H))   # (7.110)
  print(sp.simplify((res - a*g**2*k/w**2*book).rewrite(sp.exp)))  # 0: residual = (ag²k/ω²) × (7.110)
  print(sp.factor(sp.simplify((res*w**2/(a*g**2*k)).rewrite(sp.exp))))   # shows the factor (ω² − gk) explicitly
  ```
- **What it means.** One condition, two roots: a layered fluid supports two kinds of wave with the same wavelength.
- **Traps.** Cancelling (s − 1) as trivial; dropping the common factor ag²k/ω² without noting it is nonzero; sign slips
  with i² = −1.

### D31 · Roots, mode shapes and long waves: (7.111)–(7.119), $c=\sqrt{g'H}$ — ★★, 11 steps, in C14 (notebook · `two_layer_modes`)
- **Goal.** Read the two modes out of (7.110): their frequencies, how the interface moves relative to the surface, and
  their long-wave speeds.
- **Start.** (7.110) and $b=\frac a2\big(1+\frac{gk}{\omega^2}\big)e^{-kH}+\frac a2\big(1-\frac{gk}{\omega^2}\big)e^{kH}$ (7.109).
- **Plan.** (1) First factor: the surface mode. (2) Second factor: the internal mode and its amplitude ratio. (3) Deep and
  long-wave limits.
- **Tools.** small-argument forms (P168) · D29, D30.
- **Assumptions.** kH ≪ 1 for the long-wave steps (9–11); lower layer deep.
- **Steps.**
  1. *did:* Take the first factor · *tex:* $\omega^2=gk$ (7.111) · *why:* s = 1 zeroes the first factor of (7.110). · *plain:* the
     deep-water surface wave, blind to the density step.
  2. *did:* Put it into b · *tex:* $b=\frac a2(2)e^{-kH}+0=ae^{-kH}$ (7.112) · *why:* In (7.109) gk/ω² = 1: 1 + 1/s = 2 and 1 − 1/s = 0. ·
     *plain:* the interface moves in phase, smaller by e^{−kH}.
  3. *did:* Name the mode · *tex:* $\frac\eta\zeta=e^{kH}>0\ \text{(barotropic)}$ · *why:* Surfaces of constant pressure and density
     move together; the motion decays like e^{kz} from the top, as a surface wave. · *plain:* the barotropic mode. · *set
     (E8):* mode = barotropic.
  4. *did:* Take the second factor · *tex:* $\omega^2=\frac{gk(\rho_2-\rho_1)\sinh kH}{\rho_2\cosh kH+\rho_1\sinh kH}$ (7.113) · *why:* Set the
     brace of (7.110) to zero and solve for s = ω²/gk. · *plain:* the second, slow mode.
  5. *did:* Evaluate the brackets of b · *tex:* $1+\frac1s=\frac{\rho_2e^{kH}}{D},\quad1-\frac1s=-\frac{\rho_2e^{-kH}+2\rho_1\sinh kH}{D},\quad D=(\rho_2-\rho_1)\sinh kH$ ·
     *why:* 1/s = (ρ₁ sinh kH + ρ₂ cosh kH)/D; add or subtract D and use sinh ± cosh = ±e^{±kH}. · *plain:* the two brackets
     of (7.109). · *set (E8):* mode = baroclinic.
  6. *did:* Assemble b · *tex:* $b=\frac{a}{2D}\big[\rho_2-\rho_2-2\rho_1\sinh kH\,e^{kH}\big]=-\frac{a\rho_1e^{kH}}{\rho_2-\rho_1}$ · *why:*
     (a/2)(ρ₂e^{kH}/D)e^{−kH} = aρ₂/2D and (a/2)(−(ρ₂e^{−kH} + 2ρ₁ sinh kH)/D)e^{kH} = −a(ρ₂ + 2ρ₁ sinh kH e^{kH})/2D; the ρ₂ cancel and
     sinh kH cancels with D. · *plain:* the interface amplitude of the slow mode.
  7. *did:* Form the ratio · *tex:* $\eta=-\zeta\Big(\frac{\rho_2-\rho_1}{\rho_1}\Big)e^{-kH}$ (7.114) · *why:* The amplitude ratio of the two surfaces is η/ζ = a/b because both share the same factor e^{i(kx − ωt)}. · *plain:* antiphase,
     and a tiny surface when Δρ is small: the baroclinic mode.
  8. *did:* Check the deep limit · *tex:* $kH\to\infty:\ \omega^2\to gk\frac{\rho_2-\rho_1}{\rho_2+\rho_1}$ (7.95) · *why:* sinh kH ≈ cosh kH ≈
     e^{kH}/2 cancel. · *plain:* far from the surface the interface wave of C13 returns.
  9. *did:* Take long waves · *tex:* $\omega^2\approx\frac{gk(\rho_2-\rho_1)kH}{\rho_2+\rho_1kH}\approx kg\Big(\frac{\rho_2-\rho_1}{\rho_2}\Big)kH$ (7.115) ·
     *why:* sinh kH ≈ kH, cosh kH ≈ 1 (P168); ρ₁kH ≪ ρ₂ at lowest order in kH. · *plain:* long baroclinic waves. · *set (E8):*
     long waves.
  10. *did:* Read off the speed · *tex:* $c=\frac\omega k=[g'H]^{1/2},\quad g'=g\Big(\frac{\rho_2-\rho_1}{\rho_2}\Big)$ (7.116, 7.117) · *why:*
      ω² = g′H k², so c² = g′H. ⚠️ g′ has ρ₂ below the line here (Ch. 4 used ρ₁). · *plain:* the shallow-water speed with
      reduced gravity. · *live (E8):* "c = √(g′H) = … m/s".
  11. *did:* Long-wave ratio and pressure · *tex:* $\eta=-\zeta\frac{\rho_2-\rho_1}{\rho_1}$ (7.118), $p'=-\rho_1\frac{\partial\phi_1}{\partial t}=i\rho_1\omega
      (A+B)e^{i\theta}=\rho_1g\eta$ (7.119) · *why:* e^{−kH} → 1 in (7.114); in a thin layer e^{±kz} ≈ 1 so φ₁ ≈ (A + B)e^{iθ},
      and A + B = −iga/ω (D29 step 4). · *plain:* the upper layer is hydrostatic for long waves.
- **Result.** Barotropic ω² = gk, b = ae^{−kH}; baroclinic (7.113), η/ζ = −((ρ₂ − ρ₁)/ρ₁)e^{−kH}; long waves $c=\sqrt{g'H}$.
- **Check.** H = 50 m, ρ₁ = 1000, ρ₂ = 1002 kg/m³: c = 0.99 m/s vs √(gH) = 22.1 m/s ✓; a 10 m thermocline wave shows 2 cm at
  the surface ✓.
- **What it means.** The vocabulary of layered ocean and atmosphere models: barotropic (depth-uniform, fast, seen at the
  surface) and baroclinic (sheared, slow, hidden); √(g′H)/f is the baroclinic Rossby radius (Ch. 13).
- **Traps.** Which ρ in g′ (ρ₂ here); the sign of η/ζ (negative); taking kH → ∞ of the long-wave formula.

### D32 · The linear perturbation equations (7.123)–(7.131) — ★★, 10 steps, in C15 (notebook)
- **Goal.** From the Boussinesq equations, derive the linear equations for small motions about a resting stratified
  fluid.
- **Start.** $D\rho/Dt=0$ (4.9), (4.10), (7.120)–(7.122), the resting state (7.123) and $p=\bar p(z)+p'$, $\rho=\bar\rho(z)+\rho'$ (7.124).
- **Plan.** (1) Split the density equation and linearise. (2) Introduce N². (3) Subtract the hydrostatic base from the
  momentum equations.
- **Tools.** material derivative (Ch. 3) · linearisation about a base state (gloss; P68) · R16, R17, R18.
- **Assumptions.** Small amplitude (steps 4, 9); ρ̄ depends on z only (step 3); Boussinesq (ρ₀ in the inertia).
- **Steps.**
  1. *did:* Write the density equation in full · *tex:* $\frac{\partial\rho}{\partial t}+u\frac{\partial\rho}{\partial x}+v\frac{\partial\rho}{\partial y}+w\frac{\partial
     \rho}{\partial z}=0$ · *why:* D/Dt = ∂/∂t + u·∇ (Ch. 3); (4.9) says each parcel keeps its density. · *plain:* density is
     carried with the flow.
  2. *did:* Insert the split · *tex:* $\frac{\partial}{\partial t}(\bar\rho+\rho')+u\frac{\partial}{\partial x}(\bar\rho+\rho')+v\frac{\partial}{\partial y}(\bar\rho+\rho')+
     w\frac{\partial}{\partial z}(\bar\rho+\rho')=0$ (7.125) · *why:* (7.124) ρ = ρ̄(z) + ρ′. · *plain:* background plus perturbation.
  3. *did:* Use ρ̄ = ρ̄(z) · *tex:* $\frac{\partial\bar\rho}{\partial t}=\frac{\partial\bar\rho}{\partial x}=\frac{\partial\bar\rho}{\partial y}=0$ · *why:* The resting
     stratification is steady and horizontally uniform. · *plain:* three background terms vanish.
  4. *did:* Drop products of small quantities · *tex:* $\frac{\partial\rho'}{\partial t}+w\frac{d\bar\rho}{dz}=0$ (7.126) · *why:* u∂ρ′/∂x,
     v∂ρ′/∂y, w∂ρ′/∂z are small × small; w dρ̄/dz is small × O(1) and stays (linearisation). · *plain:* density changes at a
     point because vertical motion carries the background up or down.
  5. *did:* Introduce N² · *tex:* $\frac{d\bar\rho}{dz}=-\frac{\rho_0N^2}{g}\ \Rightarrow\ \frac{\partial\rho'}{\partial t}-\frac{N^2\rho_0}gw=0$ (7.131) ·
     *why:* The buoyancy frequency is defined by (7.127) N² ≡ −(g/ρ₀)dρ̄/dz (R18); substituting it names the stratification's strength. · *plain:* the density equation in terms of N.
  6. *did:* Split the vertical momentum · *tex:* $\frac{\partial w}{\partial t}=-\frac1{\rho_0}\frac{\partial(\bar p+p')}{\partial z}-\frac{(\bar\rho+\rho')g}
     {\rho_0}$ · *why:* Insert the split p = p̄(z) + p′, ρ = ρ̄(z) + ρ′ (7.124) into the vertical momentum equation (7.122) to separate background and wave. · *plain:* background and perturbation parts.
  7. *did:* Subtract the hydrostatic base · *tex:* $\frac{\partial w}{\partial t}=-\frac1{\rho_0}\frac{\partial p'}{\partial z}-\frac{\rho'g}{\rho_0}$ (7.130) ·
     *why:* (7.123): −(1/ρ₀)dp̄/dz − ρ̄g/ρ₀ = 0 cancels exactly. · *plain:* only perturbations accelerate the water vertically.
  8. *did:* Split the horizontal momentum · *tex:* $\frac{\partial u}{\partial t}=-\frac1{\rho_0}\frac{\partial p'}{\partial x},\quad\frac{\partial v}{\partial t}=
     -\frac1{\rho_0}\frac{\partial p'}{\partial y}$ (7.128, 7.129) · *why:* (7.120)–(7.121) with p = p̄(z) + p′; p̄ does not depend on x
     or y. · *plain:* horizontal acceleration by the perturbation pressure.
  9. *did:* Confirm nothing else is nonlinear · *tex:* $\mathbf u\cdot\nabla\mathbf u\ \text{already dropped in (7.120)–(7.122)}$ · *why:* The
     momentum equations of §7.8 are linear from the start (small amplitude, R12, R14); only the density equation needed
     linearising. · *plain:* the momentum side is already linear.
  10. *did:* Count equations and unknowns · *tex:* $(4.10),\ (7.128)\text{–}(7.131):\ 5\ \text{equations for}\ u,v,w,p',\rho'$ · *why:* A closed
      linear system: ready to be reduced to one equation (D33). · *plain:* the linear internal-wave equations.
- **Result.** (7.126), (7.128)–(7.131) with (4.10).
- **Check.** Units: (7.131) kg/(m³ s) on both terms ✓. `boussinesq_linear_sympy()` residuals 0 ✓. With N = 0 (uniform
  density) ρ′ stays 0 and the system is plain incompressible flow ✓.
- **What it means.** Density perturbations are made only by moving fluid across the background gradient; that is what N
  measures.
- **Traps.** Dropping w dρ̄/dz (it is the whole physics); p′ here is measured from p̄(z), not from −ρgz as in (7.30).

### D33 · The w-equation: (7.132) → (7.134) — ★★★, 12 steps, in C15 (notebook · `internal_wave_beams`)
- **Goal.** Eliminate u, v, p′ and ρ′ to get one equation for the vertical velocity.
- **Start.** (4.10), (7.128)–(7.131).
- **Plan.** (1) Time-differentiate continuity and bring in the horizontal momentum: (7.132). (2) Time-differentiate the
  vertical momentum and bring in the density equation: (7.133). (3) Eliminate p′ between them.
- **Tools.** operator elimination for linear PDEs (P177) · mixed partials commute (P121).
- **Assumptions.** Linear, inviscid, Boussinesq, no rotation; N may depend on z (step 8).
- **Steps.**
  1. *did:* Differentiate continuity in time · *tex:* $\frac{\partial}{\partial t}\Big(\frac{\partial u}{\partial x}+\frac{\partial v}{\partial y}\Big)=-\frac{\partial^2w}
     {\partial z\,\partial t}$ · *why:* (4.10) holds at all times, so its time derivative is zero; move the w-term right (P177). ·
     *plain:* prepare to bring in the momentum equations.
  2. *did:* Swap the derivatives · *tex:* $\frac{\partial}{\partial x}\Big(\frac{\partial u}{\partial t}\Big)+\frac{\partial}{\partial y}\Big(\frac{\partial v}{\partial t}\Big)=
     -\frac{\partial^2w}{\partial z\,\partial t}$ · *why:* Mixed partials of smooth fields commute (P121). · *plain:* now ∂u/∂t and
     ∂v/∂t appear. · *set (E9):* the chip "∂/∂t" lights.
  3. *did:* Insert the horizontal momentum · *tex:* $-\frac1{\rho_0}\Big(\frac{\partial^2p'}{\partial x^2}+\frac{\partial^2p'}{\partial y^2}\Big)=-\frac{\partial^2w}{\partial z\,
     \partial t}$ · *why:* (7.128)–(7.129): ∂u/∂t = −p′_x/ρ₀, ∂v/∂t = −p′_y/ρ₀. · *plain:* u and v are gone.
  4. *did:* Name the horizontal Laplacian · *tex:* $\frac1{\rho_0}\nabla_H^2p'=\frac{\partial^2w}{\partial z\,\partial t}$ (7.132) · *why:* ∇_H² ≡ ∂²/∂x² +
     ∂²/∂y²; multiply by −1. · *plain:* horizontal pressure curvature drives vertical stretching.
  5. *did:* Differentiate the vertical momentum in time · *tex:* $\frac{\partial^2w}{\partial t^2}=-\frac1{\rho_0}\frac{\partial^2p'}{\partial t\,\partial z}-\frac g
     {\rho_0}\frac{\partial\rho'}{\partial t}$ · *why:* (7.130) holds at all times; ρ₀ and g are constants. · *plain:* prepare to
     eliminate ρ′.
  6. *did:* Insert the density equation · *tex:* $\frac g{\rho_0}\frac{\partial\rho'}{\partial t}=N^2w$ · *why:* The linear density equation (7.131) gives ∂ρ′/∂t = (N²ρ₀/g)w; multiplying by g/ρ₀ turns the buoyancy term into N²w. · *plain:*
     the buoyancy acceleration in terms of w. · *set (E9):* the chip "N²w" lights.
  7. *did:* Rearrange · *tex:* $\frac1{\rho_0}\frac{\partial^2p'}{\partial t\,\partial z}=-\frac{\partial^2w}{\partial t^2}-N^2w$ (7.133) · *why:* Put step 6 into step 5 and move the pressure term to the left; this is the book's (7.133). · *plain:* a second relation between p′ and w.
  8. *did:* Apply ∇_H² to (7.133) · *tex:* $\frac1{\rho_0}\frac{\partial^2}{\partial t\,\partial z}\nabla_H^2p'=-\nabla_H^2\Big(\frac{\partial^2w}{\partial t^2}+N^2w\Big)$ ·
     *why:* ∇_H² commutes with ∂/∂t, ∂/∂z and with N(z), since it acts only in x and y. · *plain:* prepare to use (7.132). · *set
     (E9):* the chip "∇_H²" lights.
  9. *did:* Apply ∂²/∂t∂z to (7.132) · *tex:* $\frac1{\rho_0}\frac{\partial^2}{\partial t\,\partial z}\nabla_H^2p'=\frac{\partial^2}{\partial t\,\partial z}\Big(\frac{\partial^2w}
     {\partial z\,\partial t}\Big)$ · *why:* Differentiate both sides of (7.132); ρ₀ is constant. · *plain:* the same left side as step 8.
  10. *did:* Equate the right sides · *tex:* $\frac{\partial^4w}{\partial t^2\partial z^2}=-\nabla_H^2\frac{\partial^2w}{\partial t^2}-N^2\nabla_H^2w$ · *why:*
      Steps 8 and 9 have identical left sides; commute the derivatives on the left of step 9. · *plain:* p′ is eliminated.
  11. *did:* Group the time derivatives · *tex:* $\frac{\partial^2}{\partial t^2}\Big(\frac{\partial^2w}{\partial z^2}+\nabla_H^2w\Big)+N^2\nabla_H^2w=0$ · *why:*
      Move everything left and factor ∂²/∂t². · *plain:* nearly there.
  12. *did:* Use the full Laplacian · *tex:* $\frac{\partial^2}{\partial t^2}\nabla^2w+N^2\nabla_H^2w=0$ (7.134) · *why:* ∇² = ∇_H² + ∂²/∂z². Valid for
      N(z) too (step 8). · *plain:* one equation for w. · *live (E9):* "residual with your wave = …".
- **Result.** $\frac{\partial^2}{\partial t^2}\nabla^2w+N^2\nabla_H^2w=0$ (7.134).
- **Check.** Units: 1/(m² s³) per unit w ✓. N = 0: ∇²w = const in time — no waves without stratification ✓. A plane wave with
  ω from (7.138) makes the finite-difference residual ≈ 0, a wrong ω does not (C15 from scratch) ✓.
- **sympy check intent (★★★).** Re-run the construction with generic functions: build (7.128)–(7.131) for symbolic u, v, w,
  p′, ρ′ (x, y, z, t) and N(z); express ∂u/∂t, ∂v/∂t, ∂ρ′/∂t through the equations; form the combination of steps 8–10 and
  show it vanishes identically given (4.10). `check_src` sketch:
  ```python
  import sympy as sp                                        # symbolic algebra
  x, y, z, t = sp.symbols('x y z t', real=True)             # coordinates and time
  rho0, g = sp.symbols('rho0 g', positive=True)             # reference density, gravity
  N = sp.Function('N')(z)                                   # buoyancy frequency, allowed to depend on z
  p = sp.Function('p')(x, y, z, t); w = sp.Function('w')(x, y, z, t)   # perturbation pressure and vertical velocity
  lap_H = lambda f: sp.diff(f, x, 2) + sp.diff(f, y, 2)     # horizontal Laplacian ∇_H²
  # (7.132) from continuity + (7.128)-(7.129): ∇_H² p = rho0 w_zt
  eq132 = lap_H(p)/rho0 - sp.diff(w, z, t)
  # (7.133) from (7.130) + (7.131): p_tz / rho0 = -w_tt - N² w
  eq133 = sp.diff(p, t, z)/rho0 + sp.diff(w, t, 2) + N**2*w
  # steps 8-10: ∂_t∂_z of (7.132) minus ∇_H² of (7.133) must equal -(7.134)
  combo = sp.diff(eq132, t, z) - lap_H(eq133)
  w134 = sp.diff(w, t, 2, z, 2) + lap_H(sp.diff(w, t, 2)) + N**2*lap_H(w)   # (7.134), ∇² = ∇_H² + ∂²/∂z²
  print(sp.simplify(combo + w134))                          # 0: the elimination is exact, even for N(z)
  ```
  (and `ch07.boussinesq_linear_sympy()` returns every step's residual, printed in the cell.)
- **What it means.** One scalar equation governs internal waves; its operator treats horizontal and vertical directions
  differently — the source of ω = N cos θ (D34).
- **Traps.** Forgetting that ∂/∂z and N(z) do not commute (we never needed to swap them); dropping the minus sign when
  moving w_zt; losing a ρ₀.

### D34 · Internal-wave dispersion: (7.136) → (7.137) → $\omega=N\cos\theta$ (7.139) — ★★, 8 steps, in C15 (notebook · `internal_wave_beams`)
- **Goal.** Find the frequency of a plane internal wave and show that it depends only on the direction of K.
- **Start.** $\frac{\partial^2}{\partial t^2}\nabla^2w+N^2\nabla_H^2w=0$ (7.134) with constant N.
- **Plan.** (1) Note the anisotropy. (2) Insert a plane wave. (3) Solve and express through the angle of K.
- **Tools.** complex amplitudes (P176) · angle of a vector (gloss).
- **Assumptions.** Constant N (step 2); l = 0 by rotating axes (step 6).
- **Steps.**
  1. *did:* Note the anisotropy · *tex:* $\omega=\omega(k,l,m)=\omega(\mathbf K)$ (7.135) · *why:* Gravity singles out the vertical: (7.134)
     treats ∇_H² and ∂²/∂z² differently, so ω may depend on the direction of K, not just its length. · *plain:* frequency
     may depend on direction.
  2. *did:* Try a plane wave · *tex:* $w=w_0e^{i(kx+ly+mz-\omega t)}=w_0e^{i(\mathbf K\cdot\mathbf x-\omega t)}$ (7.136) · *why:* Constant N means
     constant coefficients, so exponentials solve (7.134) (P176). · *plain:* one Fourier mode.
  3. *did:* Replace the derivatives · *tex:* $\frac{\partial^2}{\partial t^2}\to-\omega^2,\quad\nabla^2\to-K^2,\quad\nabla_H^2\to-(k^2+l^2)$ · *why:* Each
     ∂/∂x brings ik, ∂/∂t brings −iω, and (ik)² = −k². · *plain:* derivatives become multiplications.
  4. *did:* Substitute · *tex:* $\omega^2K^2w-N^2(k^2+l^2)w=0\ \Rightarrow\ \omega^2K^2=N^2(k^2+l^2)$ · *why:* (−ω²)(−K²) = ω²K²; divide by
     w ≠ 0. · *plain:* an algebraic relation. · *live (E9):* "ω = N|k|/K = …".
  5. *did:* Solve for ω² · *tex:* $\omega^2=\frac{k^2+l^2}{k^2+l^2+m^2}N^2$ (7.137) · *why:* K² = k² + l² + m² (7.6). · *plain:* the
     internal-wave dispersion relation.
  6. *did:* Work in the x–z plane · *tex:* $l=0:\ \omega=\frac{\lvert k\rvert N}{\sqrt{k^2+m^2}}=\frac{\lvert k\rvert N}K$ (7.138) · *why:* The medium is the
     same in every horizontal direction, so we point x along the horizontal part of K; positive root. The book writes kN/K
     assuming k > 0. · *plain:* frequency from the components.
  7. *did:* Introduce the angle of K · *tex:* $\cos\theta=\frac{\lvert k\rvert}K,\quad\theta=\text{angle of }\mathbf K\text{ above the horizontal}$ · *why:*
     Geometry of K = (k, m) (gloss: angle of a vector, `np.arccos`). · *plain:* the direction of the wave. · *set (E9):*
     ω/N = 0.3, then 0.9.
  8. *did:* Conclude · *tex:* $\omega=N\cos\theta$ (7.139), $0<\omega<N$ · *why:* Substitute step 7 into step 6; cos θ lies between 0
     and 1. · *plain:* frequency set by direction only; N is the maximum.
- **Result.** $\omega^2=\frac{k^2+l^2}{K^2}N^2$ (7.137), $\omega=N\cos\theta$ (7.139).
- **Check.** θ = 0 (m = 0): ω = N, the Ch. 1 parcel (N145 parity) ✓. N = 0.01 rad/s, θ = 45°: T = 14.8 min ✓. Doubling K at
  fixed direction leaves ω unchanged ✓ (C15 code: three equal values).
- **What it means.** An internal wave's frequency is its direction: a source at fixed ω radiates only at cos θ = ω/N —
  the beams of C16.
- **Traps.** Using k instead of |k| for a wave going to −x; mixing θ (K from the horizontal) with the beam angle (from the
  vertical: the same number); expecting ω to grow with K as for surface waves.

### D35 · $\mathbf K\cdot\mathbf u=0$ (7.141): internal waves are transverse — ★, 5 steps, in C16 (notebook · `internal_wave_beams`)
- **Goal.** Show that in an incompressible plane wave the water moves along the crests.
- **Start.** $u=u_0e^{i(kx+ly+mz-\omega t)}$ (7.140), and the same for v, w; (4.10).
- **Plan.** Differentiate the plane wave, insert into continuity, read the dot product.
- **Tools.** complex amplitudes (P176) · dot product.
- **Assumptions.** Plane wave; incompressible.
- **Steps.**
  1. *did:* Write the velocity as a plane wave · *tex:* $u=u_0e^{i\Theta},\ v=v_0e^{i\Theta},\ w=w_0e^{i\Theta},\ \Theta=kx+ly+mz-\omega t$ (7.140) ·
     *why:* In a linear plane wave every field shares the same phase factor. · *plain:* all three components ride one wave.
  2. *did:* Differentiate · *tex:* $\frac{\partial u}{\partial x}=iku,\quad\frac{\partial v}{\partial y}=ilv,\quad\frac{\partial w}{\partial z}=imw$ · *why:* The x-derivative
     of e^{iΘ} is ik e^{iΘ}, and likewise for y and z. · *plain:* derivatives become i × wavenumber.
  3. *did:* Insert into continuity · *tex:* $i(ku+lv+mw)=0$ · *why:* (4.10) ∂u/∂x + ∂v/∂y + ∂w/∂z = 0. · *plain:* incompressibility
     for a plane wave.
  4. *did:* Recognise the dot product · *tex:* $\mathbf K\cdot\mathbf u=0$ (7.141) · *why:* Divide by i ≠ 0; ku + lv + mw = K·u. · *plain:* the
     motion is perpendicular to K — along the crests. · *watch (E9):* "the teal arrows lie along the stripes".
  5. *did:* Name the wave type · *tex:* $\text{transverse (shear) wave}$ · *why:* Only continuity and the plane-wave form were used, so
     every incompressible plane wave is transverse; surface waves are not plane in z (they decay), so they are not. ·
     *plain:* internal waves are transverse.
- **Result.** $\mathbf K\cdot\mathbf u=0$ (7.141).
- **Check.** `internal_wave_fields` returns K·u = 0 to 1e-15 ✓. At θ = 0 (K horizontal) u = 0 and w ≠ 0: vertical bobbing ✓.
- **What it means.** Particles slide along the phase lines; combined with D36, the energy also runs along them.
- **Traps.** Thinking surface waves are transverse too (their orbits are loops, not lines); forgetting this needs a plane
  wave in all three directions.

### D36 · $\mathbf c$, $\mathbf c_g$ (7.143)–(7.145) and $\mathbf c_g\cdot\mathbf c=0$ (7.146) — ★★, 9 steps, in C16 (notebook · `internal_wave_beams`)
- **Goal.** Compute the phase and group velocities of an internal wave and show they are perpendicular.
- **Start.** $\omega=\frac{kN}K$ (7.138) (k > 0), $\mathbf c=(\omega/K)\mathbf e_K$ (7.8), $\mathbf c_g=\mathbf e_x\frac{\partial\omega}{\partial k}+\mathbf e_y\frac{\partial\omega}
  {\partial l}+\mathbf e_z\frac{\partial\omega}{\partial m}$ (7.143).
- **Plan.** (1) Differentiate ω(k, m). (2) Assemble c_g and c. (3) Handle k < 0. (4) Dot product and geometry.
- **Tools.** gradient in wavenumber space (gloss; P25) · quotient rule · D34.
- **Assumptions.** Constant N; x–z plane.
- **Steps.**
  1. *did:* Write c_g as a gradient · *tex:* $\mathbf c_g=\mathbf e_x\frac{\partial\omega}{\partial k}+\mathbf e_y\frac{\partial\omega}{\partial l}+\mathbf e_z\frac{\partial
     \omega}{\partial m}$ (7.143) · *why:* In several dimensions the group velocity is the gradient of ω in wavenumber space
     (N71, C09). · *plain:* energy velocity from the slopes of ω in each direction.
  2. *did:* Write ω with k > 0 · *tex:* $\omega=NkK^{-1},\quad K=(k^2+m^2)^{1/2}$ · *why:* (7.138) in the x–z plane; k < 0 comes in step
     7. · *plain:* a function of two variables.
  3. *did:* Differentiate in k · *tex:* $\frac{\partial\omega}{\partial k}=\frac NK-\frac{Nk^2}{K^3}=\frac{N(K^2-k^2)}{K^3}=\frac{Nm^2}{K^3}$ · *why:* Product
     rule with ∂K/∂k = k/K. · *plain:* horizontal part of c_g.
  4. *did:* Differentiate in m · *tex:* $\frac{\partial\omega}{\partial m}=-\frac{Nkm}{K^3}$ · *why:* Only K depends on m, ∂K/∂m = m/K. · *plain:*
     vertical part of c_g.
  5. *did:* Assemble c_g · *tex:* $\mathbf c_g=\frac{Nm}{K^3}(m\mathbf e_x-k\mathbf e_z)$ (7.145) · *why:* Collect the two partial derivatives of steps 3 and 4 into one vector and factor out Nm/K³. · *plain:* the group velocity.
     · *live (E9):* "∂ω/∂k = …, ∂ω/∂m = …".
  6. *did:* Write c · *tex:* $\mathbf c=\frac\omega{K^2}(k\mathbf e_x+m\mathbf e_z)$ (7.144) · *why:* The phase velocity (7.8) c = (ω/K)e_K points along K; with e_K = K/K it becomes (ω/K²)(k, m). · *plain:* the phase
     velocity points along K.
  7. *did:* Allow k < 0 · *tex:* $\omega=\frac{N\lvert k\rvert}K\ \Rightarrow\ \frac{\partial\omega}{\partial k}=\mathrm{sgn}(k)\frac{Nm^2}{K^3},\quad\frac{\partial\omega}{\partial m}=
     -\frac{N\lvert k\rvert m}{K^3}$ · *why:* The printed (7.138), (7.145) assume k > 0; for K pointing up-left (Fig. 7.29) use |k|.
     The horizontal parts of c and c_g keep a common sign. · *plain:* the sign-safe form the code uses.
  8. *did:* Take the dot product · *tex:* $\mathbf c\cdot\mathbf c_g=\frac\omega{K^2}\frac{Nm}{K^3}(km-mk)=0$ (7.146) · *why:* Multiply
     components and add; the two terms cancel exactly (also in the sign-safe form). · *plain:* phase and group velocities
     are perpendicular. · *watch (E9):* "the right-angle mark".
  9. *did:* Read the geometry · *tex:* $c_z=\frac{\omega m}{K^2}=\frac{Nkm}{K^3}=-c_{g,z}$ · *why:* Use ω = Nk/K in c_z and compare with
     step 5's vertical part. · *plain:* vertical parts equal and opposite: phase up ⇔ energy down; horizontal parts the same
     way.
- **Result.** $\mathbf c=\frac\omega{K^2}(k\mathbf e_x+m\mathbf e_z)$, $\mathbf c_g=\frac{Nm}{K^3}(m\mathbf e_x-k\mathbf e_z)$ (7.144, 7.145),
  $\mathbf c_g\cdot\mathbf c=0$ (7.146).
- **Check.** N = 1 rad/s, K = 1 rad/m at 45°: c = (0.5, 0.5), c_g = (0.5, −0.5) m/s, dot 0 ✓; the central-difference
  gradient agrees to 1e-7 ✓. |c_g| = N sin θ/K → 0 as θ → 0 (ω → N) ✓.
- **What it means.** Energy leaves along the crests, at right angles to the phase: the St Andrew's cross, tidal beams
  going down while their phase goes up.
- **Traps.** Using the printed (7.145) for k < 0 (it points the wrong way); thinking c and c_g are parallel as for
  surface waves.

### D37 · Polarization (7.153), mean flux (7.158) and $\mathbf F=\mathbf c_gE$ (7.159) — ★★, 12 steps, in C16 (notebook)
- **Goal.** Express every field of an internal plane wave through w, compute the mean energy flux and show it equals group
  velocity × energy.
- **Start.** $[u,w,p',\rho']=[\hat u,\hat w,\hat p,\hat\rho]e^{i(kx+mz-\omega t)}$ with (7.128), (7.131), (7.132).
- **Plan.** (1) Polarization relations from three equations. (2) Mean flux with real parts. (3) Compare with c_g E.
- **Tools.** complex amplitudes (P176) · mean of a product of real parts (P178) · complex conjugate (P81) · (7.145),
  (7.157).
- **Assumptions.** Plane wave in the x–z plane, k > 0, real ŵ.
- **Steps.**
  1. *did:* Assume plane-wave fields · *tex:* $[u,w,p',\rho']=[\hat u,\hat w,\hat p,\hat\rho]e^{i\theta},\ \theta=kx+mz-\omega t,\ \hat w\ \text{real}$ · *why:*
     The equations are linear with constant N; choosing the phase origin makes ŵ real. · *plain:* every field is a complex
     amplitude times one wave.
  2. *did:* Get p̂ from (7.132) · *tex:* $-\frac{k^2}{\rho_0}\hat p=(im)(-i\omega)\hat w=m\omega\hat w\ \Rightarrow\ \hat p=-\frac{\omega m\rho_0}{k^2}\hat w$ · *why:*
     ∇_H² → −k² (l = 0) and ∂²/∂z∂t → (im)(−iω) = mω. · *plain:* the pressure amplitude.
  3. *did:* Get ρ̂ from (7.131) · *tex:* $-i\omega\hat\rho=\frac{N^2\rho_0}g\hat w\ \Rightarrow\ \hat\rho=\frac{iN^2\rho_0}{\omega g}\hat w$ · *why:* ∂/∂t → −iω and
     1/(−i) = i; the factor i means ρ′ is 90° out of phase with w. · *plain:* the density amplitude.
  4. *did:* Get û from (7.128) · *tex:* $-i\omega\hat u=-\frac{ik}{\rho_0}\hat p\ \Rightarrow\ \hat u=-\frac mk\hat w$ · *why:* ∂/∂t → −iω, ∂/∂x → ik;
     insert p̂. It agrees with K·u = 0 (D35): kû + mŵ = 0. · *plain:* the horizontal velocity amplitude.
  5. *did:* Collect the polarization relations · *tex:* $p'=-\frac{\omega m\rho_0}{k^2}\hat we^{i\theta},\ \rho'=\frac{iN^2\rho_0}{\omega g}\hat we^{i\theta},\
     u=-\frac mk\hat we^{i\theta}$ (7.153) · *why:* Collect the three amplitudes of steps 2–4, each multiplied by the common wave factor e^{iθ}; ŵ is the only free constant. · *plain:* everything in terms of w.
  6. *did:* Write the mean flux · *tex:* $\mathbf F=\overline{p'\mathbf u}=\mathbf e_x\overline{p'u}+\mathbf e_z\overline{p'w}$ · *why:* The energy flux
     is the pressure work p′u (7.147); average over a period. · *plain:* power per unit area.
  7. *did:* Take real parts before multiplying · *tex:* $\big\langle\mathrm{Re}(Pe^{i\theta})\,\mathrm{Re}(Qe^{i\theta})\big\rangle=\tfrac12\mathrm{Re}(PQ^*)$ ·
     *why:* The real part of a product is not the product of real parts; this rule averages correctly (P178). · *plain:* the
     right way to average a product.
  8. *did:* Compute the horizontal flux · *tex:* $\overline{p'u}=\tfrac12\mathrm{Re}(\hat p\hat u^*)=\frac{\rho_0\omega m^2\hat w^2}{2k^3}$ · *why:* p̂ û* =
     (−ωmρ₀ŵ/k²)(−mŵ/k), both real. · *plain:* horizontal energy flux.
  9. *did:* Compute the vertical flux · *tex:* $\overline{p'w}=\tfrac12\mathrm{Re}(\hat p\hat w)=-\frac{\rho_0\omega m\hat w^2}{2k^2}$ · *why:* Because ŵ is real, ŵ* = ŵ in ½Re(p̂ŵ*); the product p̂ŵ is real and negative for m > 0. ·
     *plain:* vertical energy flux (downward for m > 0).
  10. *did:* Assemble F · *tex:* $\mathbf F=\frac{\rho_0\omega m\hat w^2}{2k^2}\Big(\mathbf e_x\frac mk-\mathbf e_z\Big)$ (7.158) · *why:* Factor
      ρ₀ωmŵ²/(2k²) out of steps 8–9. · *plain:* the energy flux vector.
  11. *did:* Form c_g E · *tex:* $\mathbf c_gE=\frac{Nm}{K^3}(m\mathbf e_x-k\mathbf e_z)\frac{\rho_0}2\Big(\frac{m^2}{k^2}+1\Big)\hat w^2=\frac{\rho_0Nm\hat w^2}
      {2Kk^2}(m\mathbf e_x-k\mathbf e_z)$ · *why:* (7.145) times (7.157) E = ½ρ₀(m²/k² + 1)ŵ², and m²/k² + 1 = K²/k². · *plain:*
      group velocity times energy.
  12. *did:* Use the dispersion relation · *tex:* $\frac NK=\frac\omega k\ \Rightarrow\ \mathbf c_gE=\frac{\rho_0\omega m\hat w^2}{2k^2}\Big(\mathbf e_x\frac mk-\mathbf e_z\Big)=
      \mathbf F$ (7.159) · *why:* (7.138) ω = kN/K; take 1/k out of (m, −k). · *plain:* energy flux = group velocity × energy.
- **Result.** (7.153), (7.158) and $\mathbf F=\mathbf c_gE$ (7.159).
- **Check.** Units: ρ₀ωŵ²/k [kg/m³ · 1/s · m²/s² · m] = W/m² ✓. `internal_wave_energy` asserts F = c_g E and E_k = E_p
  ((7.154)–(7.156)) to 1e-12 ✓.
- **What it means.** The same law as (7.71) for surface waves: wave energy is transported at the group velocity —
  here at right angles to the phase. Per unit volume and area here, per unit area and crest length there.
- **Traps.** Multiplying complex fields before taking real parts; using ŵ instead of ŵ* in ⟨p′w⟩ (harmless here, ŵ real,
  but wrong in general); forgetting ρ′ is 90° out of phase (it contributes nothing to a flux with w).
