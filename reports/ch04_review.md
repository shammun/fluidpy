# ch04 — Derivation review (Phase 6)

Independent, read-only line-by-line review of `fluidpy/core/{conservation,constitutive,navier_stokes,streamfunction,rotating,curvilinear,bernoulli,interfaces,similarity,_stencil}.py`,
`fluidpy/core/thermo.py` additions and `fluidpy/ch04_conservation_laws.py` against the rendered book pages
(p127, p134, p144, p152–154, p157, p162–164, p168–170, p173–176). Suite at review time: `tests/test_ch04.py` 147 passed.

**Verdict:** no sign or factor error in the transcribed equations; two paths return wrong numbers for inputs the code
accepts and the tests do not reach. Both fixed (see Resolution).

## Must fix
1. **`cv_scenario("jet")` splits the oblique jet 50/50 for every θ** (`ch04_conservation_laws.py` ~480). For an ideal jet
   on a shear-free plate, along-plate momentum with equal sheet speeds forces $Q_{1,2} = Q(1 \pm \cos\theta)/2$. At θ = 45°
   the code gave 5.0 / 5.0 kg/s instead of 8.54 / 1.46. The E1 explainer's θ slider shows these fluxes. Fix: split by
   $(1\pm\cos\theta)/2$ and add an along-plate momentum residual.
2. **`mean_pressure` / `thermodynamic_pressure_from_stress` wrong for 2×2 input** (`constitutive.py` ~186, ~195). (4.32)
   $\tau_{ii} = -3p + (2\mu + 3\lambda)\nabla\cdot\mathbf{u}$ and (4.33) $\bar p = -\tau_{ii}/3$ need $\tau_{33}$, which a 2×2
   τ does not carry. The probe G = diag(1, 0), p = 100, μ = 1, μ_v = 0.5 gave $p - \bar p = 0.833$, but (4.34) requires
   $\mu_v \nabla\cdot\mathbf{u} = 0.5$. Fix: raise for d ≠ 3 (or accept τ₃₃) and add a 2×2 test.

## Should fix
1. `rotating.py`: a scalar `dU_dt` is silently treated as the z-component; `apparent_body_forces` crashes on a 2-vector u′.
2. `boussinesq_validity` compares M with 0.1. (4.110) makes the departure from ∇·u = 0 of size M²; §4.2 accepts
   M < 0.3. Compare M² (or use the M < 0.3 rule) so it agrees with `is_incompressible_regime`.
3. Docstring derivation IDs use the analysis numbering: rotating.py D24 → D15, D27 → D18; navier_stokes.py D31/D32/D33 → D20/D21/D22, D58 → a-D58.
4. Stale contract numbers in docstrings: bore_speed 3.3660 → 3.3661 m/s; capillary_length 2.7266 → 2.7269 mm.
5. `plane_poiseuille(y, G=…)` uses G for −dp/dx while `ch04.G` is gravity; state the clash in the docstring (signature kept stable).
6. `material_mass` for dim > 1 is a closed-form identity: label V1, not V4.
7. `core/conservation.py`: a 2-D control volume silently drops the default g = (0, 0, −g); document it.

## Verified
- Frames (4.42)–(4.45): the +2Ω×u′ and Ω×(Ω×x′) terms of (4.43), the bracket signs of (4.45), the projectile ODE (right-hand deflection in the NH), effective gravity, and the Example 4.5 rotation terms.
- Stream functions: plane $-\mathbf{e}_z\times\nabla\psi$ and axisymmetric $\rho u_R = -R^{-1}\psi_z$, $\rho u_z = R^{-1}\psi_R$ match p127.
- Cauchy: first index contracted; (4.31)/(4.37)/(4.59) consistent through λ = μ_v − ⅔μ; cube-spin torque.
- Energy: the (4.48) heat sign; (4.62)/(4.63) entropy production; the (4.112)–(4.116) coefficients; Couette ΔT_max = μU²/8k; the closed-form transient b_n.
- Bernoulli: (4.74) gauge sign; Rankine and Lamb–Oseen pressures; U-tube; accelerating sphere (added mass −½ρVU̇); pitot correction.
- Examples 4.3 (bore, any b), 4.4 (rocket face flux ṁ(b − V_e)), 4.6 (sprinkler), 4.7 (meniscus), 4.8 (wave drag λ⁻³ scaling).
- Surface tension (4.97)–(4.98): rim vector direction and the jump on the concave side.
- Boussinesq (4.86), (4.88), H_c = c²/g; similarity coefficients (4.101), (4.110).
- Every stencil in `_stencil.py` is an explicit second-order central difference.

## Resolution
Must-fix 1–2 and should-fix 1–7 sent to `concept-implementer` (fluidpy only, public signatures stable); new tests
for the along-plate jet momentum and the 2×2 pressure guard added by `math-verifier`. Re-run: tests/test_ch04.py 151 passed, full suite 478 passed; each pre-fix behaviour re-planted fails its new test.

**Review status: PASS** (all must-fix resolved).
