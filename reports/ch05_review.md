# ch05 — Derivation review (Phase 6)

Independent, read-only line-by-line review of `fluidpy/core/vorticity.py`, `fluidpy/core/biot_savart.py`,
`fluidpy/core/vortices.py` (ch05 additions), `fluidpy/core/integral_theorems.py` (`curl_theorem_box`),
`fluidpy/ch05_vorticity_dynamics.py` and `scripts/ch05_*.py` against rendered pages p201–p222. Suite at review time:
`tests/test_ch05.py` 133 passed.

**Verdict:** no must-fix. Every equation and every scenario component checked against the pages, or re-derived, is
correct. Seven should-fix items (docstrings, two guard rails used by explainers, one label, one print).

## Must fix
None.

## Should fix
1. `uniform_strain_vorticity` docstring: the `"shear_tilt"` background has vorticity (0, −s, 0), which is neither zero
   nor parallel. The result $e^{Gt}\omega_0$ is right as a passive material line element; say so.
2. `circle_image_system`: a vortex exactly at the circle centre gives an image at 0/0 = NaN, which spreads to every
   vortex through `point_vortex_rhs`. Drop the image (it goes to infinity) when $|x - c| = 0$. E8 lets users drag vortices.
3. `ring_dynamics` with `wall_z`: the thin-core model runs past its range (the ring ends 0.0046 m from the wall with
   a = 0.0115 m). Add a terminal event when the gap falls below a few core radii, or a regime warning.
4. `gaussian_tube_fields`: say that the flat-ended tube is one segment of a longer tube (V′ in Fig. 5.8) and is not a
   valid field on its own ($\nabla\cdot\omega$ has δ-sheets at the ends). C01 teaches that tubes cannot end.
5. Label honesty: Kelvin's ring speed and `vortex_pair` are published closed forms, so they are V1 form cross-checks,
   not V5 benchmarks. Make the test, the docstrings and the verification report consistent.
6. `ring_self_velocity(core="gaussian")`: state that a is the $e^{-r^2/a^2}$ radius, $a = \sqrt{4\nu t}$, for Saffman's 0.558.
7. `scripts/ch05_kelvin_baroclinic.py:148` prints "=" between 0.242259 and 0.242235. Use "≈" and name the
   mean-density hydrostatic model.

## Verified
- (5.30) budget term by term against p212: stretching acts on absolute vorticity; the planetary term is 2(Ω·∇)u; the baroclinic term is ∇ρ × ∇p/ρ² in that order. Also (5.28), (5.25) with its Lamb term, (5.20), and the (5.31) split (e_n = −Frenet N).
- (5.14) sign: +1/(4π), from ∇²u = −∇×ω and Exercise 5.9's Green's function; the book's compensating slip is confirmed, and (5.16)/(5.17) match p209.
- Sheet sign: dΓ = (u₂ − u₁)ds, counterclockwise; agrees with Exercise 5.20(b). The roll-up cot kernel was re-derived.
- Point vortices and images component by component: opposite pair, wall image (Γ/4πh = 0.397887), knife pair, the circle images of Exercise 5.14, and the channel closed form (Γ/4H)cot(πh/H), which matches the series (0.3440955).
- Rings: `ring_ring_velocity` against a 4000-gon Biot–Savart (~1e-6); m = k²; image ring −Γ; leap-frog radii behave as p218 describes; Kelvin constant ¼.
- Hill's vortex (after the F2 fix), Burgers vortex with its pressure, §5.1 vortices (σ_rθ, torque −2μΓ, dissipation budget), and Exercises 5.1 and 5.2.
- Baroclinic mechanics: D06 disc torque; lock exchange 2Δρg/((ρ₁+ρ₂)δ), counterclockwise with heavy fluid on the left.
- Kelvin scenarios: the cellular pressure satisfies Euler; rotating-loop Γ = 2ΩA₀(1 − e^{−αt}); (5.33) Γ_a = Γ + 2Ω·A; the column (ζ + f)/h.

## Resolution
Should-fix 1–4, 6 and 7 go to `concept-implementer` (fluidpy and scripts; signatures stable). Item 5, the test label, goes to `math-verifier`. Re-run: tests/test_ch05.py 135 passed, full suite 613 passed. The ring and pair tests are relabelled V1. The circle-centre image and the ring validity stop (9.4048 s, confirmed independently) are tested; `terminate` stays False by default and plots use the `valid` mask.

**Review status: PASS** (no must-fix; should-fix 1–7 resolved).
