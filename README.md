# fluidpy — *Fluid Mechanics* (Kundu, Cohen & Dowling, 5th ed.), learned through Python and interactive explainers

Each chapter of the book becomes one learning package:

* a **teaching notebook** — plain-words explanations, step-by-step mathematics, tiny worked examples, heavily commented
  Python that calls tested functions, figures, animations and interactive plotly figures;
* up to **five interactive explainers** — full-window, phone-friendly HTML pages with a guided walkthrough, sliders,
  the equations with your numbers substituted, and check-yourself questions;
* the same notebook **ready for Google Colab**, and a **web page** with every explainer shown full-window.

Read online: **https://shammun.github.io/fluidpy/** · all explainers: https://shammun.github.io/fluidpy/viz/

## Chapters

<!-- INDEX_TABLE_START -->
| Chapter | Read (HTML) | Notebook | Colab | Explainers |
|---|---|---|---|---|
| Ch. 1: Introduction | [View](https://shammun.github.io/fluidpy/notebooks/ch01_introduction.html) | [.ipynb](notebooks/ch01_introduction.ipynb) | [Open](https://colab.research.google.com/github/shammun/fluidpy/blob/main/notebooks/ch01_introduction_colab.ipynb) | [When does 'density at a point' make sense?](https://shammun.github.io/fluidpy/viz/ch01/continuum_averaging_volume.html), [How does the fluid learn that a plate moved?](https://shammun.github.io/fluidpy/viz/ch01/viscosity_momentum_diffusion.html), [Same two states: what depends on the path?](https://shammun.github.io/fluidpy/viz/ch01/heat_work_paths.html), [Push a parcel up: does it come back?](https://shammun.github.io/fluidpy/viz/ch01/parcel_stability.html), [Why can 7 pipe variables shrink to 4 numbers?](https://shammun.github.io/fluidpy/viz/ch01/buckingham_pi_machine.html) |
| Ch. 2: Cartesian Tensors | [View](https://shammun.github.io/fluidpy/notebooks/ch02_cartesian_tensors.html) | [.ipynb](notebooks/ch02_cartesian_tensors.ipynb) | [Open](https://colab.research.google.com/github/shammun/fluidpy/blob/main/notebooks/ch02_cartesian_tensors_colab.ipynb) | [Which rotates — the arrow or the ruler?](https://shammun.github.io/fluidpy/viz/ch02/rotation_of_axes.html), [Cut the point any way you like: what pushes on the cut?](https://shammun.github.io/fluidpy/viz/ch02/cauchy_traction_principal_axes.html), [Is simple shear a rotation?](https://shammun.github.io/fluidpy/viz/ch02/strain_vs_rotation_split.html), [What leaks out of a box?](https://shammun.github.io/fluidpy/viz/ch02/gauss_flux_box.html), [How much does the flow go round a loop?](https://shammun.github.io/fluidpy/viz/ch02/stokes_circulation_loop.html) |
| Ch. 3: Kinematics | [View](https://shammun.github.io/fluidpy/notebooks/ch03_kinematics.html) | [.ipynb](notebooks/ch03_kinematics.ipynb) | [Open](https://colab.research.google.com/github/shammun/fluidpy/blob/main/notebooks/ch03_kinematics_colab.ipynb) | [Three lines through one point — why do they disagree?](https://shammun.github.io/fluidpy/viz/ch03/flow_lines_unsteady.html), [Why does the station warm while the air does not?](https://shammun.github.io/fluidpy/viz/ch03/material_derivative_probe.html), [Steady or not — does the acceleration care?](https://shammun.github.io/fluidpy/viz/ch03/galilean_frames_cylinder.html), [What does each number in S measure?](https://shammun.github.io/fluidpy/viz/ch03/fluid_element_deformation.html), [Can a straight flow make a fluid element spin?](https://shammun.github.io/fluidpy/viz/ch03/spin_and_principal_axes.html), [Going round in circles or spinning?](https://shammun.github.io/fluidpy/viz/ch03/vortex_paddle_wheels.html), [What changes inside a moving box?](https://shammun.github.io/fluidpy/viz/ch03/reynolds_transport_cv.html) |
| Ch. 4: Conservation Laws | [View](https://shammun.github.io/fluidpy/notebooks/ch04_conservation_laws.html) | [.ipynb](notebooks/ch04_conservation_laws.ipynb) | [Open](https://colab.research.google.com/github/shammun/fluidpy/blob/main/notebooks/ch04_conservation_laws_colab.ipynb) | [Weigh a force by counting what flows through a box](https://shammun.github.io/fluidpy/viz/ch04/control_volume_budgets.html), [Can one number field hold a whole 2-D flow?](https://shammun.github.io/fluidpy/viz/ch04/stream_function_spacing.html), [How does a fluid decide its stress?](https://shammun.github.io/fluidpy/viz/ch04/newtonian_stress_lab.html), [Which terms of Navier–Stokes are awake here?](https://shammun.github.io/fluidpy/viz/ch04/navier_stokes_term_balance.html), [Why does a straight throw curve on a merry-go-round?](https://shammun.github.io/fluidpy/viz/ch04/rotating_frame_coriolis.html), [Bernoulli is constant along what, exactly?](https://shammun.github.io/fluidpy/viz/ch04/which_bernoulli.html), [Where does the energy go when viscosity stops a flow?](https://shammun.github.io/fluidpy/viz/ch04/viscous_dissipation_heating.html), [Density hardly changes — so why does it drive the flow?](https://shammun.github.io/fluidpy/viz/ch04/boussinesq_buoyancy.html), [When does a model behave like the real thing?](https://shammun.github.io/fluidpy/viz/ch04/dynamic_similarity_models.html) |
| Ch. 5: Vorticity Dynamics | — | — | — | — |
| Ch. 6: Ideal Flow | — | — | — | — |
| Ch. 7: Gravity Waves | — | — | — | — |
| Ch. 8: Laminar Flow | — | — | — | — |
| Ch. 9: Boundary Layers and Related Topics | — | — | — | — |
| Ch. 10: Computational Fluid Dynamics | — | — | — | — |
| Ch. 11: Instability | — | — | — | — |
| Ch. 12: Turbulence | — | — | — | — |
| Ch. 13: Geophysical Fluid Dynamics | — | — | — | — |
| Ch. 14: Aerodynamics | — | — | — | — |
| Ch. 15: Compressible Flow | — | — | — | — |
| Ch. 16: Introduction to Biofluid Mechanics | — | — | — | — |
<!-- INDEX_TABLE_END -->

## How correctness is established
The book ships no code. Every function in `fluidpy/` is written from the book's equations and then made to prove
itself — closed-form solutions, symbolic re-derivation with sympy, measured order of accuracy, conservation checks,
cited published benchmarks and asymptotic limits (`tests/`, `reports/chNN_verification.md`). Each explainer's JavaScript
physics is checked against the Python functions it mirrors.

## Run locally
```bash
python -m venv .venv
.venv/Scripts/pip install -r requirements.txt           # Windows;  .venv/bin/pip on Linux/macOS
.venv/Scripts/python -m pytest -q                        # the evidence
jupyter lab notebooks/                                   # open any chapter notebook (run its setup cell first)
```

## Layout
| Path | What |
|---|---|
| `fluidpy/chNN_<slug>.py`, `fluidpy/core/` | chapter physics; shared primitives and notebook machinery |
| `notebooks/chNN_<slug>.ipynb` · `_colab.ipynb` · `.html` | executed notebook · Colab twin · published page |
| `viz/chNN/<slug>.html` | interactive explainers (self-contained; engine in `assets/viz_lib.js`) |
| `tests/`, `reference/`, `reports/` | evidence, cited benchmark data, verification / review / explainer reports |
| `analysis/`, `knowledge/` | what each chapter contains, what we chose to teach and why, and what later chapters reuse |
| `tools/` | chapter split, page rendering, publishing, browser audits, public-repo guard |

## About the source material
This project follows Pijush K. Kundu, Ira M. Cohen and David R. Dowling, *Fluid Mechanics*, 5th edition (Academic
Press, 2012). The book's text, figures, tables and exercises are **not** included: explanations are written in our own
words, figures are generated by our code, and equations are cited by number so you can follow along with your copy.

## License
Code, notebooks, explainers and documentation: MIT (see `LICENSE`). The license does not cover the book.
