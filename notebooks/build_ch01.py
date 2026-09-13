"""Build the Chapter 1 teaching notebook: ``notebooks/ch01_introduction.ipynb``.

Source of truth: ``analysis/ch01_design.md`` Part A (storyboard, one nbkit call per row), Part E (prerequisite ledger →
primers and recaps), Part F (the 12 derivations, copied step by step) and ``analysis/ch01_curation.md`` (IDs, depths,
section coverage). Physics lives in ``fluidpy.ch01_introduction`` (imported as ``ch01``); cells only call it.

Labels shown to the reader: CORE blocks carry their id (`C06`), notes their id (`C07`), primers their number (`P26`),
so cross-references such as "(P26, C20)" in the derivations can be followed.

Run:  .venv/Scripts/python.exe notebooks/build_ch01.py
"""
from __future__ import annotations

import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "tools"))
from nbkit import ChapterNotebook  # noqa: E402

nb = ChapterNotebook("ch01")


# ---------------------------------------------------------------------------------------------------------------------
# small helpers that keep the labels consistent
# ---------------------------------------------------------------------------------------------------------------------
def core(cid: str, title: str, question: str) -> None:
    """A CORE block heading with its id."""
    nb.core(cid, f"{title} `{cid}`", question=question)


def P(pid: str, term: str, text: str, code: str | None = None) -> None:
    """A 📎 primer; ``term`` is exactly the Part E concept text (the ledger check matches it)."""
    nb.primer(term, f"`{pid}` · {text}", code=code)


def note(nid: str, title: str, text: str, equation: str | None = None, ref: str | None = None) -> None:
    """A B/C note with its curation id."""
    nb.note(f"**{title}** `{nid}` — {text}", equation=equation, ref=ref)


def D(key: str, title: str, *, ref: str, goal: str, assumptions: str, start, plan, uses, steps, result,
      interpret: str, check: str, traps: str, check_src: str | None = None) -> None:
    """A Part F derivation: goal + assumptions, plan, tools, one move per step, result, meaning, check, traps."""
    nb.derivation(key, title, ref=ref, goal=f"{goal}\n\n**Assumptions.** {assumptions}", start=start, plan=plan,
                  uses=uses, steps=[dict(did=a, tex=b, why=c, plain=d) for a, b, c, d in steps], result=result,
                  interpret=interpret, check=check, check_src=check_src)
    nb.md(f"> ⚠️ **Common confusion (traps in this derivation):** {traps}")


def see_read_change(see: str, read: str, change: str) -> None:
    nb.figure_notes(see, read, change)


# =====================================================================================================================
# A.0 front matter
# =====================================================================================================================
nb.title(
    big_idea=(
        "A fluid is matter that keeps deforming under any shear stress, and we describe it as a **continuum**: density, "
        "velocity, pressure and temperature have a value at every point. Chapter 1 builds the vocabulary that the rest of "
        "the book stands on — molecules averaged into fields (§1.4), molecules carrying momentum and heat down gradients "
        "(§1.5), pressure in a fluid at rest (§1.7), the thermodynamics of a small fluid particle (§1.8–1.9), when a "
        "stratified column of air or water is stable (§1.10), and how dimensional analysis shrinks a problem to a few "
        "numbers (§1.11)."
    ),
    roadmap=[
        "The continuum hypothesis — a value at a point is the plateau of an average (C06)",
        "Newton's law of viscosity — shear stress is a flux of momentum (C12)",
        "The hydrostatic law — pressure grows downward by the weight above (C20)",
        "The first law — heat plus work changes internal energy (C25)",
        "Entropy and the Gibbs relations (C35)",
        "The speed of sound — stiffness at constant entropy (C36)",
        "The perfect-gas law p = ρRT (C40)",
        "The isentropic law p/ρ^γ = const (C45)",
        "The displaced parcel (C50)",
        "The Brunt–Väisälä frequency N² (C51)",
        "The adiabatic lapse rate, in both sign conventions (C54)",
        "Potential temperature θ (C55)",
        "Dimensional homogeneity (C64)",
        "The dimensional matrix (C67)",
        "Buckingham's Π theorem (C69)",
    ],
    prerequisites=[
        "school physics: force, energy, SI units",
        "calculus of one variable, and a little linear algebra",
        "Python with numpy — everything else (partial derivatives, Maxwell relations, null spaces, sympy, plotly …) "
        "is primed with a 📎 box right where it is first used",
    ],
)
nb.explainer_index([
    ("continuum_averaging_volume", "When does 'density at a point' make sense?", "a value at a point is the plateau of an average"),
    ("viscosity_momentum_diffusion", "How does the fluid learn that a plate moved?", "viscosity diffuses momentum; ν sets the clock"),
    ("heat_work_paths", "Same two states: what depends on the path?", "q and w depend on the route, Δe and Δs do not"),
    ("parcel_stability", "Push a parcel up: does it come back?", "stability compares two lapse rates — in both sign conventions"),
    ("buckingham_pi_machine", "Why can 7 pipe variables shrink to 4 numbers?", "Π groups are the null space of a matrix"),
])
nb.setup()

# =====================================================================================================================
# A.1 §1.1 Fluid Mechanics
# =====================================================================================================================
nb.section("1.1", "Fluid Mechanics", intro="""
**What is this section about?** What the subject covers and the three ways it makes progress — and a map of where the
ideas of this chapter are used later in the book.
""")
note("C01", "What fluid mechanics is, and its three routes", """
Fluid mechanics studies liquids and gases at rest and in motion — from blood in a capillary to the jet stream. It moves
forward along three routes that check one another: **analysis** (equations solved with pencil and paper — most of this
book), **computation** (numerical solutions, Ch. 10) and **experiment** (laboratory and field measurements). All three
start from the same picture of a fluid as a continuum, which the C06 block of §1.4 builds. The map below shows how the
chapters lean on each other.
""")
P("P01", "matplotlib figures", """
`fig, ax = plt.subplots()` makes a figure with one set of axes; `ax.plot(x, y)` draws a line; `ax.set_xlabel(...)`
labels an axis — always with its unit in square brackets. Every figure in this notebook is followed by the same three
reading notes: *What you see*, *How to read it*, *What would change if…*.
""", code="""
import matplotlib.pyplot as plt                # matplotlib's plotting interface, used for every static figure
import logging                                 # standard library: controls library log messages
logging.getLogger("matplotlib.font_manager").setLevel(logging.ERROR)   # hide harmless font-substitution notes
x = [0, 1, 2, 3]                               # a few x values [m]
y = [0, 1, 4, 9]                               # y = x² at those points [m²]
fig, ax = plt.subplots(figsize=(3.6, 2.4))     # one small figure with one set of axes
ax.plot(x, y, marker="o")                      # draw the points joined by straight lines
ax.set_xlabel("x [m]")                         # label the horizontal axis, with its unit
ax.set_ylabel("y [m²]")                        # label the vertical axis, with its unit
plt.show()                                     # display the figure below the cell
""")
nb.figure("""
from scripts.ch01_book_map import draw_book_map   # our drawing helper (pure matplotlib, no physics) from scripts/
fig, ax = plt.subplots(figsize=(10, 4.8))          # one wide figure for 16 chapter boxes
draw_book_map(ax, highlight=("ch13",))             # boxes = chapters, arrows = "builds on"; routes into Ch. 13 highlighted
import textwrap                                    # standard library: wraps long text onto several lines
from fluidpy.core.project import chapters          # the chapter list from book.yaml (titles)
full = {c["title"][:20] + "…": c["title"] for c in chapters()}   # the helper shortens long titles to 20 characters + "…"
for label in ax.texts:                             # every text drawn on the map
    if label.get_text() in full:                   # a shortened title …
        label.set_text(textwrap.fill(full[label.get_text()], 22))   # … gets its full title back, on two lines
ax.set_title("Geophysical fluid dynamics (Ch. 13) sits at the end of the longest chains of the book")  # the message
plt.show()                                         # display it
""", see="""
Sixteen chapter boxes arranged from foundations (left) to applications (right), with arrows meaning "this chapter
needs that one". The routes that end in Chapter 13 (geophysical fluid dynamics) are highlighted.
""", read="""
Follow the highlighted arrows backwards from Ch. 13: it rests directly on the conservation laws (Ch. 4), vorticity
(Ch. 5), waves (Ch. 7), laminar flows (Ch. 8), instability (Ch. 11) and turbulence (Ch. 12) — and all of them rest on
this chapter's vocabulary.
""", change="""
Call `draw_book_map(ax, highlight=("ch15",))` instead: the compressible-flow chapter lights up, and its roots run
straight back through Ch. 4 to the thermodynamics of §1.8–1.9.
""")

# =====================================================================================================================
# A.2 §1.2 Units of Measurement
# =====================================================================================================================
nb.section("1.2", "Units of Measurement", intro="""
**What is this section about?** SI units are school knowledge; we recap them once and let the `pint` library carry
units through calculations, so every later number has its unit attached.
""")
nb.recap("R01", "SI base quantities and derived units", r"""
Four base quantities are enough for this book: length (m), mass (kg), time (s) and temperature (K). Everything else is
built from them: newton N = kg m s⁻², pascal Pa = N m⁻², joule J = N m, watt W = J s⁻¹, hertz Hz = s⁻¹. So
$[\mathrm{Pa}] = \mathrm{kg\,m^{-1}\,s^{-2}}$. **A unit is a product of powers of base units** — exactly the idea that
§1.11 turns into dimensional analysis (C64).
""", where="pre-book, school physics; used again in C64")
P("P02", "pint quantities", """
`Q_(value, "unit")` attaches a unit to a number; `.to("Pa")` converts it; `.dimensionality` shows the powers of the base
dimensions. Temperatures in °C are *offset* units (0 °C is not "no temperature"), so formulas always use kelvin.
""", code="""
from fluidpy.core.units import Q_                          # Q_(value, "unit"): a number that carries its unit (pint)
print(Q_(1, "bar").to("Pa"))                               # 1 bar in pascals → 100000 pascal
w = Q_(1000, "kg/m**3") * Q_(9.81, "m/s**2") * Q_(10, "m") # ρ g z for 10 m of water: the units multiply as well
print(w.to("Pa"), "|", w.dimensionality)                   # 98100 pascal | [mass] / [length] / [time] ** 2
""")
P("P03", "numpy arrays", """
A numpy array holds many numbers, and arithmetic acts on all of them at once ("vectorised") — no loops needed.
""", code="""
import numpy as np                     # numpy: fast arrays of numbers
T_C = np.array([0.0, 15.0, 100.0])     # three temperatures in °C, stored in one array
print(T_C + 273.15)                    # add 273.15 to every element at once → [273.15 288.15 373.15]
""")
P("P04", "f-strings", """
`f"{x:.3g} Pa"` inserts the value of `x` into a string with 3 significant figures, followed by its unit;
`{x:.1f}` gives one decimal.
""", code="""
p = 101325.0                          # standard atmospheric pressure [Pa]
print(f"p = {p/1e3:.1f} kPa")         # p = 101.3 kPa  (divide by 1000, one decimal)
print(f"p = {p:.3g} Pa")              # p = 1.01e+05 Pa (three significant figures)
""")
nb.recap("R02", "SI prefixes", """
Prefixes scale a unit by powers of ten — n 10⁻⁹, µ 10⁻⁶, m 10⁻³, k 10³, M 10⁶, G 10⁹. `pint` handles them, so hPa in
§1.10, mN/m in §1.6 and nm in §1.4 need no hand conversion. A prefix changes the number, never the dimension.
""", where="pre-book; used wherever a number is printed")
nb.recap("R03", "Celsius and kelvin", r"""
$T_{[^\circ\mathrm C]} = T_{[\mathrm K]} - 273.15$. Laws such as $p = \rho R T$ need the **absolute** temperature:
15 °C is 288.15 K, and doubling the kelvin temperature (not the Celsius one) doubles $p$ at fixed $\rho$.
""", where="pre-book; used in C40")
nb.code("""
from fluidpy import ch01_introduction as ch01                         # the tested chapter-1 module: every function cites its book section
T_K = ch01.celsius_to_kelvin(np.array([-40.0, 0.0, 15.0, 20.0]))     # four Celsius readings → kelvin in one call
print(T_K)                                                             # [233.15 273.15 288.15 293.15]
print(ch01.kelvin_to_celsius(288.15))                                  # and back: 15.0 °C
print(Q_(15, "degC").to("K"))                                          # pint knows °C is an offset unit → 288.15 kelvin
""", explain="""
1. Imports the chapter module as `ch01`; all physics of this notebook is called from it.
2. Converts four Celsius readings to kelvin at once (a numpy array in, an array out).
3. Converts 288.15 K back to °C.
4. Lets `pint` do the same conversion; it treats °C as an offset unit, so it converts correctly but refuses to multiply
   a °C value — one more reason to compute in kelvin.
""")

# =====================================================================================================================
# A.3 §1.3 Solids, Liquids, and Gases
# =====================================================================================================================
nb.section("1.3", "Solids, Liquids, and Gases", intro="""
**What is this section about?** What makes something a fluid (it cannot resist a shear stress by staying put), how
fluids respond to being squeezed and pulled, and why gases and liquids differ (the spacing of their molecules).
""")
P("P05", "stress", """
**Stress = force per area** [Pa]. On a surface it splits into a **normal** part (pushing or pulling across the surface)
and a **shear** part (sliding along it). Shear turns a small cube into a leaning one; the lean angle γ [rad] is the
**shear strain**, and how fast it grows, dγ/dt [1/s], is the **strain rate**.
""", code="""
F_slide = 2.0            # sliding (tangential) force on a lid [N]
A_lid = 0.01             # lid area [m²] (10 cm × 10 cm)
print(F_slide / A_lid)   # shear stress = force / area → 200.0 Pa
""")
note("C02", "Fluid versus solid", r"""
A solid under a small, steady shear stress τ leans to a fixed angle and stops: $\gamma = \tau/G$, where G [Pa] is its
shear modulus (*gloss:* how stiff the solid is against shearing). A fluid **never stops**: any shear stress, however
small, makes γ grow for as long as it acts — $\gamma = \tau t/\mu$ for the simple fluids of this book. The rate law behind
the second formula is Newton's law of viscosity, developed in the C12 block of §1.5.
""", equation=r"\gamma_{\rm solid} = \tau/G \qquad \gamma_{\rm fluid} = \tau\,t/\mu")
P("P06", "np.linspace and np.logspace", """
`np.linspace(a, b, n)` gives n evenly spaced numbers from a to b; `np.logspace(p, q, n)` gives n numbers from 10ᵖ to
10ᵠ that are evenly spaced on a logarithmic axis (the same factor between neighbours).
""", code="""
print(np.linspace(0, 1, 5))      # [0.   0.25 0.5  0.75 1.  ] — equal steps
print(np.logspace(-9, -6, 4).tolist())   # [1e-09, 1e-08, 1e-07, 1e-06] — equal factors of 10
""")
nb.md("""
Our drawing of the idea (in the spirit of the book's Fig. 1.1, `N02`): the same shear stress applied to a solid and to a
fluid for 1.5 s, then removed.
""")
nb.figure("""
from fluidpy.core.style import COLORS                              # the house palette shared with the explainers
t = np.linspace(0, 3, 301)                                          # time [s]; the stress acts from 0 to 1.5 s
gam_solid = ch01.shear_deformation_history(t, 10.0, "solid", G=100.0, t_off=1.5)  # γ = τ/G while loaded (τ = 10 Pa)
gam_fluid = ch01.shear_deformation_history(t, 10.0, "fluid", mu=1.0, t_off=1.5)   # γ = τt/μ while loaded (μ = 1 Pa s)
fig, (a, b) = plt.subplots(1, 2, figsize=(9, 3.2))                   # two panels: solid left, fluid right
for ax in (a, b):                                                    # same shading on both panels
    ax.axvspan(0, 1.5, color=COLORS["grid"], alpha=0.7)              # grey band = shear stress applied
    ax.set_xlabel("time t [s]")                                      # x label with unit
a.plot(t, gam_solid, color=COLORS["muted"])                          # the solid: jumps to a fixed lean, springs back
a.set_ylabel("shear strain γ [rad]")                                 # y label with unit
a.set_title("Solid (G = 100 Pa): leans 0.1 rad, springs back")      # the message of the left panel
b.plot(t, gam_fluid, color=COLORS["teal"])                           # the fluid: γ keeps growing while loaded
b.set_title("Fluid (μ = 1 Pa s): keeps deforming, keeps it")        # the message of the right panel
plt.show()                                                           # display
""", see="""
Left: the solid jumps to a lean of 0.1 rad while the stress (grey band) acts and returns to zero when it is removed.
Right: the fluid's strain climbs steadily to 15 rad during the same 1.5 s and stays there afterwards.
""", read="""
A flat line while loaded means an elastic solid; a rising line means a fluid. The slope of the rising line is the
strain rate τ/μ = 10 s⁻¹ — exactly what Newton's law of viscosity (C12) will relate to the stress.
""", change="""
Halve μ and the fluid's line rises twice as steeply (to 30 rad); the solid's panel would not change at all, because
viscosity plays no role in an elastic solid.
""")
note("N01", "In-between materials", """
Some materials sit in between: **plastic** (yield-stress) materials behave as solids until the stress passes a yield
value, and **viscoelastic** ones (paint, egg white) partly spring back — non-Newtonian fluids, where Newton's law fails
(Ch. 4 §4.5; blood in Ch. 16 §16.3).
""")
note("C03", "Normal stress: compression, tension and cavitation", r"""
A normal stress can **push** (compression) or **pull** (tension). Solids and fluids both resist compression, but a liquid
pulled hard enough breaks: when its pressure drops below the **vapour pressure** (*gloss:* the pressure at which the
liquid boils at that temperature, ≈ 2.3 kPa for water at 20 °C) vapour cavities appear — **cavitation**, which eats ship
propellers. On a surface with unit normal $\hat{\mathbf n}$ the normal stress is the part of the traction vector
$\mathbf t$ (force per area on that surface) along $\hat{\mathbf n}$ (*gloss:* the dot product $\mathbf t\cdot\hat{\mathbf n}$
picks out the component of a vector along a unit vector). Pressure (C20) is the compressive normal stress of a fluid at
rest. Number (code below): a traction (3, 4, 0) Pa on a surface facing x has σ_n = 3 Pa of tension and |τ| = 4 Pa of shear.
""", equation=r"\sigma_n = \mathbf t\cdot\hat{\mathbf n}")
P("P07", "mole, kilomole, molecular weight and Avogadro's number", """
A kilomole (kmol) is 6.022×10²⁶ molecules — Avogadro's number per kilomole, A_o. The molecular weight M_w is the mass of
one kilomole in kilograms (air 28.96 kg/kmol, water 18.02). So one molecule weighs m = M_w/A_o, and a gas or liquid of
density ρ holds n = ρA_o/M_w molecules per m³. **This book uses kmol, not mol — a factor 1000 hides here.**
""", code="""
A_o = 6.02214076e26               # Avogadro's number per kilomole [1/kmol]
print(28.9644 / A_o)              # mass of one air molecule m = M_w / A_o → 4.81e-26 kg
print(1.225 * A_o / 28.9644)      # molecules per m³ in sea-level air n = ρ A_o / M_w → 2.55e+25
""")
note("C04", "Liquid versus gas", r"""
In a liquid the molecules touch their neighbours; in a gas they sit about ten molecule-sizes apart, so a gas is mostly
empty space and easy to squeeze. A liquid poured into a container forms a **free surface**; a gas fills whatever space it
is given. The typical spacing is $n^{-1/3}$ (one molecule per little cube): water ≈ 0.31 nm, sea-level air ≈ 3.4 nm — about
11 times more. How many molecules a sampling box holds is the question of C06.
""", equation=r"\ell_{\rm spacing} = n^{-1/3}, \qquad n = \rho A_o/M_w")
nb.code("""
sig_n, tau_vec, tau_mag = ch01.traction_components([3.0, 4.0, 0.0], [1.0, 0.0, 0.0])  # C03: traction (3,4,0) Pa on a surface facing +x
print(f"normal stress {sig_n} Pa (positive = tension), shear {tau_mag} Pa")          # 3.0 Pa tension, 4.0 Pa shear
s_water = ch01.mean_molecular_spacing(998.0, 18.015)        # C04: spacing n^(-1/3) in liquid water [m]
s_air = ch01.mean_molecular_spacing(1.225, ch01.M_W_AIR)    # C04: spacing in sea-level air [m]
print(f"water {s_water:.2e} m, air {s_air:.2e} m, ratio {s_air/s_water:.0f}")        # 3.11e-10 m, 3.40e-09 m, ratio 11
""", explain="""
1. `traction_components` splits a traction vector into its normal part (positive = pulling, tension) and its shear
   part — the C03 number.
2. `mean_molecular_spacing` computes n = ρA_o/M_w (primer P07) and returns n^(−1/3), for water and for air.
3. The ratio ≈ 11 is the "ten times farther apart" of C04.
""")

# =====================================================================================================================
# A.4 §1.4 Continuum Hypothesis — C06
# =====================================================================================================================
nb.section("1.4", "Continuum Hypothesis", intro="""
**What is this section about?** Why we may speak of density, velocity and temperature *at a point* although matter is made
of molecules, and the test — the Knudsen number — that says when this picture fails.
""")
core("C06", "The continuum hypothesis: a value at a point is the plateau of an average",
     "Matter is mostly empty space between molecules — so what can 'the density at this point' possibly mean?")
nb.md("""
#### The problem in plain words
A weather model stores one temperature per grid box; a pipe-flow formula uses the velocity *at* the wall. Zoom into either
and you find molecules flying about with nothing in between. Yet every equation from Chapter 2 onward treats ρ, **u**, p
and T as smooth functions of position. We need to know why that works — and when it stops working (microchannels, the
upper atmosphere, aerosol particles).
""")
nb.md("""
#### The idea
Measure density by counting: put a box of side L around the point, weigh the molecules inside, divide by the box volume.
Now grow the box and watch the reading.

```
box side L       1 nm       10 nm      100 nm     10 µm         1 mm          10 cm
molecules N      0.03       25         2.5e4      2.5e10        2.5e16        2.5e22
relative noise    —         20 %       0.6 %      6e-6          6e-9          —
reading          0 or huge  noisy      settling   ───── plateau ─────        drifts with the flow
                 └─ molecular noise ─┘            └─ continuum window ─┘     └─ flow's own variation ─┘
```

**The density at a point is the plateau value.** It exists because a wide range of box sizes lies between the molecular
scale and the scale on which the flow itself changes.
""")
P("P08", "temperature as molecular kinetic energy", r"""
Temperature measures the mean random kinetic energy of the molecules: $\tfrac12 m\langle|\mathbf u|^2\rangle =
\tfrac32 k_B T$, with Boltzmann's constant $k_B = 1.380649\times10^{-23}$ J/K. Hotter gas = faster molecules.
""", code="""
m_air = ch01.molecular_mass(ch01.M_W_AIR)              # mass of one air molecule [kg] (primer P07)
print(np.sqrt(3 * ch01.K_B * 288.15 / m_air))          # rms speed sqrt(3 k_B T / m) at 15 °C → ≈ 498 m/s
""")
P("P09", "Newton's second law and momentum", """
Momentum = mass × velocity [kg m/s]. Newton's second law: the net force on a body equals the rate of change of its
momentum (F = ma for a fixed mass). Draw every force acting on the body (a **free-body diagram**), add them with their
signs, and set the sum equal to m a. By the third law, whatever momentum the body gains, whatever pushed it loses.
""", code="""
m_ball, v_in, v_out = 0.1, 2.0, -2.0      # a 0.1 kg ball hits a wall at 2 m/s and bounces straight back
dP = m_ball * v_in - m_ball * v_out       # momentum handed to the wall = 0.2 − (−0.2) = 0.4 kg m/s
print(dP, dP / 0.01)                      # delivered in 0.01 s → average force 40 N
""")
P("P10", "np.random.default_rng", """
`rng = np.random.default_rng(0)` makes a *seeded* random-number generator: the same "random" numbers on every run, so
results are reproducible. `rng.normal(0, s, n)` draws n Gaussian numbers with standard deviation s; `rng.poisson(lam, n)`
draws n random counts with mean lam; `rng.random((n, 3))` draws uniform numbers in [0, 1).
""", code="""
rng = np.random.default_rng(0)                  # seeded generator → reproducible numbers
print(rng.normal(0, 2, 100_000).std())          # ≈ 2.0: the spread we asked for
print(rng.poisson(100, 100_000).std())          # ≈ 10: counts with mean 100 scatter by about √100
""")
P("P11", "Gaussian velocity components and mean square speed", r"""
In a gas at rest each velocity component ($u_x$, $u_y$, $u_z$) is Gaussian with zero mean and variance $k_BT/m$. The mean
square speed adds the three: $\langle|\mathbf u|^2\rangle = \langle u_x^2\rangle + \langle u_y^2\rangle +
\langle u_z^2\rangle = 3\langle u_x^2\rangle$, because no direction is special.
""", code="""
u = ch01.maxwellian_velocities(100_000, 288.15, m_air, seed=0)   # 100 000 molecules' velocities [m/s], shape (100000, 3)
print((u**2).sum(axis=1).mean() / (u[:, 0]**2).mean())           # ⟨|u|²⟩ / ⟨u_x²⟩ → ≈ 3.0
""")
note("C05", "Pressure is momentum delivered by molecular impacts", r"""
Molecules hitting a wall bounce back, and every bounce pushes the wall. On each square metre of a wall in sea-level air
there are about 10²⁷ impacts per second, so their average push per area is perfectly steady — that average is the
pressure. The derivation below turns this picture into a formula.
""", equation=r"p = \tfrac13\,n\,m\,\langle|\mathbf u|^2\rangle")
D("D34", "Pressure from molecular impacts (kinetic theory behind §1.4; our addition)", ref="",
  goal="""Turn "pressure is the average push of molecules hitting a wall" into a formula, and see that it gives the gas
law pV = nk_BT used in §1.9.""",
  assumptions="""Ideal gas: molecules do not interact except by brief collisions (step 2: straight flight to the wall) ·
an elastic, smooth wall (step 1) · isotropic velocities, gas at rest (steps 5–7) · equilibrium, so the kinetic meaning
of T applies (step 8).""",
  start=(r"p = \dfrac{\text{momentum delivered to the wall}}{\text{area} \times \text{time}}",
         "force is momentum per time (Newton), and pressure is force per area."),
  plan=["Find the momentum one bounce gives the wall.", "Count how many molecules of one velocity hit in a time dt.",
        "Add over all velocities.", "Use isotropy and the kinetic meaning of temperature."],
  uses=["Newton's second law and momentum (P09)", "Gaussian velocity components and mean square speed (P11)",
        "temperature as molecular kinetic energy (P08)"],
  steps=[
      ("Momentum given by one bounce", r"\Delta P_{\rm wall} = 2\,m\,u_x",
       "An elastic bounce off a smooth wall normal to x reverses u_x and keeps u_y, u_z; the molecule's momentum changes "
       "by −2mu_x, so by Newton's third law the wall gains +2mu_x. This is the push we want to count.",
       "Each hit hands the wall twice the molecule's normal momentum."),
      ("Count hits of one velocity class in time dt", r"dN = n_{u_x}\,A\,u_x\,dt \quad (u_x > 0)",
       "Only molecules moving toward the wall (u_x > 0) and within a distance u_x dt of it arrive in dt: they fill a slab "
       "of volume A u_x dt, and n_{u_x} is the number per volume with that velocity.",
       "Faster molecules reach the wall from farther away."),
      ("Multiply hits by momentum per hit", r"dP = 2\,m\,n_{u_x}\,u_x^2\,A\,dt",
       "Total momentum = number of hits × momentum per hit (steps 1 and 2); u_x appears twice, once from how many arrive "
       "and once from how hard each hits.",
       "The delivered momentum grows with the square of the speed."),
      ("Divide by area and time", r"p_{u_x} = 2\,m\,n_{u_x}\,u_x^2",
       "Force = momentum per time (P09), pressure = force per area; dividing by A dt gives this velocity class's share "
       "of the pressure.",
       "One velocity class adds this much pressure."),
      ("Add all classes moving toward the wall",
       r"p = 2m\sum_{u_x>0} n_{u_x}u_x^2 = 2m\cdot\tfrac{n}{2}\langle u_x^2\rangle",
       "By symmetry half of the n molecules per volume move toward the wall, and those have the same mean square u_x as "
       "the whole gas. We sum because pressure is the total push.",
       "The half of the molecules heading for the wall carry the average square speed."),
      ("Cancel the 2 and the ½", r"p = n\,m\,\langle u_x^2\rangle",
       "2 × ½ = 1; nothing else changes, and the formula now uses only averages over the whole gas.",
       "Pressure is number density × mass × mean square normal velocity."),
      ("Replace ⟨u_x²⟩ using isotropy", r"p = \tfrac13\,n\,m\,\langle|\mathbf u|^2\rangle",
       "In a gas at rest no direction is special, so ⟨u_x²⟩ = ⟨u_y²⟩ = ⟨u_z²⟩ and their sum is ⟨|u|²⟩ (P11); hence "
       "⟨u_x²⟩ = ⅓⟨|u|²⟩. We want the speed, not one component.",
       "Pressure is one third of density times mean square speed."),
      ("Use the kinetic meaning of temperature", r"p = n\,k_B\,T",
       "½m⟨|u|²⟩ = (3/2)k_BT (P08), so m⟨|u|²⟩ = 3k_BT and the ⅓ cancels the 3. This links the molecular picture to a "
       "measurable temperature.",
       "Pressure = molecules per volume × k_B × temperature."),
  ],
  result=(r"p = \tfrac13\,n\,m\,\langle|\mathbf u|^2\rangle = n\,k_B T",
          "multiplying by a volume V with n = N/V gives pV = Nk_BT, the molecular gas law (1.21) stated in C38."),
  interpret="""A steady pressure is the average of an enormous number of tiny, random pushes — which is why it is well
defined only for boxes holding many molecules (D35 below). Heating a gas at fixed n raises p because molecules hit harder
and more often. The formula fails for dense gases and liquids, where molecules interact during their whole flight.""",
  check="""Units: m⁻³ · kg · m² s⁻² = kg m⁻¹ s⁻² = Pa ✓. Number: sea-level air, n = 2.55×10²⁵ m⁻³, m = 4.81×10⁻²⁶ kg,
⟨|u|²⟩ = 3k_BT/m = 2.48×10⁵ m²/s² → p = ⅓ × 2.55×10²⁵ × 4.81×10⁻²⁶ × 2.48×10⁵ = 1.013×10⁵ Pa ✓. The next code cell
samples molecules and gets within 2 % of nk_BT both ways.""",
  traps="losing the factor 2 of the bounce; forgetting that only half the molecules move toward the wall; using "
        "⟨u_x²⟩ = ⟨|u|²⟩ instead of ⅓⟨|u|²⟩.")
nb.code("""
n_air = ch01.number_density(1.225, ch01.M_W_AIR)           # molecules per m³ in sea-level air, n = ρ A_o / M_w [1/m³]
N_mol = 200_000 if not FAST else 50_000                    # how many molecules we sample (fewer in FAST mode)
u = ch01.maxwellian_velocities(N_mol, 288.15, m_air, seed=0)  # their random velocities at 288.15 K [m/s], shape (N, 3)
p_kin = ch01.molecular_pressure(n_air, m_air, u)           # D34 step 7: p = ⅓ n m ⟨|u|²⟩ from the sample [Pa]
p_wall = ch01.wall_impact_pressure(u, m_air, n_air)        # D34 steps 1–5: count wall hits in a short time, add 2 m u_x [Pa]
p_nkT = n_air * ch01.K_B * 288.15                          # D34 step 8: p = n k_B T [Pa]
print(f"n k_B T = {p_nkT:,.0f} Pa")                        # ≈ 101 327 Pa: sea-level pressure from molecules
print(f"⅓ n m <|u|²> = {p_kin:,.0f} Pa  (ratio {p_kin/p_nkT:.3f})")      # within sampling scatter of 1
print(f"wall impacts = {p_wall:,.0f} Pa  (ratio {p_wall/p_nkT:.3f})")    # within sampling scatter of 1
""", explain="""
1. `number_density` and `molecular_mass` (primer P07) give n and m for air.
2. `maxwellian_velocities` draws Gaussian velocity components with variance k_BT/m (primers P10, P11).
3. `molecular_pressure` evaluates step 7 of D34 on the sample.
4. `wall_impact_pressure` literally counts the molecules that would reach a wall in a short time and adds their 2mu_x —
   steps 1–5.
5. Step 8, p = nk_BT, is the exact value; both sampled estimates agree with it to within their random scatter.
""")
P("P12", "Poisson counting statistics", r"""
When very many independent things each have a small chance to land in a box, the count N is **Poisson**: its variance
equals its mean $\bar N$, so its standard deviation is $\sqrt{\bar N}$ and its **relative** scatter is $1/\sqrt{\bar N}$.
""", code="""
rng = np.random.default_rng(0)            # seeded generator (primer P10)
c = rng.poisson(25, 100_000)              # 100 000 boxes that expect 25 molecules each
print(c.std() / c.mean(), 1 / np.sqrt(25))   # measured relative scatter ≈ 0.20 = 1/√25
""")
P("P13", "power laws and log–log plots", r"""
$y = a\,x^k$ is a straight line of slope k on log–log axes (`ax.loglog`), because $\log y = \log a + k\log x$. So the
slope you read off a log–log plot is the exponent; `np.polyfit(log x, log y, 1)[0]` measures it.
""", code="""
x = np.logspace(0, 3, 4)                                   # 1, 10, 100, 1000
y = 5 * x**-1.5                                            # a power law with exponent −1.5
print(np.polyfit(np.log10(x), np.log10(y), 1)[0])          # slope of the straight line on log–log axes → −1.5
""")
P("P43", "exponent rules", r"""
$a^m a^n = a^{m+n}$; $(a^m)^n = a^{mn}$; $a^{-m} = 1/a^m$; $(a/b)^m = a^m/b^m$; and $(ab)^m = a^m b^m$. For example
$(L^3)^{-1/2} = L^{-3/2}$.
""", code="""
print(2.0**0.4 * 2.0**1.0, 2.0**1.4)            # a^m a^n = a^(m+n): both 2.639
print((2.0**0.5)**2, 1 / 2.0**-1.4, 2.0**1.4)   # (a^m)^n = a^(mn) and a^(−m) = 1/a^m
print((3.0 * 4.0)**-0.5, 3.0**-0.5 * 4.0**-0.5) # (ab)^k = a^k b^k: both 0.2887
""")
nb.md(r"""
#### The maths, step by step
1. Around a point **x** take a small volume δV [m³] holding a mass δm [kg].
2. Define the continuum density as the ratio $\rho(\mathbf x) = \delta m/\delta V$ (§1.4).
3. The ratio is useful only if δV holds so many molecules that the scatter of the count is negligible (**lower bound**,
   worked out in D35 below) **and** δV is so small that the flow's own density does not change across it (**upper
   bound**, set by the flow's length scale L_flow [m]).
4. Every field variable of the book — **u**, p, T — is defined the same way, as a plateau average.

Symbols: n [m⁻³] molecules per volume; m [kg] mass of one molecule; $\bar N = n\,\delta V$ the expected count;
ε [–] the fractional amplitude of the flow's own density variation.
""")
D("D35", "How noisy is a density measured in a small box? (§1.4; our addition)", ref="",
  goal="Find how the scatter of a box's density reading depends on the box size — the lower end of the continuum window.",
  assumptions="""Molecules placed independently (ideal gas; step 1) · the box holds a tiny fraction of all molecules
(step 1: Poisson limit) · the macroscopic density does not vary across the box (step 3; the flow's variation is the
separate upper bound).""",
  start=(r"\bar N = n\,\delta V = n\,L^3", "a box of side L in a gas with n molecules per m³ holds N̄ molecules on average."),
  plan=["The molecule count is Poisson.", "Scale its scatter into a density scatter.", "Divide by the mean.",
        "Write the result in terms of L."],
  uses=["Poisson counting statistics (P12)", "power laws and log–log plots (P13)",
        "exponent rules (P43)"],
  steps=[
      ("Use Poisson statistics for the count", r"\operatorname{Var}(N) = \bar N",
       "Each of very many independent molecules has a tiny chance to be inside the box; such counts are Poisson, whose "
       "variance equals the mean (P12). This is where the randomness enters.",
       "The count scatters by about √N̄."),
      ("Take the square root", r"\sigma_N = \sqrt{\bar N}",
       "The standard deviation is by definition the square root of the variance, and it has the units of the count.",
       "A box expecting 100 molecules typically sees 90 to 110."),
      ("Turn counts into densities", r"\sigma_\rho = \dfrac{m\,\sigma_N}{L^3},\qquad \rho = \dfrac{m\bar N}{L^3}",
       "The reading is δm/δV = mN/L³ (C06); multiplying a random number by the constant m/L³ multiplies both its mean and "
       "its standard deviation by that constant.",
       "The density reading inherits the count's scatter."),
      ("Divide by the mean", r"\dfrac{\sigma_\rho}{\rho} = \dfrac{\sigma_N}{\bar N} = \bar N^{-1/2}",
       "m/L³ cancels; √N̄/N̄ = N̄^(−1/2). A *relative* scatter is what decides whether a value is usable.",
       "The relative noise is one over the square root of the expected count."),
      ("Write it in terms of the box side", r"\dfrac{\sigma_\rho}{\rho} = (nL^3)^{-1/2} = n^{-1/2}L^{-3/2}",
       "Substitute N̄ = nL³ and use (ab)^k = a^k b^k and (L³)^(−1/2) = L^(−3/2). The exponent is the slope on a log–log "
       "plot (P13).",
       "Every factor 10 in box side cuts the noise by 10^1.5 ≈ 32."),
  ],
  result=(r"\dfrac{\sigma_\rho}{\rho} = (n\,L^3)^{-1/2}",
          "the relative noise of a box density falls as the box side to the power −3/2."),
  interpret="""The continuum's lower limit is statistical: a density exists once a box holds, say, 10⁶ molecules
(noise 0.1 %), i.e. L ≈ 0.3 µm in sea-level air. In a liquid the molecules are not independent (they touch), so the true
scatter is smaller than Poisson and this estimate is an upper bound. In the explainer below, the grey band's half-width
at the marker is exactly this number.""",
  check="""Dimensionless ✓ (nL³ is a count). Limits: L → 0 gives infinite noise (a box with no molecules), L → ∞ gives
zero ✓. Number: a 10 µm cube of air, N̄ = 2.55×10¹⁰ → 6.3×10⁻⁶ ✓; the sweep below follows a line of slope −1.5, and a
from-scratch molecule count agrees within its statistical scatter.""",
  traps="writing L^(−1/2) (it is the volume, not the side, under the root); confusing the absolute scatter σ_ρ with "
        "the relative scatter σ_ρ/ρ.")
nb.worked_example("a 10 µm cube of sea-level air", """
1. n = ρA_o/M_w = 1.225 × 6.022×10²⁶ / 28.96 = 2.55×10²⁵ m⁻³.
2. δV = (10⁻⁵ m)³ = 10⁻¹⁵ m³.
3. N̄ = nδV = 2.55×10¹⁰ molecules.
4. δm = N̄m = 2.55×10¹⁰ × 4.81×10⁻²⁶ kg = 1.225×10⁻¹⁵ kg, so δm/δV = 1.225 kg/m³.
5. Relative noise 1/√N̄ = 1/(1.6×10⁵) = 6.3×10⁻⁶ — utterly negligible.
6. The same air in a 10 nm cube: N̄ = 25, noise = 1/5 = 20 % — no useful density at all.
""")
P("P14", "tuple unpacking", """
A function may return several values packed in a tuple; `a, b = f()` unpacks them into two names in one line.
""", code="""
result = (1.225, 0.01)       # a function result with two values: (mean, standard deviation)
mean, std = result           # unpack them into two names
print(mean, std)             # 1.225 0.01
""")
nb.code("""
L = np.logspace(-9, 0, 37 if not FAST else 19)                  # box sides from 1 nm to 1 m, evenly spaced in log [m]
eps, L_flow = 0.2, 1.0                                          # the flow's density varies ±20 % over 1 m
mean, std = ch01.sample_density(L, n_air, m_air, n_samples=200 if not FAST else 60,   # sampled box densities (continued below)
                                variation=eps, L_flow=L_flow, seed=0)  # 200 random boxes per size: mean and std [kg/m³]
rel = std / mean                                                # measured relative scatter [–]
noise = ch01.density_noise_expected(L, n_air)                   # D35 prediction (n L³)^(-1/2) [–]
rho_pt = 1.225 * (1 + eps)                                      # the density AT the point (the box sits on a crest) [kg/m³]
box_avg = ch01.box_average_density(L, 1.225, eps, L_flow)       # the noiseless box average: blurs the variation [kg/m³]
drift = np.abs(box_avg / rho_pt - 1)                            # how far the box average has drifted from the point value
for Lq in (1e-8, 1e-5, 1e-1):                                   # three sizes: tiny, window, too big
    i = np.argmin(np.abs(np.log10(L / Lq)))                     # index of the nearest box size in the sweep
    print(f"L = {L[i]:.0e} m: measured noise {rel[i]:.2e}, D35 {noise[i]:.2e}, drift {drift[i]:.1e}")   # one table row
""", explain="""
1. Builds 37 box sizes from 1 nm to 1 m (primer P06).
2. `sample_density` places 200 random boxes per size in air whose density varies by ±20 % over 1 m, draws each box's
   molecule count as in D35 around the box average, and returns the mean and standard deviation of the readings
   (tuple unpacking, primer P14).
3. `rel` is the measured relative scatter and `noise` the D35 prediction — they agree.
4. `box_average_density` is the average of the flow's smooth variation over the box (a sinusoid averaged over a box
   gives the factor sinc(πL/L_flow) = sin(πL/L_flow)/(πL/L_flow)); `drift` measures how far it has moved from the value at the point.
5. The table: at 10 nm the noise dominates, at 10 µm both are negligible (the plateau), at 10 cm the drift has grown.
""")
P("P15", "assert np.allclose", """
`assert np.allclose(a, b, rtol=…)` stops the notebook with an error if two numbers or arrays differ by more than the
relative tolerance `rtol`. We use it as proof that a hand-written version matches the tested library.
""", code="""
assert np.allclose([1.0, 2.0], [1.0, 2.0 + 1e-9])   # passes silently: equal to within the default tolerance
print("the two arrays agree")                        # reached only if the assert passed
""")
nb.md("""
#### From scratch: place molecules and count them
The library draws counts from the Poisson law. Here we do it the naive way — scatter molecules uniformly in a unit cube,
drop random boxes into it and count what falls inside — and check that the relative scatter follows D35.
""")
nb.check_agree("""
rng = np.random.default_rng(1)                               # seeded generator (primer P10)
n_mol = 100_000 if not FAST else 40_000                      # molecules in the unit cube → number density n = n_mol per m³
pos = rng.random((n_mol, 3))                                 # their positions, uniform in [0, 1)³ [m]
n_boxes = 300 if not FAST else 120                           # how many random boxes per size
for b in (0.05, 0.1, 0.2):                                   # three box sides [m]
    corners = rng.random((n_boxes, 3)) * (1 - b)             # random lower corners that keep the box inside the cube
    counts = np.array([np.all((pos >= c) & (pos < c + b), axis=1).sum() for c in corners])  # molecules inside each box
    rel_mine = counts.std(ddof=1) / counts.mean()            # measured relative scatter of the count (= of the density)
    rel_d35 = ch01.density_noise_expected(b, n_mol)          # D35: (n b³)^(-1/2)
    m_lib, s_lib = ch01.sample_density(b, n_mol, 1.0, n_samples=n_boxes)   # the library's sampler, molecule mass 1 kg
    print(f"b = {b}: N̄ = {n_mol*b**3:.0f}, mine {rel_mine:.3f}, D35 {rel_d35:.3f}, library {s_lib/m_lib:.3f}")   # compare the three estimates
    assert np.allclose(rel_mine, rel_d35, rtol=0.25)         # both random: agreement within the statistical scatter
    assert np.allclose(rel_mine, s_lib / m_lib, rtol=0.25)   # the library's sampler behaves like real counting
""")
nb.md("""
Statistical agreement (to 25 %) is the right test here: both sides are random, and the scatter of a standard deviation
estimated from a few hundred boxes is several per cent. The deterministic law is D35.
""")
nb.figure("""
from fluidpy.core.style import savefig                                  # saves a copy to outputs/ch01/ (not committed)
fine = np.logspace(-9, 0, 400)                                          # a fine grid of box sides for the smooth curves [m]
noise_f = ch01.density_noise_expected(fine, n_air)                      # D35 line (slope −3/2 on log–log axes)
drift_f = np.abs(ch01.box_average_density(fine, 1.225, eps, L_flow) / rho_pt - 1)  # drift of the box average
ok = (noise_f < 1e-3) & (drift_f < 1e-3)                                # the continuum window: both below 0.1 %
fig, ax = plt.subplots(figsize=(7.5, 4.2))                              # one log–log panel
ax.loglog(L, rel, "o", color=COLORS["teal"], ms=5, label="measured relative scatter (200 boxes)")   # the samples
ax.loglog(fine, noise_f, "--", color=COLORS["muted"], label=r"D35: $(nL^3)^{-1/2}$")                # prediction
ax.loglog(fine[drift_f > 1e-16], drift_f[drift_f > 1e-16], color=COLORS["orange"], label="drift of the box average")   # drift (zeros hidden on a log axis)
ax.axvspan(fine[ok][0], fine[ok][-1], color=COLORS["teal"], alpha=0.12, label="continuum window (both < 0.1 %)")   # shade the window
ax.axvline(s_air, color=COLORS["accent"], lw=1, ls=":")                 # molecular spacing of air (3.4 nm)
ax.text(s_air * 1.3, 1e-11, "molecular\\nspacing", color=COLORS["accent"], fontsize=9)   # label the spacing line
ax.axvline(L_flow, color=COLORS["orange"], lw=1, ls=":")                # the flow's length scale (1 m)
ax.set_ylim(1e-13, 30)                                                  # fixed range so all three regimes show
ax.set_xlabel("box side L [m]")                                         # x label with unit
ax.set_ylabel("relative error of the density [–]")                      # y label (dimensionless)
ax.set_title(f"A density 'at a point' exists between {fine[ok][0]*1e6:.1f} µm and {fine[ok][-1]*100:.0f} cm")   # the message, with the window's edges
ax.legend(loc="upper center", fontsize=8.5)                             # legend in the empty upper middle
savefig(fig, "ch01", "c06_continuum_window")                            # keep a copy for review
plt.show()                                                              # display
""", see="""
Teal dots (measured scatter) lying on a grey dashed line of slope −3/2, an orange curve rising at the right, and a shaded
window about five decades wide where both are below 10⁻³.
""", read="""
Inside the shaded window any box gives the same density to better than 0.1 % — that shared number is ρ(**x**). Left of
it the reading is dominated by molecular noise; right of it the box is so big that it averages over the flow's own
variation, and the "value at the point" gets blurred. The window's width, not molecules being absent, is what makes the
continuum description work.
""", change="""
In water n is ≈ 1300× larger, so the grey line drops by √1300 ≈ 36 and the window starts about 11× smaller boxes. In a
thin boundary layer (the thin layer next to a wall where the flow speed changes quickly; L_flow = 1 mm) the orange curve moves three decades to the left and the window shrinks from the
right.
""")
P("P16", "animate and show_animation", """
`animate(update, frames=n, fig=fig)` calls your function `update(i)` once for every frame number i; `update` moves
artists that are already drawn (with `set_data`, `set_offsets`, `set_text` …). `show_animation(anim, player="frames")`
embeds it with ◀ ▮▮ ▶ step buttons; `player="video"` makes a smooth MP4.
""", code="""
from fluidpy.core.anim import animate               # builds a matplotlib animation from an update function
fig, ax = plt.subplots(figsize=(3.2, 1.4))          # a tiny figure
(dot,) = ax.plot([0], [0], "o", color=COLORS["accent"])   # draw the first frame once: a dot at x = 0
ax.set_xlim(-0.5, 9.5); ax.set_ylim(-1, 1)          # fixed axes, so only the dot moves
def update(i):                                      # called for frames i = 0 … 9
    dot.set_data([i], [0])                          # move the dot to x = i
    return (dot,)                                   # the changed artist
show_animation(animate(update, frames=10, fig=fig, interval=200), player="frames")   # step through the 10 frames
""")
nb.animation("""
from fluidpy.core.anim import animate                                     # animation helper (primer P16)
frames = 37 if not FAST else 19                                           # one frame per quarter decade of box size
L_anim = np.logspace(-9, 0, frames)                                       # box side for each frame [m]
import warnings                                                           # standard library: silence expected warnings
with warnings.catch_warnings():                                           # a 1-box "standard deviation" is undefined …
    warnings.simplefilter("ignore", RuntimeWarning)                       # … so numpy warns; we only use the mean
    reading = np.array([ch01.sample_density(Li, n_air, m_air, n_samples=1, variation=eps, L_flow=L_flow, seed=i)[0]   # (continued)
                        for i, Li in enumerate(L_anim)])                  # ONE random box reading per size [kg/m³]
fig, (axl, axr) = plt.subplots(1, 2, figsize=(8.4, 3.5), gridspec_kw=dict(width_ratios=[1, 1.35]))   # left: the gas, right: the readings
axl.set_xticks([]); axl.set_yticks([]); axl.set_aspect("equal")          # left: a slice of the gas, 4 box sides wide
field = axl.imshow(np.zeros((2, 60)), extent=(-2, 2, -2, 2), vmin=1.0, vmax=1.5, cmap="Blues", origin="lower")   # the smooth density field (blank at first)
dots = axl.scatter([], [], s=6, color=COLORS["ink"])                      # molecules (only while few enough to draw)
axl.add_patch(plt.Rectangle((-0.5, -0.5), 1, 1, fill=False, ec=COLORS["accent"], lw=2))   # the sampling box
label = axl.set_title("")                                                 # shows L and the expected count
axr.axhline(rho_pt, color=COLORS["muted"], ls="--", lw=1)                 # the density at the point (1.47 kg/m³)
axr.text(2e-9, rho_pt + 0.05, "value at the point", color=COLORS["muted"], fontsize=8)   # label the dashed line
(trace,) = axr.semilogx([], [], "-", color=COLORS["teal"], lw=1)          # readings so far
(now,) = axr.semilogx([], [], "o", color=COLORS["orange"], ms=9)          # the current reading, enlarged
axr.set_xlim(5e-10, 2); axr.set_ylim(0, 3)                                # fixed axes (readings above 3 are clipped)
axr.set_xlabel("box side L [m]"); axr.set_ylabel("one box reading δm/δV [kg/m³]")   # axis labels with units
rng_anim = np.random.default_rng(3)                                       # seeded positions for the dots
def update(i):                                                            # frame i: box side L_anim[i]
    Li = L_anim[i]                                                        # current box side [m]
    n_view = n_air * (4 * Li)**2 * Li                                     # expected molecules in the drawn slice (4L × 4L × L)
    if n_view < 2000:                                                     # few enough: draw them as dots
        k = rng_anim.poisson(n_view)                                      # how many fall in the slice this time
        dots.set_offsets(rng_anim.uniform(-2, 2, (k, 2)))                 # random positions in box units
        field.set_data(np.full((2, 60), np.nan))                          # no smooth field yet
    else:                                                                 # too many: show the smooth density field instead
        dots.set_offsets(np.empty((0, 2)))                                # hide the dots
        x = L_flow / 4 + np.linspace(-2, 2, 60) * Li                      # positions across the view [m] (centred on a crest)
        field.set_data(np.tile(1.225 * (1 + eps * np.sin(2 * np.pi * x / L_flow)), (2, 1)))   # ρ(x) of the flow
    label.set_text(f"L = {Li:.1e} m, expected N = {n_air*Li**3:.2g}")    # title with the numbers
    trace.set_data(L_anim[:i + 1], np.minimum(reading[:i + 1], 2.95))     # readings so far (clipped at the top)
    now.set_data([Li], [min(reading[i], 2.95)])                           # the current one
    return dots, field, trace, now, label   # the artists that changed
show_animation(animate(update, frames=frames, fig=fig, interval=350), player="frames")   # step decade by decade
""", explain="""
1. Frame i sets the box side L = 10^(−9 + i/4) m; `sample_density(..., n_samples=1)` gives one random reading per size.
2. Left panel: while the drawn slice would hold fewer than 2000 molecules they are drawn as dots; after that the slice
   shows the flow's smooth density field (it only varies visibly once L approaches the 1 m flow scale).
3. Right panel: the reading so far, with the current one enlarged; readings larger than 3 kg/m³ are clipped at the top.
""")
nb.md("""
**What you see.** Dots and a wildly jumping reading (zero, or far off the top) for nanometre boxes; then the dots become
too many to draw and the reading settles onto the dashed line; in the last frames the reading slides down as the box
spans the flow's own variation.

**How to read it.** The moment the molecules become too many to draw is roughly where the reading stops jumping —
beyond it you are in the continuum window.

**What would change if…** the flow had no density variation (ε = 0): the right end would stay flat on the line forever,
and the only limit would be the molecular one.
""")
P("P17", "slider_figure", """
`slider_figure(fn, name, values)` calls `fn(value)` for every slider position **now**, in Python, and stores all the
curves in one plotly figure. Dragging the slider just switches which curves are visible — so it keeps working on the
published web page without any Python running. `fn` returns `{trace name: (x, y)}`.
""", code="""
from fluidpy.core.interact import slider_figure                       # precomputed-slider plotly figures
xs = np.linspace(0, 1, 11)                                            # x values
fig = slider_figure(lambda a: {"y = a·x": (xs, a * xs)}, "a", [1, 2, 3],   # fn, name, values …
                    xlabel="x [–]", ylabel="y [–]", title="Drag a: the slope changes", height=320)   # 3 slider steps
fig.show()                                                            # interactive figure below the cell
""")
nb.plotly("""
z_alt = np.linspace(0, 80e3, 17)                                          # 17 geopotential altitudes from 0 to 80 km [m]
L_body = np.logspace(-9, 1, 120)                                          # body or channel sizes from 1 nm to 10 m [m]
def kn_curves(z_km):                                                      # called once per slider position (altitude in km)
    T, p, rho = ch01.standard_atmosphere(z_km * 1e3)                      # US Standard Atmosphere 1976 at that geopotential height
    l = ch01.mean_free_path_air(T, p)                                     # mean free path of air there [m] (Jennings)
    return {"Kn = l/L": (L_body, ch01.knudsen_number(l, L_body)),        # Knudsen number for every body size
            "Kn = 0.01 (continuum limit)": (L_body[[0, -1]], [0.01, 0.01]),   # guide line
            "Kn = 0.1 (slip → transition)": (L_body[[0, -1]], [0.1, 0.1])}    # guide line
fig = slider_figure(kn_curves, "altitude", z_alt / 1e3, unit="km", xlabel="body size L [m]", ylabel="Kn [–]",   # one slider position per altitude
                    title="The continuum fails for smaller and smaller bodies as the air thins", active=0)   # start at sea level
fig.update_xaxes(type="log", range=[-9, 1])                               # log axes: ranges are given in powers of ten
fig.update_yaxes(type="log", range=[-11, 4])                              # Kn from 1e-11 to 1e4
fig.show()                                                                # drag the altitude slider
""", explain="""
1. For each of 17 altitudes, `standard_atmosphere` returns the US Standard Atmosphere 1976 temperature, pressure and
   density (a cited public table, see reference/ch01/SOURCES.md). Its height argument is the *geopotential* altitude,
   which corrects for g weakening with height; at 80 km it is about 1 km less than the geometric height, which does not
   matter for the orders of magnitude read here.
2. `mean_free_path_air` computes the mean free path l with Jennings' formula and Sutherland's viscosity.
3. `knudsen_number` gives Kn = l/L for every body size; the two flat lines mark Kn = 0.01 and 0.1.
""")
see_read_change(
    'A falling straight line of Kn against body size on log–log axes, with two flat guide lines at Kn = 0.01 and 0.1; the slider moves the line up as the altitude grows.',
    'Where the line lies below the lower guide, a body of that size sees a continuum; between the guides, slip flow; above them the gas behaves as separate molecules. Read the body size where the line crosses Kn = 0.01 — that is the smallest object the continuum description suits at that altitude: about 7 µm at sea level, while at 80 km (l ≈ 4.4 mm) even a 40 cm instrument is outside it.',
    'At a fixed altitude, doubling the temperature at the same pressure would halve the density and roughly double l, moving the whole line up by a factor of 2 or more (warmer air is also more viscous) — a small shift on these axes.',
)
note("C07", "The Knudsen number", r"""
The same scale question arises for *motion*: molecules fly a **mean free path** l between collisions. If the body or
channel size L is not much larger than l, the gas cannot smooth out velocity differences and the continuum description
fails. The test is the Knudsen number: continuum for Kn ≲ 0.01, slip flow for 0.01–0.1, transition up to about 10,
free-molecular flow beyond. Air at 300 K around a 1 µm particle: Kn = 0.067 — already slip flow (code below).
""", equation=r"Kn = l/L", ref="§1.4")
note("N03", "The mean free path of air and water", """
For air at 300 K and 1 atm, l ≈ 67 nm (Jennings' formula with Sutherland's viscosity — our value, the published order of
magnitude). In liquid water the molecules are in contact, so l is comparable with the 0.3 nm molecular size and Kn is
negligible for any practical L.
""")
note("C24", "The fluid particle (book §1.8, taught here)", """
A **fluid particle** is the plateau box followed as the flow carries it along. It (1) always contains the same molecules;
(2) is large enough that its properties are well defined — the lower end of the window; and (3) is small enough, with
collisions frequent enough, that it **relaxes** to equilibrium (settles after any disturbance; the relaxation time is
defined with C23 in §1.8) much faster than the flow changes it. A 10 µm particle of air holds 2.5×10¹⁰ molecules (noise
6×10⁻⁶) and each molecule collides about every 1.4×10⁻¹⁰ s. Chapter 3 turns this into the material volume.
""")
nb.code("""
l_300 = ch01.mean_free_path_air(300.0)                          # N03: mean free path of air at 300 K, 1 atm [m]
print(f"l = {l_300:.2e} m")                                     # ≈ 6.7e-08 m = 67 nm
print(f"Kn for a 1 µm particle = {ch01.knudsen_number(l_300, 1e-6):.3f}")          # C07: 0.067 → slip flow
print(f"collision time = {ch01.collision_time(l_300, 300.0, m_air):.2e} s")        # C24: l / mean speed ≈ 1.4e-10 s
T80, p80, _ = ch01.standard_atmosphere(80e3)                    # temperature and pressure at 80 km geopotential height (USSA-1976)
print(f"l at 80 km = {ch01.mean_free_path_air(T80, p80):.1e} m")                   # ≈ 4.4e-03 m: millimetres
rho80 = ch01.standard_atmosphere(80e3)[2]                                          # density at 80 km [kg/m³]
print(f"sea-level density / density at 80 km = {ch01.standard_atmosphere(0.0)[2] / rho80:,.0f}")   # how much thinner the air is
""", explain="""
1. `mean_free_path_air` at room conditions (N03).
2. `knudsen_number` for a 1 µm aerosol particle (C07).
3. `collision_time` = mean free path ÷ mean molecular speed: the relaxation clock of a fluid particle (C24).
4. The same mean free path 80 km up; the last line prints how many times less dense the air is there (tens of
   thousands), which is why l has grown from nanometres to millimetres.
""")
P("P18", "show_viz", """
`show_viz("ch01", slug)` embeds one of this chapter's interactive explainers: in Jupyter it loads the file from the
repository, in Colab from the project's web site, and on the published page the explainer fills the whole window. Each
explainer has a Walkthrough, Explore, Explain (every number worked out with your settings), Derivation, Equations, Code and
Check-yourself tabs.
""", code="""
print(show_viz.__doc__.splitlines()[0])   # the one-line description of the embedding helper imported in the setup cell
""")
nb.explainer("continuum_averaging_volume", heading="When does 'density at a point' make sense?", why="""
The continuum window is a sweep across nine decades of box size with randomness in it. Dragging the box, resampling and
switching the medium shows noise, plateau and drift in a way a single static curve cannot — and the Derivation tab walks
through D35 with your box.
""", tries=[
    "Press ▶ in the Walkthrough and watch the teal samples squeeze into the grey band as the box grows.",
    "Switch the medium to liquid water and find the box size where the reading first stays within 1 %.",
    "Set L_flow = 1 mm (a thin boundary layer): does a continuum window still exist?",
    "Choose '80 km air' and read Kn for a 1 cm sensor.",
])
nb.code("""
z_grid = np.linspace(0, 50e3, 5001)                                         # geopotential heights every 10 m [m]
rho_ratio_z = ch01.standard_atmosphere(z_grid)[2] / ch01.standard_atmosphere(0.0)[2]   # density relative to sea level [–]
z_tenth = np.interp(0.1, rho_ratio_z[::-1], z_grid[::-1])                  # height where ρ/ρ0 = 1/10 (ratio falls with z, so reverse)
print(f"air is ten times less dense than at sea level at about {z_tenth/1e3:.1f} km (geopotential height, USSA-1976)")
""")
nb.md("""
**What would change if…** the gas were ten times thinner (the height printed above)? n drops tenfold, so the noise line
rises by √10 and l grows tenfold: the window shrinks from both ends. **Next:** once properties exist at points, they can differ
from point to point — and molecules carry them across. That is §1.5 (C12).
""")

# =====================================================================================================================
# A.5 §1.5 Molecular Transport Phenomena — C12
# =====================================================================================================================
nb.section("1.5", "Molecular Transport Phenomena", intro="""
**What is this section about?** Random molecular motion carries salt, heat and momentum from where there is more to where
there is less. Three look-alike laws describe it; the one for momentum defines viscosity.
""")
core("C12", "Newton's law of viscosity: shear stress is a flux of momentum",
     "Hold one plate still and drag the other: how does the fluid in between find out, and what force does it exert on the plates?")
nb.md("""
#### The problem in plain words
Stir honey and it resists; stir water and it barely does. Oil in a bearing, wind dragging the ocean surface along, the
atmospheric boundary layer over a field (the lowest few hundred metres of air, slowed by the ground) — all are fluid layers
sliding over each other. We want the force per area between layers, and how fast motion spreads into still fluid.
""")
P("P19", "ordinary derivative as a slope", """
du/dy is the slope of the profile u(y): how much u changes per metre of y — for a velocity its unit is (m/s)/m = 1/s.
""", code="""
y2 = np.array([0.0, 1e-3])          # two heights 1 mm apart [m]
u2 = np.array([0.0, 1.0])           # velocities there [m/s]
print((u2[1] - u2[0]) / (y2[1] - y2[0]))   # slope du/dy = 1 m/s over 1 mm → 1000.0 1/s
""")
P("P25", "partial derivative", r"""
For a quantity that depends on several variables, such as f(y, t) or p(x, y, z), the **partial derivative** ∂f/∂y is its
slope in y with every other variable held fixed; if p depends on z alone, ∂p/∂z is just the ordinary dp/dz. Taking the
slope twice, **∂²f/∂y²** is the slope of the slope — the *curvature* of the profile: positive where the profile bends up
(a dip), negative where it bends down (a bump). The **gradient** ∇f = (∂f/∂x, ∂f/∂y, ∂f/∂z) collects the three slopes into
one vector that points where f increases fastest (Ch. 2 §2.9 treats it properly); "down the gradient" means along −∇f.
""", code="""
def f(y, t):                                   # a profile that also changes in time: f = y² e^(−t)
    return y**2 * np.exp(-t)
h_ = 1e-3                                      # a small step
print((f(1 + h_, 0) - f(1 - h_, 0)) / (2 * h_))                     # ∂f/∂y at y = 1, t = 0 (t held fixed) → 2.0
print((f(1 + h_, 0) - 2 * f(1, 0) + f(1 - h_, 0)) / h_**2)          # ∂²f/∂y²: slope of the slope → 2.0 (bends up)
print((f(1, h_) - f(1, -h_)) / (2 * h_))                            # ∂f/∂t at y = 1 (y held fixed) → −1.0
""")
nb.md("""
#### The idea
Molecules hop between neighbouring layers and carry their own layer's property with them:

```
y ↑   fast layer   u + du    ● → ● → ●      a molecule hopping DOWN brings extra x-momentum
      ─ ─ ─ ─ ─ ─ ─ ─ ─ A ─ ─ B ─ ─ ─ ─      a molecule hopping UP brings a deficit
      slow layer   u         ● → ●          net: x-momentum flows from fast to slow = a drag τ on the surface AB
```

| what spreads | its gradient | law | diffusivity [m²/s] |
|---|---|---|---|
| a species (mass fraction Y) | ∇Y | Fick (1.1) | κ_m |
| heat | ∇T | Fourier (1.2) | κ = k/ρC_p |
| x-momentum | du/dy | Newton (1.3) | ν = μ/ρ |

(∇ is the gradient of P25; k is the thermal conductivity and C_p the heat needed per kg per kelvin, both met properly in
C11 and §1.8.)

**Viscosity is diffusion of momentum.** In §1.3 (C02) a fluid kept deforming under any shear stress; the law below says
how fast.
""")
note("C08", "Transport down a gradient", r"""
Where a property (salt, heat, momentum) is unevenly spread, random molecular motion carries more of it from high to low
than back: a **flux** down the gradient, proportional to the gradient when the gradient is gentle. A profile that is not
straight therefore relaxes — its bumps fill in. For a profile f(y) this is the model diffusion equation with diffusivity
D [m²/s] (ours in Ch. 1; derived for heat in Ch. 4 §4.8). A bump 1 mm wide in water (D = ν ≈ 10⁻⁶ m²/s) smooths out in
about h²/D ≈ 1 s.
""", equation=r"\frac{\partial f}{\partial t} = D\,\frac{\partial^2 f}{\partial y^2}")
note("C09", "Mass fraction", """
The amount of a constituent is measured by its mass fraction Y (kg of it per kg of mixture; its partial density is ρY) —
salinity S in §1.10 is one, and salt and heat diffusing at different rates drive double-diffusive instability (Ch. 11 §11.5).
""")
note("C10", "Fick's law", r"""
The mass flux of a constituent points down its mass-fraction gradient, with mass diffusivity κ_m [m²/s] (*gloss:* ∇Y is
the vector of slopes (∂Y/∂x, ∂Y/∂y, ∂Y/∂z), pointing where Y increases fastest; Ch. 2 §2.9 treats it properly). Number:
water vapour in air, κ_m = 2.6×10⁻⁵ m²/s, ρ = 1.20 kg/m³, Y increasing upward by 0.01 per metre → J = −3.1×10⁻⁷ kg m⁻² s⁻¹,
i.e. a downward flux.
""", equation=r"\mathbf J_m = -\rho\,\kappa_m\,\nabla Y", ref="1.1")
note("N04", "The species picture", """
Drawn as a profile with a downward flux arrow, this is the "species" mode of the explainer at the end of this block.
""")
note("C11", "Fourier's law", r"""
Heat flux points down the temperature gradient, with thermal conductivity k [W m⁻¹ K⁻¹]. Divided by the heat stored per
kelvin per unit volume, ρC_p, it gives the thermal diffusivity κ = k/(ρC_p) [m²/s] (*gloss:* C_p, the heat needed per kg
per kelvin at constant pressure, is defined properly in §1.8). Number: water, k = 0.60 W m⁻¹ K⁻¹, temperature falling
upward by 100 K/m → q = +60 W/m² (upward); κ = 1.44×10⁻⁷ m²/s.
""", equation=r"\mathbf q = -k\,\nabla T", ref="1.2")
P("P20", "boundary conditions", """
Boundary conditions say what the solution must do at the edges. At a solid wall a viscous fluid moves with the wall — the
**no-slip** condition. **Steady** means nothing changes in time any more; **transient** is the approach to it.
""", code="""
u0 = np.zeros(5)      # initial velocity profile across a gap, 5 grid points [m/s]: fluid at rest
u0[-1] = 1.0          # no-slip at the top plate, which suddenly moves at 1 m/s
print(u0)             # [0. 0. 0. 0. 1.]
""")
nb.md(r"""
#### The maths, step by step
1. Layers slide in x with speed u(y) [m/s]; y [m] runs across the gap.
2. Across a surface AB of constant y, hopping molecules carry x-momentum per area per time — a **momentum flux**
   [kg m s⁻¹ per m² per s = kg m⁻¹ s⁻² = N/m² = Pa], i.e. a stress.
3. As in (1.1)–(1.2), the flux is proportional to the gradient du/dy. The constant is the **dynamic viscosity**
   μ [Pa s = kg m⁻¹ s⁻¹]:
   $$\tau = \mu\,\frac{du}{dy} \qquad (1.3)$$
   τ is the drag that the faster fluid above AB exerts on the fluid below it.
4. Units check: Pa s × s⁻¹ = Pa ✓.
5. Divide by the density to get a diffusivity, in the same units as κ_m and κ:
   $$\nu \equiv \mu/\rho \qquad (1.4) \quad [\mathrm{m^2/s}]$$
6. A diffusivity D spreads a disturbance across a distance h in a time of order h²/D. Why that combination: it is the
   only time you can build from h [m] and D [m²/s] (a preview of dimensional analysis, C64), and the animation below
   uses that clock: it runs h = 1 mm, but on the clock t/(h²/ν) a wider gap would give exactly the same frames.

Symbols: τ [Pa] shear stress, μ [Pa s], ν [m²/s], h [m] gap width, U [m/s] plate speed.
""")
note("C13", "Viscosity depends on temperature — in opposite ways", """
A **gas** gets more viscous when hot (faster molecules carry momentum farther; roughly μ ∝ T^(1/2)); a **liquid** gets
less viscous (the cohesive forces that resist sliding weaken). Numbers (code below): air 1.79×10⁻⁵ Pa s at 288 K falls to
1.42×10⁻⁵ at 217 K (Sutherland's law); the √T rule gives 1.55×10⁻⁵, 9 % high. Water: 1.00×10⁻³ Pa s at 20 °C falls to
3.5×10⁻⁴ at 80 °C.
""")
note("C14", "Kinematic viscosity", """
ν = μ/ρ is the diffusivity of momentum, so it — not μ — decides how fast motion spreads, and the ranking can flip: air's
μ is 55× smaller than water's, but its ν is 15× larger, so motion diffuses through still air 15× faster.
""", equation=r"\nu \equiv \mu/\rho", ref="1.4")
nb.worked_example("water between plates 1 mm apart", """
1. Top plate U = 1 m/s, bottom plate fixed, gap h = 1 mm: the steady slope is du/dy = U/h = 1000 s⁻¹.
2. Water, μ = 1.0×10⁻³ Pa s → τ = μU/h = 1.0 Pa.
3. On a 10 cm × 10 cm plate: F = τA = 1.0 × 0.01 = 0.01 N.
4. ν = μ/ρ = 1.0×10⁻³/998 = 1.0×10⁻⁶ m²/s → h²/ν = 10⁻⁶/10⁻⁶ = 1.0 s to set up the flow.
5. Air: μ = 1.81×10⁻⁵ → τ = 0.018 Pa (55× less), but ν = 1.51×10⁻⁵ → h²/ν = 0.066 s (15× faster).
""")
P("P23", "Python dictionaries", """
A dict maps names to values, `d["mu"]`. fluidpy returns several named results this way.
""", code="""
fp = {"rho": 998.2, "mu": 1.0e-3}     # a small dictionary: density [kg/m³] and viscosity [Pa s] of water
print(fp["mu"] / fp["rho"])           # look values up by name → ν = 1.0e-06 m²/s
""")
nb.code("""
h, U = 1e-3, 1.0                                          # gap width [m] and top-plate speed [m/s]
for name in ("water", "air"):                             # the two fluids of the worked example, at 20 °C
    fp = ch01.fluid_properties(name)                      # dict with rho [kg/m³], mu [Pa s], nu [m²/s] from cited correlations
    tau = ch01.newton_shear_stress(fp["mu"], U / h)       # Eq. (1.3): steady Couette stress τ = μ U/h [Pa]
    nu = ch01.kinematic_viscosity(fp["mu"], fp["rho"])    # Eq. (1.4): ν = μ/ρ [m²/s]
    t_d = ch01.diffusion_time(h, nu)                      # time for momentum to cross the gap, h²/ν [s]
    print(f"{name:5s}: τ = {tau:.4g} Pa, ν = {nu:.4g} m²/s, h²/ν = {t_d:.3g} s")   # one line per fluid
""", explain="""
1. `fluid_properties` returns density and viscosity at 20 °C (primer P23 for the dictionary).
2. `newton_shear_stress` is Eq. (1.3) with the steady slope U/h.
3. `kinematic_viscosity` is Eq. (1.4); `diffusion_time` is h²/ν.
4. Water: τ ≈ 1.00 Pa and 1 s; air: τ ≈ 0.018 Pa but only 0.066 s — the μ-versus-ν surprise in two lines.
""")
P("P21", "finite differences", r"""
A derivative from samples: the **central difference** $(f_{i+1} - f_{i-1})/(2\Delta y)$ is second-order accurate (its error
shrinks like Δy²). Stepping $\partial f/\partial t = D\,\partial^2 f/\partial y^2$ forward in time with
$f_i^{\rm new} = f_i + r\,(f_{i+1} - 2f_i + f_{i-1})$, $r = D\Delta t/\Delta y^2$, is the **FTCS** scheme (forward in time,
centred in space); it is stable only for r ≤ ½ (Ch. 10 §10.2 explains why).
""", code="""
d = 0.01                                               # sample spacing Δ
approx = (np.sin(0.5 + d) - np.sin(0.5 - d)) / (2 * d) # central difference of sin at 0.5
print(approx - np.cos(0.5))                            # error vs the exact cos(0.5): ≈ −1.5e-05 = −cos(0.5)·Δ²/6 (∝ Δ²)
""")
nb.code("""
N = 101 if not FAST else 51                               # grid points across the gap
y = np.linspace(0, h, N)                                  # y from the fixed plate (0) to the moving plate (h) [m]
dy = y[1] - y[0]                                          # grid spacing [m]
mu_w, rho_w = 1.0e-3, 998.2                               # water viscosity [Pa s] and density [kg/m³]
nu_w = ch01.kinematic_viscosity(mu_w, rho_w)              # Eq. (1.4) [m²/s]
dt = ch01.stable_time_step(nu_w, dy)                      # FTCS step with r = 0.9 × ½ (primer P21) [s]
nsteps = int(1.5 * h**2 / nu_w / dt)                      # march to t = 1.5 h²/ν
save_every = max(nsteps // 60, 1)                         # keep about 60 profiles for the animation
u0 = np.zeros(N); u0[-1] = U                              # fluid at rest, top plate suddenly moving (no-slip, P20)
F = ch01.ftcs_diffusion_1d(u0, nu_w, dy, dt, nsteps, values=(0.0, U), save_every=save_every)  # u(y, t), shape (nsave, N)
times = np.arange(F.shape[0]) * save_every * dt           # the time of every saved profile [s]
tau_b, tau_t = ch01.wall_shear_history(F, dy, mu_w)       # Eq. (1.3) at the bottom and top walls for every saved time [Pa]
u_exact = ch01.couette_startup_profile(y, times[-1], U, h, nu_w)   # analytic series (derived in Ch. 8; used here only as a check)
print(f"{nsteps} steps, {F.shape[0]} saved profiles, last t = {times[-1]:.2f} s")   # size of the run
print(f"wall stresses at the end: bottom {tau_b[-1]:.3f} Pa, top {tau_t[-1]:.3f} Pa (steady value {mu_w*U/h:.3f})")   # both walls reach μU/h
print(f"largest difference from the series: {np.abs(F[-1] - u_exact).max():.1e} m/s")   # numerical vs exact
""", explain="""
1. A grid across the gap and a time step that satisfies the FTCS limit (primer P21) — `stable_time_step`.
2. The initial profile: still water with the top plate suddenly at U (boundary conditions, primer P20).
3. `ftcs_diffusion_1d` marches ∂u/∂t = ν∂²u/∂y² (the C08 diffusion picture with D = ν), saving about 60 profiles.
4. `wall_shear_history` applies (1.3) at both walls for every saved profile: both stresses approach μU/h = 1 Pa.
5. `couette_startup_profile` is the exact series solution of the same problem (derived in Ch. 8); the numerical profile
   agrees with it.
""")
P("P22", "np.gradient", """
`np.gradient(f, y, edge_order=2)` returns central differences inside the array and second-order one-sided differences at
the two ends. The default `edge_order=1` is only first-order accurate at the ends.
""", code="""
print(np.gradient(np.array([0.0, 1.0, 4.0, 9.0]), 1.0, edge_order=2))   # slopes of y² at y = 0,1,2,3 → [0. 2. 4. 6.], exact
""")
nb.md("""
#### From scratch: the stress profile by hand
`shear_stress_profile` applies (1.3) to a sampled profile. Here is the same thing written out as a loop, compared with the
library and with `np.gradient`.
""")
nb.check_agree("""
u_mid = ch01.couette_startup_profile(y, 0.1 * h**2 / nu_w, U, h, nu_w)   # a curved profile early in the start-up [m/s]
tau_mine = np.empty_like(u_mid)                                          # our stress profile [Pa]
for i in range(1, N - 1):                                                # interior points: central difference × μ
    tau_mine[i] = mu_w * (u_mid[i + 1] - u_mid[i - 1]) / (2 * dy)   # τ = μ du/dy with a central difference
tau_mine[0] = mu_w * (-3 * u_mid[0] + 4 * u_mid[1] - u_mid[2]) / (2 * dy)     # bottom wall: one-sided, 2nd order
tau_mine[-1] = mu_w * (3 * u_mid[-1] - 4 * u_mid[-2] + u_mid[-3]) / (2 * dy)  # top wall: one-sided, 2nd order
assert np.allclose(tau_mine, ch01.shear_stress_profile(u_mid, y, mu_w))       # same numbers as the library
assert np.allclose(tau_mine, mu_w * np.gradient(u_mid, y, edge_order=2))      # and as numpy's gradient (P22)
print("stress profile: hand loop = library = np.gradient")   # reached only if both asserts passed
fmt = lambda a: ", ".join(f"{v:.3g}" for v in np.atleast_1d(a))                 # print small numbers with 3 significant figures
print("C13 air μ at 288.15, 216.65 K (Sutherland):", fmt(ch01.sutherland_viscosity(np.array([288.15, 216.65]))), "Pa s")   # gas: μ falls as T falls
print("C13 air μ at 216.65 K from the √T rule:", fmt(ch01.viscosity_power_law(216.65, ch01.sutherland_viscosity(288.15), 288.15)), "Pa s")   # the rough rule, ≈ 9 % high
print("C13 water μ at 20 °C and 80 °C:", fmt(ch01.water_viscosity(np.array([293.15, 353.15]))), "Pa s")   # falls when heated
print("C10 Fick flux (x, y, z):", fmt(ch01.fick_mass_flux(1.20, 2.6e-5, np.array([0.0, 0.01, 0.0]))), "kg m⁻² s⁻¹")  # downward
print("C11 Fourier flux (x, y, z):", fmt(ch01.fourier_heat_flux(0.60, np.array([0.0, -100.0, 0.0]))), "W/m²")        # upward
print("C11 thermal diffusivity of water:", fmt(ch01.thermal_diffusivity(0.60, 998.0, 4182.0)), "m²/s")              # κ = k/(ρ C_p)
""")
nb.animation("""
frames_c12 = 40 if not FAST else 20                                      # frames spread over the saved profiles
idx = np.linspace(0, F.shape[0] - 1, frames_c12).astype(int)             # which saved profiles to show
fig, (axu, axt) = plt.subplots(1, 2, figsize=(8, 3.4), sharey=True)      # left u(y), right τ(y); same y axis
axu.plot([0, 1], [0, 1], "--", color=COLORS["muted"], lw=1, label="steady Couette")   # the final straight profile
(line_u,) = axu.plot([], [], color=COLORS["teal"], lw=2.5, label="FTCS u(y, t)")      # the numerical profile
(line_s,) = axu.plot([], [], color=COLORS["orange"], lw=1, label="series solution")   # the exact check (Ch. 8)
axu.set_xlim(-0.02, 1.02); axu.set_ylim(0, 1)                            # fixed axes
axu.set_xlabel("u / U [–]"); axu.set_ylabel("y / h [–]"); axu.legend(loc="upper left", fontsize=8)   # dimensionless axes, legend
axt.axvline(mu_w * U / h, color=COLORS["muted"], ls="--", lw=1)          # steady stress μU/h = 1 Pa
(line_t,) = axt.plot([], [], color=COLORS["rose"], lw=2.5)               # τ(y, t) from (1.3)
axt.set_xlim(0, 6); axt.set_xlabel("shear stress τ [Pa] (clipped at 6)")   # fixed stress axis
title = fig.suptitle("")                                                  # shows the time in units of h²/ν
def update(k):                                                            # frame k shows saved profile idx[k]
    j = idx[k]                                                            # index of the saved profile
    line_u.set_data(F[j] / U, y / h)                                      # numerical profile
    line_s.set_data(ch01.couette_startup_profile(y, max(times[j], 1e-6), U, h, nu_w) / U, y / h)   # series at that time
    line_t.set_data(ch01.shear_stress_profile(F[j], y, mu_w), y / h)      # stress profile, Eq. (1.3)
    title.set_text(f"t = {times[j]*nu_w/h**2:.2f} h²/ν = {times[j]:.2f} s (water, h = 1 mm)")   # time in units of h²/ν and in seconds
    return line_u, line_s, line_t, title   # the artists that changed
show_animation(animate(update, frames=frames_c12, fig=fig, interval=80), player="video")   # smooth MP4
""", explain="""
1. Each frame takes one saved FTCS profile: the teal curve is u(y, t)/U, the thin orange curve the exact series at the
   same time, the grey dashed line the final straight (Couette) profile.
2. The right panel applies (1.3) to the same profile with `shear_stress_profile`; the grey line is the steady μU/h.
3. The title gives the time in units of h²/ν, the natural clock of diffusion.
""")
nb.md("""
**What you see.** A boundary layer creeps down from the moving plate into still water; the stress, huge at the top wall
at first, decays and spreads until it is the same at every height.

**How to read it.** A uniform τ means every layer passes on exactly the momentum it receives — the steady Couette flow.
By t ≈ 0.5 h²/ν the profile is within about 1 % of straight.

**What would change if…** the gap were twice as wide: everything would happen 4× more slowly in seconds, but identically
in units of h²/ν — which is why the title uses that clock.
""")
note("N05", "The relaxation picture", """
The animation is our version of the book's relaxation sketch (Fig. 1.3): a sheared profile relaxing toward the straight
one, with the stress on any surface AB equal to μ times the local slope.
""")
note("N06", "Linear, first-derivative laws", """
All three transport laws are **linear** in the gradient and use **first** derivatives only — accurate because molecular
steps are tiny compared with the scale of the gradients; Ch. 4 §4.5 generalises (1.3) to a tensor law, and the diffusion
equations appear in Ch. 4 §4.8 and Ch. 8 §8.4.
""")
nb.explainer("viscosity_momentum_diffusion", heading="How does the fluid learn that a plate moved?", why="""
The phenomenon is a transient: a profile creeping into the gap and settling. Playing, scrubbing and swapping fluids shows
h²/ν at work and the μ-versus-ν surprise, which no single frame can.
""", tries=[
    "Play with water, then with air, and compare the settling times on the τ_w(t) panel.",
    "Pick honey with a 5 mm gap and predict first: faster or slower than water?",
    "Switch to heat mode and read the Prandtl-number comparison in the Explain tab (Pr = ν/κ: how much faster momentum "
    "diffuses than heat; ≈ 7 for water, ≈ 0.7 for air).",
    "Click the profile to see the τ(y) arithmetic at that height.",
])
nb.md("""
**What would change if…** the top plate stopped again? Momentum would diffuse out through both walls and the profile
would decay on the same clock h²/ν. **Next:** when nothing moves at all, every shear stress vanishes and only the normal
stress — pressure — remains: fluid statics (§1.7, C20). First, one more force that acts at interfaces (§1.6).
""")

# =====================================================================================================================
# A.6 §1.6 Surface Tension
# =====================================================================================================================
nb.section("1.6", "Surface Tension", intro="""
**What is this section about?** An interface between two fluids pulls like a stretched sheet; where it is curved, the
pressure on its concave side is higher. The C20 block of §1.7 uses this for capillary rise.
""")
note("C15", "Surface tension σ", """
Surface tension is the pull per unit length along any line drawn in the interface [N/m] — equivalently the energy needed
to create one extra square metre of interface [J/m² = N/m]. It depends on the pair of fluids, the temperature and traces
of surfactants (soap lowers it). Water against air at 20 °C: σ = 72.7 mN/m (IAPWS).
""")
note("N07", "The spherical drop", r"""
For a drop or bubble of radius R, the pull σ × 2πR around a cut through its middle balances the pressure excess acting on
the cut's area πR² (stated, not derived). A 1 mm drop of water has 145 Pa extra inside; a 1 µm fog droplet has
1.45×10⁵ Pa — about 1.4 atmospheres.
""", equation=r"p_i - p_o = \frac{2\sigma}{R}")
note("C16", "Laplace's pressure jump", r"""
For any gently curved interface the jump uses its two principal radii of curvature (*gloss:* at a point on a surface, the
circles that best fit it in two perpendicular directions; their radii R₁ and R₂ are taken negative when the centre lies on
the other side, as on a saddle). The concave side has the higher pressure. Cylinder (R₂ → ∞): σ/R. Sphere (R₁ = R₂): 2σ/R.
Symmetric saddle (R₂ = −R₁): no jump at all.
""", equation=r"p_i - p_o = \sigma\left(\frac{1}{R_1} + \frac{1}{R_2}\right)", ref="1.5")
nb.code("""
sigma = ch01.surface_tension_water(293.15)                          # C15: water–air surface tension at 20 °C [N/m]
print(f"σ = {sigma:.5f} N/m")                                       # ≈ 0.07274 N/m
print(ch01.laplace_pressure_jump(sigma, np.array([1e-3, 1e-6])))    # N07: spheres of 1 mm and 1 µm → [145.5, 1.455e5] Pa
print(ch01.laplace_pressure_jump(sigma, 1e-3, np.inf))              # C16: a 1 mm cylinder → 72.7 Pa
print(ch01.laplace_pressure_jump(sigma, 1e-3, -1e-3))               # C16: a symmetric saddle → 0.0 Pa
""", explain="""
1. `surface_tension_water` evaluates the IAPWS correlation at 293.15 K.
2. `laplace_pressure_jump` with one radius is a sphere (both radii equal): the N07 numbers.
3. An infinite second radius is a cylinder; opposite radii are a saddle, whose two curvatures cancel.
""")
note("N08", "Capillarity", """
Narrow tubes where these jumps matter are **capillary** tubes, and the phenomena **capillarity**; surface tension becomes a
boundary condition in Ch. 4 §4.10 and drives capillary waves in Ch. 7 §7.3.
""")
note("N09", "Principal radii in pictures", """
The book sketches the hemisphere force balance and a surface patch with its two principal radii; those radii return with
capillary waves (Ch. 7 §7.3).
""")

# =====================================================================================================================
# A.7 §1.7 Fluid Statics — C20
# =====================================================================================================================
nb.section("1.7", "Fluid Statics", intro="""
**What is this section about?** Pressure in a fluid at rest: it acts equally in all directions, does not change sideways,
and grows downward by exactly the weight of the fluid above.
""")
core("C20", "The hydrostatic law: pressure grows with depth by the weight above",
     "A cube of lake water sits perfectly still. What must the pressure do between its top and bottom faces so that nothing moves?")
nb.md("""
#### The problem in plain words
Your ears hurt at the bottom of a pool; a dam is thicker at its base; a barometer weighs the atmosphere; the deep ocean
crushes submarines. The base state of every stratified flow, wave and weather system in this book is this one balance.
""")
nb.md("""
#### The idea
```
            (p + dp)·dx dy   ↓        the top face is pushed down
          ┌────────────┐
          │  weight    │  ρ g dx dy dz ↓
          └────────────┘
            p·dx dy          ↑        the bottom face is pushed up
at rest:  (push up on the bottom) − (push down on the top) = weight   ⇒   pressure must fall going up
```
""")
nb.figure("""
from scripts.ch01_drawings import draw_cube_forces      # our drawing of the fluid cube (geometry only, no physics)
fig, ax = plt.subplots(figsize=(4.2, 3.6))               # a small square figure
draw_cube_forces(ax)                                     # the cube, its top and bottom pressure forces and its weight
ax.set_title("Forces on a still fluid cube (z up)")      # the message
plt.show()                                               # display
""", see="""
A small cube of fluid with two orange pressure forces — up on the bottom face, down on the top face — and its rose weight
pointing down.
""", read="""
The two orange arrows are not equal: the difference between them is exactly what holds up the cube's weight. The
derivation D05 below makes this exact.
""", change="""
In a denser fluid (mercury instead of water) the same cube is 13.6× heavier, so the pressure difference between its faces
must be 13.6× larger.
""")
note("N11", "The cube balance", """
The pressure difference between a fluid element's top and bottom faces balances its weight — the picture above, made
exact in D05.
""")
P("P24", "weight and gravitational acceleration", """
Weight W = mg points down; standard gravity is g = 9.80665 m/s². For fluid of density ρ filling a volume V, W = ρVg.
(Newton's second law was primed in C06, P09.)
""", code="""
rho_water, V_litre = 1000.0, 1e-3              # density of water [kg/m³] and one litre [m³]
print(rho_water * V_litre * 9.80665)           # its weight W = ρ V g → 9.81 N
""")
nb.md("""
> 🔁 Partial derivatives ∂p/∂x (slope in one variable with the others held fixed) were primed in C12 (P25); below they
> tell us how pressure may vary sideways and upward in still water.
""")
P("P26", "first-order Taylor expansion", r"""
Close to a point, a smooth function is its value plus slope × step: $p(z + dz) \approx p(z) + (dp/dz)\,dz$. The neglected
terms shrink like dz², so after dividing by dz they vanish as dz → 0 — they are of a higher **order of smallness**.
""", code="""
step = 0.1                                   # a small step
print(np.exp(step), 1 + step)                # e^0.1 = 1.10517 vs the first-order estimate 1.1
print(np.exp(step) - (1 + step), step**2/2)  # the error 0.00517 ≈ step²/2: second order
""")
note("C17", "Absolute and gauge pressure", r"""
In a fluid at rest the normal stress is the **pressure** — the compressive normal stress of C03 (§1.3). **Absolute**
pressure is measured from vacuum; **gauge** pressure from the local atmosphere. The standard atmosphere is 101 325 Pa and
1 bar = 10⁵ Pa. A tyre at 2.0 bar gauge holds 3.01 bar absolute.
""", equation=r"p_{\rm gauge} = p - p_{\rm atm}")
note("C18", "Pressure is the same in every direction", r"""
(Stated, not derived.) Take a tiny wedge of fluid at rest. The forces on its faces scale with face area (∝ size²) but its
weight with its volume (∝ size³), so as the wedge shrinks the weight becomes negligible and the face pressures must be
equal, whatever the wedge's orientation. For a 1 mm tall wedge of water the pressure difference is only ½ρg dz ≈ 4.9 Pa,
and it vanishes as dz → 0 (printed in the code below). So pressure at a point is a single number, a scalar.
""", equation=r"p_1 = p_2 = p_3", ref="1.6")
note("N10", "Towards the stress tensor", """
This direction-independent normal stress becomes the −pδᵢⱼ part of the stress tensor (Ch. 2 §2.6, Ch. 4 §4.5).
""")
note("C19", "Pascal's law", r"""
(Stated.) On a fluid cube at rest the two faces normal to x feel p(x) and p(x + dx); nothing else acts sideways (gravity
is vertical), so they must be equal and ∂p/∂x = 0 — likewise in y. Points at the same height in one connected, resting
fluid have the same pressure, which is why water levels in connected vessels line up.
""", equation=r"\frac{\partial p}{\partial x} = \frac{\partial p}{\partial y} = 0", ref="1.7")
D("D05", "The hydrostatic law", ref="1.8",
  goal="""Find how pressure must change with height inside a fluid that is not moving, so that every little piece of it
stays at rest. This is the base state of everything stratified in the book.""",
  assumptions="""Fluid at rest (step 6: a = 0; step 4: p = p(z) only) · z points up (steps 1–2: signs of the face forces) ·
g uniform over the tiny cube (step 5) · dz small enough for a first-order Taylor step (step 3).""",
  start=(r"\sum F_z = m\,a_z = 0",
         "a small cube of fluid at rest has zero acceleration, so the vertical forces on it add to zero (Newton's second law, P09)."),
  plan=["List the vertical forces on a cube dx × dy × dz: pressure on the bottom, pressure on the top, weight.",
        "Relate the top pressure to the bottom one with a Taylor step.",
        "Add the forces, cancel, divide by the volume."],
  uses=["Newton's second law (P09, C06)", "weight W = mg (P24)",
        "pressure acts normal to a surface and is the same in every direction (C18)",
        "Pascal's law, p independent of x and y at rest (C19)", "partial derivative (P25)",
        "first-order Taylor expansion (P26)"],
  steps=[
      ("Write the push on the bottom face", r"F_{\rm bottom} = +\,p(z)\,dx\,dy",
       "Pressure pushes into the fluid normal to a face; on the bottom face that is upward (positive z). We start the force "
       "list with it because the cube is held up from below.",
       "The fluid underneath pushes the cube up with pressure times area."),
      ("Write the push on the top face", r"F_{\rm top} = -\,p(z+dz)\,dx\,dy",
       "On the top face the fluid above pushes down, so the force is negative with z up; it uses the pressure at the top "
       "face's height z + dz.",
       "The fluid above pushes the cube down with the pressure found one step higher."),
      ("Expand the top pressure to first order", r"p(z+dz) \approx p(z) + \dfrac{\partial p}{\partial z}\,dz",
       "First-order Taylor expansion (P26), an approximation (≈): the dropped terms are proportional to dz², so they vanish "
       "when we divide by the volume and let the cube shrink, and the relation becomes exact as dz → 0. We need the two face "
       "pressures in terms of one value.",
       "The top pressure is the bottom pressure plus slope times height."),
      ("Use Pascal's law to make it an ordinary derivative", r"\dfrac{\partial p}{\partial z} = \dfrac{dp}{dz}",
       "At rest p does not change with x or y (C19, Eq. 1.7), so p depends on z alone and its partial derivative in z is "
       "the ordinary one. The book writes dp directly; this is why it may.",
       "In still fluid pressure depends only on height."),
      ("Write the weight", r"W = -\,\rho g\,dx\,dy\,dz",
       "Mass = density × volume and weight = mass × g, pointing down (P24); ρ and g are treated as constant across the "
       "tiny cube.",
       "The cube's weight pulls it down."),
      ("Add the three forces and set the sum to zero",
       r"p\,dx\,dy - \Big(p + \dfrac{dp}{dz}dz\Big)dx\,dy - \rho g\,dx\,dy\,dz = 0",
       "Newton's second law with zero acceleration (the fluid is at rest), using steps 1–5. The horizontal forces already "
       "cancel by Pascal's law, so only this vertical balance is left.",
       "Push up, push down and weight exactly cancel."),
      ("Cancel the p dx dy terms", r"-\dfrac{dp}{dz}\,dx\,dy\,dz - \rho g\,dx\,dy\,dz = 0",
       "The term p dx dy appears once with a plus and once with a minus sign, so the two cancel; what is left are the "
       "parts that differ between the faces, plus the weight.",
       "Only the pressure difference between the faces and the weight remain."),
      ("Divide by the volume dx dy dz", r"-\dfrac{dp}{dz} - \rho g = 0",
       "Every remaining term contains the volume dx dy dz ≠ 0, so we may divide by it. This removes the arbitrary size of "
       "the cube, and the dropped Taylor terms (∝ dz) vanish as the cube shrinks.",
       "Per unit volume, the upward pressure-gradient force balances the weight."),
      ("Move the weight term across", r"\dfrac{dp}{dz} = -\rho g",
       "Add ρg to both sides, then multiply both sides by −1 (an equation, not an inequality, so nothing else changes). We "
       "isolate dp/dz because it is what we want to know.",
       "Pressure falls with height by ρg per metre."),
  ],
  result=(r"\dfrac{dp}{dz} = -\rho g",
          "going up one metre in still fluid, pressure drops by the weight of a one-metre column of unit cross-section."),
  interpret="""Pressure in still fluid is the weight of everything above, per unit area. It holds for liquids and gases
alike — ρ may vary with z or with p, and then we integrate numerically, as `integrate_hydrostatic` does below. It fails
when the fluid accelerates (Ch. 4 adds ρ Du/Dt) or when g varies over the column (planetary scales).""",
  check="""Units: Pa/m = kg m⁻² s⁻² and ρg = kg m⁻³ · m s⁻² = kg m⁻² s⁻² ✓. Limit: ρ → 0 (vacuum) gives constant p ✓.
Number: water, ρg = 1000 × 9.81 = 9810 Pa/m, so 10 m of depth adds 9.81×10⁴ Pa ≈ 1 atm ✓ (the code below prints
`hydrostatic_pressure_uniform` at −10 m).""",
  traps="with z pointing *down* (depth) the sign flips, dp/dh = +ρg; forgetting that the top-face force points down "
        "gives dp/dz = +ρg; writing dp/dz before knowing that p depends on z only.")
P("P27", "definite integral", r"""
$\int_a^b f\,dz$ adds up the pieces f dz from z = a to z = b; for a constant, $\int_a^b c\,dz = c\,(b - a)$. The trapezoid
rule `np.trapezoid(f, z)` approximates it from samples (exact for straight lines).
""", code="""
zs = np.linspace(-10, 0, 11)                 # heights from −10 m to 0 m
print(np.trapezoid(np.full(11, 9807.0), zs)) # ∫ ρg dz of a constant 9807 Pa/m over 10 m → 98070.0 Pa
""")
note("C21", "Uniform density", r"""
(Stated.) For uniform density, integrating (1.8) from z = 0, where the pressure is p₀, up or down to z gives a straight
line: pressure grows by ρgh at a depth h = −z. Ten metres of water add 9.81×10⁴ Pa ≈ one atmosphere. And because the
pressure pushes harder on the bottom of a submerged body than on its top, the body feels an upward net force —
**buoyancy**, derived next.
""", equation=r"p = p_0 - \rho g z", ref="1.9")
P("P28", "net force from pressure", r"""
Pressure p pushes on each bit of surface dA along the inward normal, giving a force $-p\,\hat{\mathbf n}\,dA$ with
$\hat{\mathbf n}$ pointing out of the body. The net force on a closed body is the sum over its whole surface,
$\mathbf F = -\oint p\,\hat{\mathbf n}\,dA$ — for a box, just six faces.
""", code="""
F_uniform = ch01.net_pressure_force_on_box(lambda x, y, z: 1e5 + 0 * z, (0, 1, 0, 1, 0, 1))  # uniform 1e5 Pa on a 1 m cube
print(F_uniform)                                                                             # → [0. 0. 0.] N: pushes cancel
""")
D("D37", "Buoyancy: the net pressure force on a submerged body (from Eq. 1.9; our addition)", ref="",
  goal="""Show that still fluid pushes up on a submerged body with a force equal to the weight of the fluid the body
displaces — whatever the body is made of. The parcel argument (D18, §1.10) needs this.""",
  assumptions="""Fluid at rest with uniform density near the body (step 5) · the body is fully submerged and does not
disturb the fluid's pressure (steps 2–9) · g uniform.""",
  start=(r"p(z) = p_0 - \rho g z", "in fluid of uniform density ρ the pressure grows linearly with depth (1.9)."),
  plan=["Put a box of horizontal area A between heights z₁ (bottom) and z₂ (top) into the fluid.",
        "Show that the side forces cancel.", "Subtract the top force from the bottom force and use (1.9).",
        "Build any shape from thin boxes."],
  uses=["(1.9) (C21)", "Pascal's law (C19)", "net force from pressure as a surface sum (P28)", "definite integral (P27)"],
  steps=[
      ("Pair up the side faces", r"F_{x,\rm left} + F_{x,\rm right} = p\,A_x - p\,A_x = 0",
       "Opposite side faces sit at the same heights, so by Pascal's law (C19) they feel the same pressure at each height, "
       "pushing in opposite directions. We clear the horizontal forces first to be left with a one-dimensional balance.",
       "Sideways pushes cancel."),
      ("Write the upward force on the bottom face", r"F_{\rm bottom} = p(z_1)\,A",
       "Pressure acts normal to the face, which on a bottom face means upward, over the face's area A.",
       "The fluid below pushes up."),
      ("Write the downward force on the top face", r"F_{\rm top} = p(z_2)\,A",
       "The same rule on the top face, where the normal push points down; the top is higher, z₂ > z₁.",
       "The fluid above pushes down, less hard because it is higher."),
      ("Take the net upward force", r"F_{\rm net} = \big[p(z_1) - p(z_2)\big]\,A",
       "Up minus down: this is the vertical part of the surface sum of P28, since the sides contribute nothing (step 1).",
       "What is left is the difference of the two face pressures."),
      ("Substitute the hydrostatic pressure", r"p(z_1) - p(z_2) = \rho g\,(z_2 - z_1)",
       "From (1.9): (p₀ − ρgz₁) − (p₀ − ρgz₂); p₀ cancels. This is where the fluid's density enters — and only the fluid's.",
       "The pressure difference is ρg times the box's height."),
      ("Recognise the volume", r"F_{\rm net} = \rho g\,A\,(z_2 - z_1) = \rho\,g\,V",
       "A (z₂ − z₁) is the volume V of the box; multiplying area by height is the definition of a box's volume.",
       "The net push up equals the weight of fluid that would fill the box."),
      ("Project a slanted piece of surface onto the horizontal", r"dF_z = p\,|n_z|\,dA_s = p\,dA",
       "A curved body has slanted top and bottom surfaces. On a small slanted piece of area dA_s the pressure force is "
       "p dA_s along the normal; its vertical part is p |n_z| dA_s, and |n_z| dA_s is exactly the piece's shadow dA on a "
       "horizontal plane. So for vertical forces a slanted cap acts like a flat face of area dA.",
       "Only the horizontal footprint of a surface matters for its vertical pressure force."),
      ("Apply the box result to one thin vertical column", r"dF_{\rm net} = \rho g\,h(x,y)\,dA",
       "Slice the body into vertical columns of footprint dA and height h(x, y). By step 7 its caps act like the flat faces of "
       "a box, and steps 4–6 give the net upward force ρg × height × area. Side forces between neighbouring columns cancel "
       "in pairs, as in step 1.",
       "Each thin column is buoyed up by the weight of fluid that would fill it."),
      ("Add all the columns", r"F_{\rm net} = \displaystyle\int_A \rho g\,h(x,y)\,dA = \rho g V",
       "Summing infinitely many thin columns is the definite integral of P27 over the body's footprint A; ρg is a constant "
       "and comes outside, and ∫h dA is the body's volume V.",
       "Any submerged body feels an upward force equal to the weight of the displaced fluid."),
  ],
  result=(r"F_b = \rho_{\rm fluid}\,g\,V", "Archimedes' principle, derived from the hydrostatic pressure."),
  interpret="""Buoyancy is nothing but the pressure difference between a body's bottom and top. It does not care what the
body is made of; the body's own weight decides whether it rises or sinks. In a stratified fluid (D18) the relevant density
is the environment's at the body's current height. It fails if the body sits on the bottom (no fluid underneath to push)
or is only partly submerged (then use the submerged volume).""",
  check="""Units: kg m⁻³ · m s⁻² · m³ = N ✓. Special cases: ρ_fluid → 0 gives no buoyancy ✓; a body with the fluid's own
density feels weight = buoyancy and is neutral ✓. Number: 1 L in water → 9.81 N; `net_pressure_force_on_box` integrates the
six faces numerically and returns (0, 0, 9.807) N ✓ (code below).""",
  traps="thinking the side forces add up; putting the body's density into the buoyancy; forgetting that the net force "
        "points up.")
note("C03", "Reminder: liquids under tension (from §1.3)", """
A liquid cannot be pulled below its vapour pressure: a sealed water column hanging below a suction pump cannot be lifted
more than about 10 m, because by (1.9) the pressure at its top would have to fall below ≈ 2.3 kPa.
""")
note("C22", "Capillary rise (Example 1.1, stated)", r"""
In a thin tube dipped into water the curved meniscus (C15, C16 in §1.6) lowers the pressure just under it, so water rises
until its weight balances the surface-tension pull around the rim: σ·2πR·sin α = ρgh·πR². **The book's α is measured from
the horizontal** — the complement of the usual contact angle θ_c (α = 90° − θ_c; α = 90° is perfect wetting;
`ch01.alpha_from_contact_angle` converts). Clean glass and water, R = 1 mm: h = 14.9 mm; R = 0.1 mm: 149 mm.
""", equation=r"h = \frac{2\sigma\sin\alpha}{\rho g R}")
note("N12", "Pressure along the tube axis", """
Along the axis the pressure falls linearly from atmospheric at the outside water level to p_atm − ρgh just under the
meniscus (inset of the p(z) figure below, where F marks the axis at the outside water level and E the point just under the meniscus); pressure jumps at free surfaces return in Ch. 7 §7.3.
""")
nb.worked_example("10 m of water and a 1-litre block", """
(g ≈ 9.81 m/s²)
1. p(−10 m) = 101 325 + 1000 × 9.81 × 10 = 199 425 Pa absolute; gauge 98 100 Pa.
2. A 1-litre block: buoyancy ρ_water g V = 1000 × 9.81 × 0.001 = 9.81 N, whatever the block is made of.
3. Steel (7850 kg/m³) weighs 77.0 N → net 67.2 N down: it sinks.
4. Pine (600 kg/m³) weighs 5.89 N → net 3.92 N up: it floats.
5. 1 m of oil (800 kg/m³) over 9 m of water: p = 101 325 + 800 × 9.81 × 1 + 1000 × 9.81 × 9 = 197 463 Pa
   (197 430 Pa with g = 9.80665).
""")
P("P29", "functions as arguments and lambda", """
Python functions are values: you can pass one to another function. `lambda z, p: 1000.0` is a one-line function (here a
density that ignores its inputs). `integrate_hydrostatic(z, rho_fn, p0)` takes the density as such a function.
""", code="""
apply = lambda f, x: f(x)          # a function that calls whatever function it is given
print(apply(lambda z: 2 * z, 3.0)) # pass a doubling function → 6.0
""")
P("P30", "explicit stepping", r"""
To follow dy/dx = f(x, y), take small steps: $y_{k+1} = y_k + f(x_k, y_k)\,\Delta x$ (**Euler's method**; its error shrinks
like Δx). For oscillators, update the velocity first and use the *new* velocity to move the position (**Euler–Cromer**;
keeps the energy bounded).
""", code="""
yk, dx = 1.0, 1e-3                  # start y(0) = 1, step 0.001
for _ in range(1000):               # 1000 Euler steps of dy/dx = −y to x = 1
    yk = yk + (-yk) * dx            # y_{k+1} = y_k + f(y_k) Δx
print(yk, np.exp(-1))               # 0.3677 vs the exact e^(−1) = 0.3679
""")
P("P31", "scipy.integrate.solve_ivp", """
`integrate_hydrostatic` hands dp/dz = −ρ(z, p) g to scipy's adaptive Runge–Kutta solver `solve_ivp`, which chooses its own
step sizes to reach a requested tolerance (here 10⁻¹⁰).
""", code="""
from scipy.integrate import solve_ivp                                              # scipy's ODE solver
print(solve_ivp(lambda x, yy: -yy, (0, 1), [1.0], rtol=1e-10, atol=1e-12).y[0, -1])  # dy/dx = −y to x = 1 → 0.367879
""")
nb.code("""
z = np.linspace(0, -10, 201)                                         # heights from the surface to 10 m down [m] (z up)
p_w = ch01.hydrostatic_pressure_uniform(z, ch01.P_ATM, 1000.0)       # Eq. (1.9): uniform water [Pa]
p_layers = ch01.layered_pressure(z, [1.0, 9.0], [800.0, 1000.0])     # 1 m of oil over 9 m of water, (1.9) layer by layer [Pa]
rho_fn = lambda zz, pp: 1000.0 - 0.5 * zz                            # a salty lake: 0.5 kg/m³ denser per metre down (P29)
p_strat = ch01.integrate_hydrostatic(z, rho_fn, ch01.P_ATM)          # Eq. (1.8) integrated numerically (solve_ivp, P31) [Pa]
p_exact = ch01.P_ATM + ch01.G0 * (1000.0 * (-z) + 0.25 * z**2)       # the same integral done by hand for this ρ(z) [Pa]
assert np.allclose(p_strat, p_exact, rtol=1e-10)                     # the integrator reproduces the exact answer
print(f"10 m down: uniform {p_w[-1]:.1f} Pa, oil over water {p_layers[-1]:.1f} Pa, salty lake {p_strat[-1]:.1f} Pa")   # the three profiles at the bottom
print(f"gauge pressure 10 m down in water: {ch01.gauge_pressure(p_w[-1]):.1f} Pa")             # C17
box = (0, 0.1, 0, 0.1, -1.1, -1.0)                                   # a 10 cm cube 1.0–1.1 m below the surface [m]
F_box = ch01.net_pressure_force_on_box(lambda xx, yy, zz: ch01.hydrostatic_pressure_uniform(zz, ch01.P_ATM, 1000.0), box)   # −∮p n dA over the cube's six faces [N]
print(f"net pressure force on the cube: {F_box} N;  ρ g V = {ch01.buoyancy_force(1000.0, 1e-3):.3f} N")   # D37 check
print("C18 wedge, 1 mm tall:", ch01.wedge_pressure_difference(1000.0, 1e-3, 0.6))              # differences of ≈ 5 Pa
h_cap = ch01.capillary_rise(ch01.surface_tension_water(293.15), np.pi / 2, 998.2, np.array([1e-3, 1e-4]))  # C22
print(f"capillary rise for R = 1 mm and 0.1 mm: {h_cap} m")                                     # [0.0149 0.149] m
""", explain="""
1. `hydrostatic_pressure_uniform` is (1.9); `layered_pressure` applies it layer by layer with a kink at the oil–water
   interface.
2. `integrate_hydrostatic` integrates (1.8) for a density that depends on z (passed as a lambda, P29) with `solve_ivp`
   (P31); the assert compares it with the integral done by hand.
3. `gauge_pressure` subtracts the atmosphere (C17).
4. `net_pressure_force_on_box` adds −p n̂ dA over the six faces of a submerged cube (P28): the horizontal parts cancel and
   the vertical part equals `buoyancy_force` = ρgV — derivation D37 checked numerically.
5. `wedge_pressure_difference` gives the tiny face-pressure differences of a 1 mm wedge (C18); `capillary_rise` gives
   the C22 heights.
""")
nb.md("""
#### From scratch: march down the salty lake
Euler's method (P30) applied to (1.8): start at the surface and step down 1 mm at a time, adding ρ(z)g|Δz| each step.
""")
nb.check_agree("""
dz_step = -1e-3                                         # step downward, 1 mm [m]
nz = 10_000                                             # 10 000 steps reach 10 m
z_mine = np.arange(nz + 1) * dz_step                    # heights 0, −0.001, …, −10 m
p_mine = np.empty(nz + 1)                               # pressure at every height [Pa]
p_mine[0] = ch01.P_ATM                                  # start: atmospheric pressure at the surface
for k in range(nz):                                     # Euler step of dp/dz = −ρ g (Eq. 1.8)
    p_mine[k + 1] = p_mine[k] - rho_fn(z_mine[k], p_mine[k]) * ch01.G0 * dz_step   # p_{k+1} = p_k − ρ g Δz
assert np.allclose(p_mine[::50], p_strat, rtol=1e-6)    # every 50th point matches the 201-point library profile
print(f"Euler at −10 m: {p_mine[-1]:.3f} Pa; library: {p_strat[-1]:.3f} Pa")   # hand march vs library
""")
nb.md("""
Euler's error here is about ½ρ′g|Δz| per metre of depth — roughly 0.02 Pa after 10 m, far inside a relative tolerance of
10⁻⁶.
""")
nb.figure("""
from scripts.ch01_drawings import draw_capillary                        # our capillary-tube drawing (geometry only)
fig, ax = plt.subplots(figsize=(7.2, 4.4))                              # one panel: pressure across, height up
ax.plot(p_w / 1e3, z, color=COLORS["teal"], label="uniform water (1000 kg/m³)")          # straight, slope −ρg
ax.plot(p_layers / 1e3, z, color=COLORS["orange"], label="1 m oil (800) over water")    # kink at z = −1 m
ax.plot(p_strat / 1e3, z, color=COLORS["blue"], ls="--", label="salty lake, ρ = 1000 − 0.5 z")   # slightly curved
ax.axhline(-1.0, color=COLORS["muted"], lw=0.8, ls=":")                 # the oil–water interface
ax.set_xlabel("absolute pressure p [kPa]")                              # bottom axis: absolute
ax.set_ylabel("height z [m] (z = 0 at the surface)")                    # z up
top = ax.secondary_xaxis("top", functions=(lambda P: P - ch01.P_ATM / 1e3, lambda G: G + ch01.P_ATM / 1e3))   # same axis in gauge
top.set_xlabel("gauge pressure [kPa]")                                  # top axis: gauge = absolute − 101.325 kPa
ax.set_title("Still water: pressure falls with height at slope −ρg")    # the message
ax.legend(loc="upper right", fontsize=8.5)                              # legend where the lines are not
ax_in = ax.inset_axes([0.02, 0.03, 0.36, 0.5])                           # small inset in the empty lower-left corner
draw_capillary(ax_in, R=1.0, alpha=np.radians(90), h=3.0)                # N12: water rises in a thin tube
ax_in.text(0.5, -0.02, "capillary rise (C22): E under the meniscus, F at the outside level", transform=ax_in.transAxes,
           ha="center", va="top", fontsize=7)                           # inset label below the drawing (σ arrow sits on top)
savefig(fig, "ch01", "c20_pressure_profiles")                           # keep a copy for review
plt.show()                                                              # display
""", see="""
Three nearly straight lines running from 101 kPa at the surface to about 198–200 kPa at 10 m depth: teal for uniform
water, orange with a kink at 1 m where oil meets water, dashed blue for the salty lake curving very slightly. The top axis
reads the same pressures as gauge values; the inset shows water standing higher inside a thin tube.
""", read="""
The slope of each line is −ρg, so the kink marks a density jump and the gentle curvature a density that grows with depth.
Gauge and absolute pressure differ only by a constant, so the physics (the slope) is the same on both axes.
""", change="""
Mercury (13 600 kg/m³) would make the line 13.6 times steeper — that is why a barometer column of mercury is only 76 cm
tall, while a water barometer would need about 10 m.
""")
nb.plotly("""
z_tank = np.linspace(0, -10, 201)                                     # heights in a 10 m tank [m]
z_faces = np.array([-1.45, -1.55])                                    # top and bottom faces of a 10 cm cube centred 1.5 m down
def tank(rho_top):                                                    # called once per slider value (upper-layer density)
    p_g = ch01.layered_pressure(z_tank, [2.0, 8.0], [rho_top, 1000.0]) - ch01.P_ATM   # gauge p(z) of 2 m oil over 8 m water
    p_f = ch01.layered_pressure(z_faces, [2.0, 8.0], [rho_top, 1000.0]) - ch01.P_ATM  # gauge pressure on the two faces
    return {"p(z) [kPa gauge]": (p_g / 1e3, z_tank), "cube faces": (p_f / 1e3, z_faces)}   # two traces: profile and face markers
rhos = np.linspace(600, 1000, 25)                                      # 25 slider positions [kg/m³]
fig = slider_figure(tank, "ρ upper layer", rhos, unit="kg/m³", xlabel="gauge pressure [kPa]", ylabel="height z [m]",   # one slider position per density
                    title="A lighter upper layer: gentler slope above the interface, less buoyancy on the cube",   # the message
                    modes={"cube faces": "markers"}, active=10)   # faces as dots; start mid-range
for i, r in enumerate(rhos):                                          # add a per-position title: plotly "update" steps take
    Fb = r * ch01.G0 * 1e-3                                           #   (data changes, layout changes); buoyancy ρ g V [N]
    st = fig.layout.sliders[0].steps[i]                               # the slider step for this density
    st.args = [st.args[0], {"title.text": f"upper layer {r:.0f} kg/m³: face-pressure gap × area = ρgV = {Fb:.2f} N"}]   # keep the visibility, add a title
fig.show()                                                            # drag the slider
""", explain="""
1. For each of 25 upper-layer densities, `layered_pressure` gives the gauge pressure profile and the pressures on the top
   and bottom faces of a 10 cm cube sitting in the upper layer.
2. The loop adds a title to every slider position: plotly's slider steps (method "update") accept a second dictionary of
   layout changes, here the title text with the buoyancy ρgV of the cube.
""")
see_read_change(
    'The gauge pressure profile of a 2 m oil layer over 8 m of water, with a kink at the interface, and two markers for the pressure on the top and bottom faces of a small cube in the upper layer.',
    'The slope of each straight piece is −ρg of that layer; the horizontal gap between the two markers is the pressure difference across the cube, and the title turns it into the buoyancy force.',
    'Slide the upper density toward 1000 kg/m³: the kink disappears (one fluid), and the face-pressure gap and the buoyancy grow to the values for water.',
)
nb.md("""
**What would change if…** the fluid were air over 1 km? Its density falls with height because air is compressible, so
(1.8) still holds but p(z) is no longer a straight line — that needs the equation of state (§1.9, C40). **Next:** a fluid
particle also has a temperature and an energy; §1.8 is its thermodynamics (C25).
""")

# =====================================================================================================================
# A.8 §1.8 Classical Thermodynamics — C25, C35, C36
# =====================================================================================================================
nb.section("1.8", "Classical Thermodynamics", intro="""
**What is this section about?** The energy bookkeeping of a fluid particle: what heat and work do to it (the first law),
which quantities describe its state, entropy and the Gibbs relations, and two material properties that the rest of the
book needs — the speed of sound and the thermal expansion coefficient. (The fluid particle itself, C24, was already met
in §1.4.)
""")
# ---------------------------------------------------------------------------------------------------- C25
core("C25", "The first law: heat in plus work done changes the internal energy",
     "Take a kilogram of air from one state to another by two different routes. What is the same at the end, and what depends on the route?")
nb.md("""
#### The problem in plain words
A bicycle pump gets warm; an engine turns heat into work cycle after cycle; air rising over a mountain cools although
nobody removes any heat. To follow the temperature of a fluid particle we need a bookkeeping rule for energy — Chapter 4
turns it into the energy equation.
""")
nb.md("""
#### The idea
Think of a bank account. The **balance** (the internal energy e) depends only on where you are now. Deposits arrive
through two channels — **heat** q and **work** w — and how much came through each channel depends on the route you took.

```
p ↑   1 ●
      │  ╲   route A: along the isotherm (heat flows in while the gas expands)
      │ │ ╲
      │ │  ╲________● 2
      │ └──────────┘      route B: cool at fixed volume, then expand at fixed pressure
      └────────────────────► v
same start, same end  ⇒  same Δe;   different routes  ⇒  different q and w
```
""")
note("C23", "System, equilibrium and relaxation", """
A **thermodynamic system** is a lump of matter we choose to follow whose molecules stay the same: energy may cross its
boundary as heat or as work, but matter may not. After a disturbance it **relaxes** back to **equilibrium**, where p, T and the other properties are well defined; the time
this takes is the **relaxation time** — in air about the time between molecular collisions, ≈ 10⁻¹⁰ s (C24). That is why a
fluid particle can be treated as always being in equilibrium.
""")
P("P32", "internal and kinetic energy per unit mass", """
Internal energy e [J/kg] is the energy of the molecules' random motion (and of their interactions); the particle's bulk
motion adds kinetic energy u²/2 [J/kg]. Thermodynamics here uses only e.
""", code="""
print(0.5 * 10.0**2)        # kinetic energy of a 10 m/s breeze → 50 J/kg
print(717.6 * 288.15)       # internal energy stored in air at 15 °C, e = C_v T (working model of the next primer) → ≈ 2.07e5 J/kg
""")
P("P33", "perfect-gas law as a working model", """
Until §1.9 derives it, we treat air as a perfect gas so that the examples have numbers: p = ρRT (equivalently pv = RT)
with R = 287 J kg⁻¹ K⁻¹, and e = C_vT with C_v = 718 J kg⁻¹ K⁻¹ (C40 and C41 show why).
""", code="""
v1 = ch01.R_AIR * 300.0 / 1e5       # specific volume of air at 300 K and 1 bar, v = RT/p [m³/kg]
print(v1)                           # → 0.861 m³/kg
""")
P("P34", "differentials", r"""
dx is a tiny change of a quantity that has a value in every state — an **exact** differential: adding up the pieces
gives x₂ − x₁ whatever the route. δq and δw are tiny amounts *transferred* along a route; they are not changes of
anything (**inexact** differentials), so their sums depend on the route.
""", code="""
steps_a = [50.0, 50.0, 50.0]                     # warm 300 → 450 K in three equal steps [K]
steps_b = [200.0, -100.0, 25.0, 25.0]            # or overshoot, cool, then creep up [K]
print(sum(steps_a), sum(steps_b))                # the sum of dT is 150 K both ways: T is a state function
""")
nb.md(r"""
#### The maths, step by step
1. Per unit mass, energy is conserved: the heat added δq plus the work done **on** the particle δw equals the rise of
   its internal energy,
   $$\delta q + \delta w = \Delta e \qquad (1.10) \quad [\mathrm{J/kg}].$$
2. Sign convention: work done *on* the gas is positive (compression), so an expanding gas has w < 0.
3. For a slow, frictionless (reversible) change the only work is the boundary moving: δw = −p dv (C27 below).
4. Along a route from state 1 to state 2 the work is $w = -\int_1^2 p\,dv$ — minus the area under the route on a p–v
   diagram.
5. Then the heat follows from the first law: q = Δe − w.
""")
note("N13", "Specific volume", """
v = 1/ρ [m³/kg] is the specific volume; it appears in (1.11)–(1.18) and in the energy equation of Ch. 4 §4.8.
""")
note("C26", "Path functions and state functions", """
**Heat and work are path functions** — energy in transit across the boundary. The internal energy is a **state function**
— a property with one value in each state, so Δe = e₂ − e₁ ignores the route. A **reversible** process is carried out so
slowly and without friction that the system is always in equilibrium; only then is the route a line on a p–v diagram.
Number (code below): routes A, B and C between the same states all have Δe = 0, but w = −59.7, −43.1 and −86.1 kJ/kg.
""")
note("C27", "The reversible first law", r"""
For a reversible change with boundary work only, the work done on unit mass is −p dv (a piston face of area A pushed out
by dx does work pA dx on the surroundings, and A dx is the volume increase), so (stated) the first law reads as below.
Friction work is excluded. Doubling the volume of air at a constant 0.5 bar: w = −p Δv = −43.1 kJ/kg.
""", equation=r"de = dq - p\,dv", ref="1.11")
note("C28", "Equations of state", r"""
An **equation of state** ties state functions together; for a simple one-component fluid **two** properties fix the state,
for example p = p(v, T) (thermal) or e = e(p, T) (caloric). Seawater needs a third, salinity (§1.10). Air at
v = 0.861 m³/kg and 300 K has p = 1.00×10⁵ Pa (code below).
""", equation=r"p = p(v, T), \qquad e = e(p, T)", ref="1.12")
P("P35", "line integral along a path", r"""
$\int p\,dv$ along a route adds p × (tiny volume change) step by step along that route; on a p–v diagram it is the area
under the route (negative if v decreases along it). A different route between the same end points encloses a different
area.
""", code="""
print(5e4 * (1.722 - 0.861))      # a horizontal route at p = 5e4 Pa from v = 0.861 to 1.722 m³/kg → area 43050 J/kg
""")
P("P36", "natural logarithm and exponential", r"""
ln x undoes eˣ; $\int dv/v = \ln(v_2/v_1)$; ln 2 = 0.693. e⁻¹ ≈ 0.368 — a quantity that shrinks by that factor has decayed
by one **e-folding**.
""", code="""
print(np.log(2.0), np.exp(-1.0))   # ln 2 = 0.6931, e^(−1) = 0.3679
""")
nb.worked_example("doubling the volume of 1 kg of air at 300 K, three ways", """
State 1: T = 300 K, p = 10⁵ Pa → v₁ = RT/p = 287.06 × 300/10⁵ = 0.861 m³/kg. State 2: v₂ = 1.722 m³/kg, T = 300 K.
Δe = C_v(300 − 300) = 0 on **every** route.

- **A (isothermal):** p = RT/v, so w = −RT ln(v₂/v₁) = −287.06 × 300 × 0.693 = −59.7 kJ/kg and q = Δe − w = +59.7 kJ/kg.
- **B (cool at fixed v to 150 K, then expand at fixed p):** leg 1: w = 0, q = C_vΔT = 717.6 × (−150) = −107.6 kJ/kg.
  Leg 2: p = R·150/v₁ = 5.0×10⁴ Pa, w = −p(v₂ − v₁) = −43.1 kJ/kg, Δe = +107.6, q = +150.7.
  Totals: w = −43.1, q = +43.1 kJ/kg.
- **C (expand at fixed p to 600 K, then cool at fixed v):** w = −p₁(v₂ − v₁) = −86.1 kJ/kg, q = +86.1 kJ/kg.
""")
P("P37", "trapezoid rule", r"""
Approximates $\int f\,dx$ by adding up trapezoids between neighbouring samples; the error shrinks like Δx².
`np.trapezoid(f, x)` returns the total and `scipy.integrate.cumulative_trapezoid` the running sum.
""", code="""
xs = np.linspace(1, 2, 101)                          # 101 samples between 1 and 2
print(np.trapezoid(1 / xs, xs), np.log(2.0))         # ∫ dx/x from 1 to 2: 0.693153 vs ln 2 = 0.693147
""")
nb.code("""
routes = {"isothermal": "A: along the 300 K isotherm",   # route A
          "isochoric-isobaric": "B: cool at fixed v, then expand at fixed p",   # route B
          "isobaric-isochoric": "C: expand at fixed p, then cool at fixed v"}      # the three routes of the example
tot = {}                                                                          # exact totals per route
for kind, label in routes.items():                                                # same start and end state every time
    tot[kind] = ch01.path_heat_work_totals(kind, v1, 300.0, 2 * v1, 300.0)        # leg-by-leg q, w, Δe [J/kg]
    print(f"{label:45s} q = {tot[kind]['q']/1e3:7.2f}  w = {tot[kind]['w']/1e3:7.2f}  Δe = {tot[kind]['de']/1e3:5.2f} kJ/kg")   # q, w, Δe of this route [kJ/kg]
path_B = ch01.process_path("isochoric-isobaric", (v1, 300.0), (2 * v1, 300.0))    # route B sampled: arrays v, T, p
run_B = ch01.process_heat_work(path_B["v"], path_B["T"])                          # running w = −∫p dv, Δe, q along it
print(f"route B by trapezoid sums: w = {run_B['w'][-1]/1e3:.2f} kJ/kg")          # same as the exact total
print("C28: state (p, ρ, T) of air at v = 0.861 m³/kg, 300 K:", ch01.perfect_gas_state(rho=1 / v1, T=300.0))   # two properties fix the third
""", explain="""
1. `path_heat_work_totals` computes q, w and Δe leg by leg with exact formulas for the three routes between the same two
   states: Δe is 0 each time, q and w are not.
2. `process_path` samples route B as arrays of v, T and p, including its corner.
3. `process_heat_work` accumulates w = −∫p dv with the trapezoid rule (P37), Δe = C_vΔT and q = Δe − w along the samples.
4. `perfect_gas_state` fills in the missing one of (p, ρ, T): two properties fix the state (C28).
""")
nb.md("""
#### From scratch: trapezoids along each route
""")
nb.check_agree("""
for kind in routes:                                                        # the three routes again
    path = ch01.process_path(kind, (v1, 300.0), (2 * v1, 300.0))           # sampled route
    v_s, T_s = path["v"], path["T"]                                        # specific volume [m³/kg] and temperature [K]
    p_s = ch01.R_AIR * T_s / v_s                                           # pressure along the route, p = RT/v [Pa]
    w = 0.0                                                                # work done ON the gas [J/kg]
    for k in range(len(v_s) - 1):                                          # one trapezoid per pair of neighbouring samples
        w -= 0.5 * (p_s[k] + p_s[k + 1]) * (v_s[k + 1] - v_s[k])           # δw = −p dv (Eq. 1.11)
    de = ch01.CV_AIR * (T_s[-1] - T_s[0])                                  # Δe = C_v ΔT (perfect gas)
    q = de - w                                                             # first law (1.10): q = Δe − w
    assert np.allclose([q, w, de], [tot[kind]["q"], tot[kind]["w"], tot[kind]["de"]], rtol=1e-4, atol=1e-6)   # = exact totals
    assert np.isclose(w, ch01.process_heat_work(v_s, T_s)["w"][-1])        # = the library's running sum
    print(f"{kind:20s} mine: q = {q/1e3:.3f}, w = {w/1e3:.3f} kJ/kg")   # our totals
""")
nb.figure("""
fig, ax = plt.subplots(figsize=(7.2, 4.2))                                        # one p–v diagram
vv = np.linspace(0.7, 1.9, 200)                                                   # specific volumes for the isotherms [m³/kg]
for T_iso in (150.0, 300.0, 600.0):                                               # faint isotherms p = RT/v
    ax.plot(vv, ch01.R_AIR * T_iso / vv / 1e3, color=COLORS["grid"], lw=1.2)    # one isotherm [kPa]
    ax.text(1.88, ch01.R_AIR * T_iso / 1.88 / 1e3, f"{T_iso:.0f} K", color=COLORS["muted"], fontsize=8, va="bottom")  # its label
colors = {"isothermal": COLORS["blue"], "isochoric-isobaric": COLORS["orange"], "isobaric-isochoric": COLORS["rose"]}  # route colours
for kind, label in routes.items():                                                # draw the three routes
    path = ch01.process_path(kind, (v1, 300.0), (2 * v1, 300.0))                  # sampled route (v, T, p)
    ax.plot(path["v"], path["p"] / 1e3, color=colors[kind], lw=2.5, label=f"{label[:1]}: w = {tot[kind]['w']/1e3:.1f} kJ/kg")  # the route
    mid = len(path["v"]) * 3 // 4                                                 # put a direction arrow 3/4 along the route
    ax.annotate("", xy=(path["v"][mid + 1], path["p"][mid + 1] / 1e3), xytext=(path["v"][mid - 3], path["p"][mid - 3] / 1e3),
                arrowprops=dict(arrowstyle="-|>", color=colors[kind], lw=2))     # arrow head showing the direction
leg2 = path_B["corner"]                                                           # index where route B turns the corner
ax.fill_between(path_B["v"][leg2:], 0, path_B["p"][leg2:] / 1e3, color=COLORS["blue"], alpha=0.15,   # shade between p = 0 and the route
                label="area under B = work done by the gas")                      # fill_between shades between two curves
ax.plot([v1, 2 * v1], [1e2, 50.0], "o", color=COLORS["ink"], ms=7)                # state 1 and state 2
ax.text(v1 - 0.03, 101, "1", ha="right", fontsize=11); ax.text(2 * v1 + 0.02, 52, "2", fontsize=11)   # label the two states
ax.set_xlim(0.7, 1.95); ax.set_ylim(0, 220)                                       # fixed axes
ax.set_xlabel("specific volume v [m³/kg]"); ax.set_ylabel("pressure p [kPa]")    # axis labels with units
ax.set_title("Same two states, three routes, three different amounts of work")   # the message
ax.legend(loc="upper right", fontsize=8.5)                                        # legend with the work of each route
savefig(fig, "ch01", "c25_pv_routes")                                             # keep a copy for review
plt.show()                                                                        # display
""", see="""
Three routes from state 1 (0.861 m³/kg, 100 kPa) to state 2 (1.722 m³/kg, 50 kPa): the blue isotherm, the orange route
that drops first and then runs across, and the rose route that runs across at 100 kPa and then drops. The shaded area
lies under route B's horizontal leg.
""", read="""
The area under a route is the work done *by* the gas (−w). Because Δe = 0 between these two states, the heat taken in
equals that same area — so the three routes need 43, 60 and 86 kJ/kg of heat to arrive at the same place.
""", change="""
Run route C backwards (2 → 1) and every w and q changes sign. Going out along C and back along A encloses an area: the
net work of a heat engine cycle, paid for by the difference in heat.
""")
nb.animation("""
nfr = 20 if not FAST else 12                                                   # frames per route
legs = []                                                                      # (route label, v, T, p, running q, w, Δe)
for kind in ("isothermal", "isochoric-isobaric"):                              # route A, then route B
    path = ch01.process_path(kind, (v1, 300.0), (2 * v1, 300.0))               # sampled route
    run = ch01.process_heat_work(path["v"], path["T"])                         # running q, w, Δe along it [J/kg]
    sel = np.linspace(0, len(path["v"]) - 1, nfr).astype(int)                  # the samples shown as frames
    legs += [(routes[kind][:1], path["v"][j], path["T"][j], path["p"][j], run["q"][j], run["w"][j], run["de"][j],   # one tuple per frame (continued)
              run["q"][min(j + 1, len(path["v"]) - 1)] - run["q"][j]) for j in sel]   # last item: sign of the next δq
fig, (axc, axd) = plt.subplots(1, 2, figsize=(8.4, 3.3), gridspec_kw=dict(width_ratios=[1, 1.2]))   # left: cylinder, right: p–v plane
axc.set_xlim(0, 2.3); axc.set_ylim(-0.6, 1.3); axc.axis("off")                # left: the cylinder drawing
axc.plot([0, 2.2, 2.2, 0, 0], [0, 0, 1, 1, 0], color=COLORS["ink"], lw=1.5)   # cylinder walls
gas = plt.Rectangle((0, 0), 1, 1, color=plt.cm.coolwarm(0.5))                  # the gas; its colour follows T
axc.add_patch(gas)   # put the gas on the drawing
piston = plt.Rectangle((1, -0.05), 0.08, 1.1, color=COLORS["muted"])          # the piston face
axc.add_patch(piston)   # put the piston on the drawing
heat_txt = axc.text(0.6, -0.35, "", ha="center", fontsize=10)                  # heater / cooler label
count_txt = axc.text(0.0, 1.15, "", fontsize=8.5, family="monospace")          # running counters
for T_iso in (150.0, 300.0, 600.0):                                            # right: faint isotherms on the p–v plane
    axd.plot(vv, ch01.R_AIR * T_iso / vv / 1e3, color=COLORS["grid"], lw=1)   # one isotherm [kPa]
(trace,) = axd.plot([], [], color=COLORS["accent"], lw=2.5)                   # the route so far
(dot,) = axd.plot([], [], "o", color=COLORS["orange"], ms=8)                   # the current state
axd.set_xlim(0.7, 1.9); axd.set_ylim(0, 220)   # fixed axes
axd.set_xlabel("v [m³/kg]"); axd.set_ylabel("p [kPa]")   # axis labels with units
def update(i):                                                                 # frame i: one state along A, then along B
    lab, v, T, p, q, w, de, dq_next = legs[i]   # unpack this frame's state (P14)
    start = 0 if i < nfr else nfr                                              # the current route's first frame
    x_p = 2.0 * v / (2 * v1)                                                   # piston position ∝ volume
    gas.set_width(x_p); piston.set_x(x_p)                                      # move the piston
    gas.set_color(plt.cm.coolwarm((T - 150.0) / 450.0))                        # 150 K blue … 600 K red
    heat_txt.set_text("heat flows IN" if dq_next > 1e-9 else ("heat flows OUT" if dq_next < -1e-9 else ""))   # sign of the next δq
    heat_txt.set_color(COLORS["orange"] if dq_next > 0 else COLORS["blue"])   # orange in, blue out
    count_txt.set_text(f"route {lab}  T = {T:5.1f} K\\nq = {q/1e3:6.1f}  w = {w/1e3:6.1f}  Δe = {de/1e3:6.1f} kJ/kg")   # running totals
    trace.set_data([s[1] for s in legs[start:i + 1]], [s[3] / 1e3 for s in legs[start:i + 1]])   # route so far on the p–v plane
    dot.set_data([v], [p / 1e3])   # current state
    return gas, piston, heat_txt, count_txt, trace, dot   # the artists that changed
show_animation(animate(update, frames=len(legs), fig=fig, interval=150), player="frames", dpi=64)   # step through; dpi 64 keeps it small
""", explain="""
1. For routes A and B, `process_path` samples the states and `process_heat_work` gives the running q, w and Δe.
2. Each frame draws the piston at a position proportional to v, colours the gas by its temperature (blue cold, red hot)
   and says whether heat is flowing in or out next.
3. The right panel traces the route on the p–v plane; the counters show the running totals.
""")
nb.md("""
**What you see.** Both routes end with the piston in the same place and the gas the same colour; on route B the gas first
turns blue (cooling at fixed volume) and then warms back while it expands.

**How to read it.** The counters for q and w end at different values on the two routes, while Δe comes back to zero on
both — the state function ignores the route, the path functions do not.

**What would change if…** the route went through 600 K instead (route C): the expansion would happen at twice the
pressure, so it would need twice the work and twice the heat of route B.
""")
nb.md("""
**What would change if…** we divided each little δq by the temperature at which it arrives before adding? The next block
shows that this sum is the same on every reversible route — a new state function, entropy (C35).
""")
# ---------------------------------------------------------------------------------------------------- C35
core("C35", "Entropy and the Gibbs relations: a heat-like quantity that is a state function",
     "Heat depends on the route. Can we divide it by something so that the result no longer does?")
nb.md("""
#### The problem in plain words
Why does a hot cup cool down but never warm itself up? Why can the potential temperature of rising air stay constant
(§1.10)? Both need a property that tracks heat but belongs to the state.
""")
nb.md("""
#### The idea
The three routes of C25 once more:

| route | q [kJ/kg] | ∫ dq/T [J kg⁻¹ K⁻¹] |
|---|---|---|
| A | 59.7 | 199.0 |
| B | 43.1 | 199.0 |
| C | 86.1 | 199.0 |

**Heat depends on the route; heat divided by temperature, summed along a reversible route, does not.** That sum is the
entropy change.
""")
P("P39", "partial derivative with a variable held fixed", r"""
In thermodynamics every property depends on two others, so a partial derivative must say which one is held fixed:
$(\partial h/\partial T)_p$ is the rate of change of h with T **at constant p**. Changing what is held fixed changes the
answer — $(\partial e/\partial T)_v$ and $(\partial e/\partial T)_p$ are different numbers for most fluids.
""", code="""
h_fun = lambda T, p: 1004.7 * T                                    # enthalpy of a perfect gas, h = C_p T [J/kg]
print(ch01.partial_derivative(h_fun, "T", {"T": 300.0, "p": 1e5})) # (∂h/∂T) at fixed p → 1004.7 J/(kg K)
""")
note("C29", "Enthalpy", r"""
h = e + pv [J/kg]: internal energy plus the "flow work" pv needed to push the particle into its place. It is the natural
energy for constant-pressure processes and for flows. Air at 300 K: h − e = pv = RT = 86.1 kJ/kg.
""", equation=r"h \equiv e + pv", ref="1.13")
note("C30 · C31", "Specific heats", r"""
How fast h and e grow with temperature — C_p at constant pressure, C_v at constant volume [J kg⁻¹ K⁻¹]. They are
properties (defined from other properties, not from heat). Air: C_p = 1004.7, C_v = 717.6 J kg⁻¹ K⁻¹.
""", equation=r"C_p \equiv \left(\frac{\partial h}{\partial T}\right)_p, \qquad C_v \equiv \left(\frac{\partial e}{\partial T}\right)_v",
     ref="1.14, 1.15")
note("C32", "Heat per degree", """
Only along a reversible constant-volume (constant-pressure) route with nothing but p dv work is the heat per kelvin equal
to C_v (C_p) — stirring heats a gas with no heat at all; this is why h is written C_pT in Ch. 4 §4.8 and in stagnation
enthalpy, Ch. 15 §15.4.
""")
note("C33", "Second law (i): entropy exists", r"""
There is a state function s [J kg⁻¹ K⁻¹] whose change is the sum of dq/T along **any reversible** route between the two
states. Routes A, B and C all give 199.0 J kg⁻¹ K⁻¹ = R ln 2.
""", equation=r"s_2 - s_1 = \int_1^2 \frac{dq_{\rm rev}}{T}", ref="1.16")
note("N15", "The differential form", "For a reversible change the heat received is T times the entropy change:",
     equation=r"T\,ds = dq", ref="1.17")
P("P38", "product rule for differentials", r"""
d(pv) = p dv + v dp — both factors may change, and each change contributes. Forgetting one of the two terms is the classic
slip.
""", code="""
p_a, v_a, p_b, v_b = 1e5, 0.861, 1.001e5, 0.8615           # two nearby states (pressure [Pa], volume [m³/kg])
exact = p_b * v_b - p_a * v_a                              # the true change of pv [J/kg]
rule = p_a * (v_b - v_a) + v_a * (p_b - p_a)               # p dv + v dp
print(exact, rule, exact - rule)                           # they differ only by Δp·Δv = 0.05 J/kg — second order
""")
D("D10", "The Gibbs relations", ref="1.18",
  goal="""Eliminate heat — which depends on the route — from the reversible first law, leaving relations between state
functions only. Those relations then hold for *any* process.""",
  assumptions="""A reversible route is used as a device only (steps 1–2) · simple compressible substance: only p dv work
(step 1) · states 1 and 2 are equilibrium states (step 7).""",
  start=(r"T\,ds = dq \quad (1.17) \qquad\text{and}\qquad de = dq - p\,dv \quad (1.11)",
         "for a reversible change the heat received is T ds, and the first law has boundary work −p dv."),
  plan=["Solve the first law for dq and substitute it into T ds = dq.",
        "Rewrite with enthalpy using the product rule.",
        "Explain why a route-independent relation holds for every process."],
  uses=["(1.17) (N15)", "(1.11) (C27, in C25)", "h = e + pv (C29)", "differentials (P34)",
        "product rule for differentials (P38)"],
  steps=[
      ("Solve the reversible first law for dq", r"dq = de + p\,dv",
       "Add p dv to both sides of (1.11). We isolate dq because it is the route-dependent quantity we want to replace.",
       "Heat in = rise of internal energy + work done by the gas."),
      ("Substitute into T ds = dq", r"T\,ds = de + p\,dv",
       "(1.17) says the same dq equals T ds on a reversible route; equal things can replace each other in an equation.",
       "The first Gibbs relation — no heat left in it."),
      ("Differentiate the definition of enthalpy", r"dh = de + d(pv)",
       "h = e + pv (1.13); the differential of a sum is the sum of the differentials. We want a second form written with dh.",
       "A change in h is a change in e plus a change in pv."),
      ("Apply the product rule", r"dh = de + p\,dv + v\,dp",
       "d(pv) = p dv + v dp (P38): both p and v may change, and each change contributes its own term.",
       "Both pressure and volume changes feed pv."),
      ("Isolate de + p dv", r"de + p\,dv = dh - v\,dp",
       "Subtract v dp from both sides of step 4; the left side is now exactly the right side of step 2.",
       "The same combination, written with enthalpy."),
      ("Substitute into step 2", r"T\,ds = dh - v\,dp",
       "Replace de + p dv in step 2 by its equal from step 5; substituting equals for equals keeps the equation true.",
       "The second Gibbs relation."),
      ("Drop the reversible restriction", r"T\,ds = de + p\,dv = dh - v\,dp\quad\text{(any process)}",
       "T, s, e, h, p and v are all state functions, so a relation between their changes for neighbouring equilibrium "
       "states cannot know which route took the fluid there; any irreversible change can be replaced by a reversible route "
       "between the same states. Only the *interpretation* T ds = dq is lost.",
       "The relations hold for friction and stirring too — but then T ds is no longer the heat."),
  ],
  result=(r"T\,ds = de + p\,dv, \qquad T\,ds = dh - v\,dp",
          "entropy changes are fixed by changes of energy (or enthalpy), volume and pressure, for any process."),
  interpret="""Entropy is computable from measurable state properties alone. For stirring at constant volume (dv = 0) the
first form gives T ds = de: the stirring work shows up as entropy although no heat entered — the Clausius–Duhem inequality
(C34 below). The isentropic law (D14) and the adiabatic lapse rate (D19) both start from these relations. In the explainer
`heat_work_paths` (C45 block), switch on the stirring mode: the Gibbs relation still gives Δs = 23.5 J/(kg K) with δq = 0.""",
  check="""Units: J/kg on every term ✓. Perfect gas heated from 300 to 600 K at constant p: form 1 gives C_v ln 2 + R ln 2
= 696.4 J/(kg K), form 2 gives C_p ln 2 − R ln 1 = 696.4 ✓ (worked example and code below). The sympy cell below shows
T ds − de − p dv = 0 and T ds − dh + v dp = 0 symbolically for s = C_v ln T + R ln v.""",
  traps="concluding that the relations hold only for reversible processes; dropping v dp in the product rule.")
note("C34", "Second law (ii): the Clausius–Duhem inequality", r"""
For *any* process the entropy change is at least the actual heat received divided by the temperature at which it arrives,
with equality only when the process is reversible. (We write the actual heat δq; with dq_rev on the right it would simply
equal Δs by (1.16).) Free expansion of air into vacuum to twice the volume: q = w = Δe = 0, yet Δs = R ln 2 = 199 J kg⁻¹ K⁻¹
> 0. Stirring from 300 to 310 K at fixed v: δq = 0, Δs = C_v ln(310/300) = 23.5 J kg⁻¹ K⁻¹ > 0.
""", equation=r"s_2 - s_1 \ge \int_1^2 \frac{\delta q}{T}")
note("N14", "Second law (iii)", """
μ and k must be positive — otherwise momentum or heat would un-mix by itself — so viscous dissipation is never negative
(Ch. 4 §4.8).
""")
nb.worked_example("heating air from 300 K to 600 K at 1 bar", """
Constant p, so v doubles too.
- First Gibbs form, integrated for a perfect gas: Δs = C_v ln(T₂/T₁) + R ln(v₂/v₁) = 717.6 × 0.693 + 287.1 × 0.693 = 497.4 + 199.0
  = 696.4 J kg⁻¹ K⁻¹.
- Second form: Δs = C_p ln(T₂/T₁) − R ln(p₂/p₁) = 1004.7 × 0.693 − 0 = 696.4 J kg⁻¹ K⁻¹.

Same number, two routes through the algebra.
""")
nb.code("""
print(ch01.perfect_gas_entropy_change(300.0, v1, 600.0, 2 * v1))       # first Gibbs form integrated: Δs [J/(kg K)] → 696.4
print(ch01.perfect_gas_entropy_change_p(300.0, 1e5, 600.0, 1e5))       # second form, with p: → 696.4
path_p = ch01.process_path("isobaric", (v1, 300.0), (2 * v1, 600.0))   # the same heating as a sampled isobaric route
run_p = ch01.process_heat_work(path_p["v"], path_p["T"])               # running q along it [J/kg]
print(ch01.entropy_change_reversible(run_p["q"], path_p["T"]))         # ∫ dq/T along the route (1.16) → ≈ 696.4
print(ch01.irreversible_process("free_expansion", 300.0, v1, v2=2 * v1))   # C34: q = w = Δe = 0 but Δs = 199
print(ch01.irreversible_process("stirring", 300.0, v1, T2=310.0))          # C34: δq = 0 but Δs = 23.5
""", explain="""
1. `perfect_gas_entropy_change` and `perfect_gas_entropy_change_p` are the two Gibbs forms integrated for a perfect gas.
2. `entropy_change_reversible` sums dq/T along the sampled isobaric route — the definition (1.16) — and agrees.
3. `irreversible_process` returns the totals for free expansion and for stirring: in both, `int_dq_over_T` is 0 while
   `ds` is positive — the Clausius–Duhem inequality with actual heat.
""")
P("P40", "sympy", """
sympy does algebra with symbols: `sp.symbols` makes symbols, `sp.Function("T")(t)` an unknown function, `sp.diff`
differentiates, `sp.simplify` tidies up. An expression that simplifies to 0 is an identity — true for every value of the
symbols.
""", code="""
import sympy as sp                                  # symbolic mathematics
x = sp.symbols("x")                                 # a symbol, not a number
print(sp.simplify(sp.diff(x**3, x) - 3 * x**2))     # d(x³)/dx − 3x² → 0: an identity
""")
nb.code("""
t = sp.symbols("t")                                        # a parameter that moves along an arbitrary route
cv_s, R_s = sp.symbols("c_v R", positive=True)             # constant specific heat and gas constant
T_t = sp.Function("T")(t)                                  # temperature along the route: any function of t
v_t = sp.Function("v")(t)                                  # specific volume along the route: any function of t
s_t = cv_s * sp.log(T_t) + R_s * sp.log(v_t)               # perfect-gas entropy (up to a constant)
e_t = cv_s * T_t                                           # internal energy e = C_v T
p_t = R_s * T_t / v_t                                      # pressure from p v = R T
h_t = e_t + p_t * v_t                                      # enthalpy (1.13)
d = lambda f: sp.diff(f, t)                                # "d" along the route = derivative with respect to t
print(sp.simplify(T_t * d(s_t) - d(e_t) - p_t * d(v_t)))   # first Gibbs form residual → 0
print(sp.simplify(T_t * d(s_t) - d(h_t) + v_t * d(p_t)))   # second Gibbs form residual → 0
""", explain="""
1. The route is left completely general: T(t) and v(t) are arbitrary functions.
2. With the perfect-gas s, e, p and h, both residuals T ds − de − p dv and T ds − dh + v dp simplify to zero — the Gibbs
   relations hold along every route, not only reversible ones we picked.
""")
nb.figure("""
fig, ax = plt.subplots(figsize=(7.2, 4.2))                                        # one T–s diagram
s_rel = lambda T, v: ch01.perfect_gas_entropy_change(300.0, v1, T, v)             # entropy relative to state 1 [J/(kg K)]
ss = np.linspace(-450, 700, 200)                                                  # entropy range for the guide lines
ax.plot(ss, 300.0 * np.exp(ss / ch01.CV_AIR), color=COLORS["grid"], lw=1.2)        # isochore through state 1: T ∝ e^{s/C_v}
ax.plot(ss, 300.0 * np.exp(ss / ch01.CP_AIR), color=COLORS["grid"], lw=1.2, ls="--")  # isobar through state 1: T ∝ e^{s/C_p}
for kind, label in routes.items():                                                # the three routes of C25
    path = ch01.process_path(kind, (v1, 300.0), (2 * v1, 300.0))                  # sampled route (v, T)
    ax.plot(s_rel(path["T"], path["v"]), path["T"], color=colors[kind], lw=2.5, label=f"{label[:1]}: q = {tot[kind]['q']/1e3:.1f} kJ/kg")  # T against s
sB, TB = s_rel(path_B["T"], path_B["v"]), path_B["T"]                             # route B in the T–s plane
ax.fill_between(sB[leg2:], 0, TB[leg2:], color=COLORS["orange"], alpha=0.18, label="heat in on B's second leg")     # ∫T ds > 0
ax.fill_between(sB[:leg2 + 1], 0, TB[:leg2 + 1], color=COLORS["blue"], alpha=0.12, label="heat out on B's first leg")  # ∫T ds < 0
ax.plot([0, s_rel(300.0, 2 * v1)], [300, 300], "o", color=COLORS["ink"], ms=7)    # states 1 and 2
ax.set_xlim(-450, 700); ax.set_ylim(0, 650)                                       # fixed axes, T from 0 so areas are honest
ax.set_xlabel("entropy change s − s₁ [J/(kg K)]"); ax.set_ylabel("temperature T [K]")   # axis labels with units
ax.set_title("On a T–s diagram heat is the area under the route; every route ends at the same s")   # the message
ax.legend(loc="upper left", fontsize=8.5)                                         # legend
plt.show()                                                                        # display
""", see="""
The three routes of C25 redrawn with temperature against entropy: all start at s − s₁ = 0 and all end at the same point,
199 J/(kg K) to the right. The faint solid line is the isochore through state 1, the faint dashed line the isobar.
""", read="""
On a T–s diagram the heat of a reversible route is the area under it (∫T ds). Route B gives heat away on its first leg
(blue area) and takes more in on its second (orange area); the net is its 43 kJ/kg. Different areas mean different q, yet
all routes share one Δs.
""", change="""
A reversible adiabatic route has ds = 0, so it would be a vertical line on this diagram — the isentrope of C45.
""")
nb.md("""
**What would change if…** the process were not reversible (friction)? The Gibbs relations would still give Δs, because
they only link state functions; only T ds = dq fails. **Next:** how stiff a fluid is when squeezed at constant entropy sets
the speed of sound (C36).
""")
# ---------------------------------------------------------------------------------------------------- C36
core("C36", "The speed of sound: stiffness at constant entropy",
     "Why does a sound pulse cross 340 m of air in one second, but about 1480 m of water?")
nb.md("""
#### The problem in plain words
Thunder arrives after the lightning; sonar measures ocean depth; weather models filter sound waves out. More important for
this book: comparing a flow's speed with the speed of sound decides whether density changes can be ignored — the
incompressible approximation of Ch. 4 §4.9 — and the speed of sound is where Ch. 15 starts.
""")
nb.md("""
#### The idea
A pressure push squeezes the next layer of fluid, which pushes the next… The signal travels faster the **stiffer** the
fluid (more pressure change per density change) and slower the more **inertia** it has. The squeeze is far too quick for
heat to flow, so the stiffness is measured at constant entropy.

| squeeze by Δp = 1 kPa | fractional density change Δρ/ρ |
|---|---|
| air at 1 atm | ≈ 0.7 % |
| water | ≈ 4.5×10⁻⁷ |
""")
nb.md(r"""
#### The maths, step by step
1. **Which speed can there be?** The only properties that matter are the stiffness ∂p/∂ρ [Pa per kg/m³ = m²/s²] and the
   density. A stiffness already has the units of a speed squared, so the signal speed must be a number times
   √(∂p/∂ρ) (dimensional reasoning, previewing C64).
2. **Why the number is 1 — a piston push (plausibility, not a proof).** Push a piston of area A into still fluid at a small
   speed u. After a time t a pressure front has travelled ct and all the fluid behind it moves at u.
   - *Mass:* the fluid that filled the length ct now fills ct − ut, so its density has risen by Δρ ≈ ρu/c.
   - *Momentum:* the extra pressure Δp acting on A for the time t gave the mass ρAct the speed u, so Δp·A·t = ρAct·u,
     i.e. Δp = ρcu.
   - *Stiffness:* Δp = (∂p/∂ρ)Δρ = (∂p/∂ρ)ρu/c. Equating the two expressions for Δp gives c² = ∂p/∂ρ.
   The squeeze is too fast for heat to flow, so the stiffness is measured with entropy held fixed (P39). The full proof is in
   Ch. 15 §15.2:
   $$c^2 = \left(\frac{\partial p}{\partial \rho}\right)_s \qquad (1.19) \quad [\mathrm{m^2/s^2}]$$
3. Define the bulk modulus $K \equiv \rho\,(\partial p/\partial\rho)_s$ [Pa]: how many pascals per *fractional* density
   change. Then c² = K/ρ.
4. For air squeezed isentropically, $p = p_0(\rho/\rho_0)^\gamma$ (derived in C45), so $\partial p/\partial\rho =
   \gamma p_0\rho^{\gamma-1}/\rho_0^\gamma = \gamma p/\rho$ and K = γp.
5. For water K ≈ 2.2 GPa (the modified Tait model, `tait_pressure`).
6. **Incompressible limit:** if ∂ρ/∂p → 0 then ∂p/∂ρ → ∞ and c → ∞ — a pressure change is felt everywhere at once.
""")
nb.worked_example("air and water", """
- Air at 101 325 Pa, 1.225 kg/m³, γ = 1.4: c² = γp/ρ = 1.4 × 101 325/1.225 = 115 800 m²/s² → c = 340.3 m/s.
- Water: c = √(K/ρ) = √(2.2×10⁹/1000) = 1483 m/s.
- Compressing water by 1 % needs Δp = 0.01K = 22 MPa — the pressure about 2.2 km down in the ocean.
""")
nb.code("""
P_air = lambda rho, s: ch01.isentropic_pressure(rho, ch01.P_ATM, 1.225)   # air along its isentrope p = p0 (ρ/ρ0)^γ [Pa]
P_water = lambda rho, s: ch01.tait_pressure(rho)                         # water: modified Tait equation of state [Pa]
c_air = ch01.sound_speed_from_eos(P_air, 1.225)                          # Eq. (1.19): sqrt of dp/dρ at fixed s [m/s]
c_water = ch01.sound_speed_from_eos(P_water, 1000.0)                     # the same for water [m/s]
print(f"air {c_air:.2f} m/s, water {c_water:.1f} m/s, perfect-gas formula {ch01.perfect_gas_sound_speed(288.15):.2f} m/s")   # three speeds of sound
""", explain="""
1. Each equation of state is passed as a function p(ρ, s) (primer P29); the entropy argument is not used because both
   functions already describe one isentrope.
2. `sound_speed_from_eos` differentiates p with respect to ρ at fixed s and takes the square root — (1.19).
3. `perfect_gas_sound_speed` is the closed form √(γRT) that C47 derives; it agrees with the air value.
""")
nb.md("#### From scratch: a central difference of the equation of state")
nb.check_agree("""
c_mine = []                                                        # our speeds of sound [m/s]
for P, rho0 in ((P_air, 1.225), (P_water, 1000.0)):                # the two fluids
    drho = 1e-6 * rho0                                             # a tiny density step [kg/m³]
    slope = (P(rho0 + drho, 0.0) - P(rho0 - drho, 0.0)) / (2 * drho)   # central difference (P21) of p(ρ) = (∂p/∂ρ)_s
    c_mine.append(np.sqrt(slope))                                  # c = sqrt(slope), Eq. (1.19)
assert np.allclose(c_mine, [c_air, c_water], rtol=1e-6)            # same as the library
assert np.isclose(c_mine[0], ch01.perfect_gas_sound_speed(288.15), rtol=1e-4)   # and as √(γRT) for air at 15 °C
print(f"mine: air {c_mine[0]:.2f} m/s, water {c_mine[1]:.1f} m/s")   # ours
""")
nb.figure("""
K = np.logspace(4, 11, 200)                                               # bulk moduli from 10 kPa to 100 GPa [Pa]
fig, ax = plt.subplots(figsize=(7, 4))                                    # one log–log panel
for rho_c, col, lab in ((1.2, COLORS["teal"], "ρ = 1.2 kg/m³ (air-like)"), (1000.0, COLORS["blue"], "ρ = 1000 kg/m³ (water-like)")):  # two densities
    ax.loglog(K, np.sqrt(K / rho_c), color=col, label=lab)               # c = sqrt(K/ρ): slope ½ on log–log axes
K_air = ch01.GAMMA_AIR * ch01.P_ATM                                        # air's isentropic bulk modulus γp [Pa]
ax.loglog([K_air], [c_air], "o", color=COLORS["orange"], ms=8)             # air marked
ax.annotate(f"air: K = γp = {K_air/1e3:.0f} kPa, c = {c_air:.0f} m/s", (K_air, c_air), (2e4, 2e3), fontsize=9,   # label text and position
            arrowprops=dict(arrowstyle="->", color=COLORS["muted"]))      # arrow from the label to the point
ax.loglog([2.2e9], [c_water], "o", color=COLORS["orange"], ms=8)           # water marked
ax.annotate(f"water: K = 2.2 GPa, c = {c_water:.0f} m/s", (2.2e9, c_water), (1e7, 3e0), fontsize=9,   # label text and position
            arrowprops=dict(arrowstyle="->", color=COLORS["muted"]))      # arrow from the label to the point
ax.annotate("K → ∞: incompressible, c → ∞", (8e10, 9e5), (2e7, 9e5), fontsize=9, color=COLORS["rose"], va="center",   # the limit
            arrowprops=dict(arrowstyle="->", color=COLORS["rose"]))       # arrow pointing to ever larger K
ax.set_xlabel("bulk modulus K = ρ(∂p/∂ρ)_s [Pa]"); ax.set_ylabel("speed of sound c [m/s]")   # axis labels with units
ax.text(3e9, 60, f"time to cross 1 km:\\nair {1000/c_air:.2f} s, water {1000/c_water:.2f} s",   # what c means in practice
        fontsize=9, va="top", bbox=dict(facecolor="white", edgecolor=COLORS["grid"]))
ax.set_ylim(1, 4e6)                                                        # headroom above the lines for the limit label
ax.set_title("Stiffer fluids carry sound faster: c = √(K/ρ)")             # the message
ax.legend(loc="upper left", fontsize=8.5)                                  # legend
plt.show()                                                                 # display
""", see="""
Two parallel straight lines of slope ½ on log–log axes, one for a light and one for a heavy fluid, with air and water
marked, an arrow toward infinite stiffness, and a box giving the time sound needs to cross 1 km in each (about 3 s in air,
under 0.7 s in water — the "count the seconds after lightning" rule).
""", read="""
Along a line (same density) a fluid 100× stiffer carries sound 10× faster. Water is about 15 000× stiffer than air but also
800× denser, so its sound speed is only about √(15 000/800) ≈ 4.4× larger.
""", change="""
Heating air raises p at fixed ρ, hence K = γp and c; for a perfect gas c ∝ √T (C47). Squeezing water harder barely
changes its K, which is why the ocean's sound speed varies by only a few per cent.
""")
note("C37", "The thermal expansion coefficient", r"""
α [1/K] is the same kind of held-fixed derivative, at constant pressure: the fractional density drop per kelvin of
warming. Water is strange: α(20 °C) = +2.07×10⁻⁴ K⁻¹ but α(2 °C) = −3.3×10⁻⁵ K⁻¹ — below about 4 °C water contracts when
warmed, which is why lakes freeze from the top. α returns in the adiabatic lapse rate (C54) and in the Boussinesq
approximation (Ch. 4).
""", equation=r"\alpha \equiv -\frac{1}{\rho}\left(\frac{\partial \rho}{\partial T}\right)_p", ref="1.20")
nb.code("""
alpha_w = ch01.thermal_expansion_coefficient(lambda T, p: ch01.water_density(T), np.array([275.15, 293.15]), 1e5)  # (1.20)
print(f"water: α(2 °C) = {alpha_w[0]:.2e} 1/K, α(20 °C) = {alpha_w[1]:.2e} 1/K")   # negative below 4 °C, positive above
""", explain="""
`thermal_expansion_coefficient` differentiates a density function ρ(T, p) at fixed p and divides by −ρ, Eq. (1.20); the
water density correlation ignores p, which is fine at 1 bar.
""")
nb.md("""
**What would change if…** the fluid were a perfect gas? Then K = γp = γρRT and c = √(γRT) depends on temperature only —
§1.9 makes this and α = 1/T explicit. **Next:** the perfect-gas law itself (C40).
""")

# =====================================================================================================================
# A.9 §1.9 Perfect Gas — C40, C45
# =====================================================================================================================
nb.section("1.9", "Perfect Gas", intro="""
**What is this section about?** The equation of state of air and of most gases at ordinary conditions, and what follows
from it: R from the molecular weight, C_p − C_v = R, the isentropic law, and the speed of sound and expansion coefficient
of a gas. (The isothermal atmosphere and the scale height of book §1.10 are stated here, in C40, because they are just
p = ρRT combined with hydrostatics.)
""")
# ---------------------------------------------------------------------------------------------------- C40
core("C40", "The perfect-gas law in continuum form: p = ρRT",
     "The molecular gas law counts molecules in a box. How does it become a law for the density at a point?")
nb.md("""
#### The problem in plain words
Aircraft performance, weather-balloon soundings, the air density in a drag formula, altimeters, every atmospheric model:
all need density from pressure and temperature. The molecular law is written for a container full of molecules; we need
it for a fluid particle.
""")
nb.md("""
#### The idea
```
p = (N/V) k_B T  →  multiply and divide by m  →  ρ (k_B/m) T  →  m = M_w/A_o  →  ρ (R_u/M_w) T  =  ρ R T
    molecules per volume           continuum density (C06)       per kilomole (P07)      one gas constant per gas
```
T is always the absolute temperature in kelvin here (recap R03 in §1.2).
""")
note("C38", "The molecular perfect-gas law", r"""
For n non-interacting molecules in a volume V at absolute temperature T — the result D34 (C06) derived from wall impacts.
It holds when attractions between molecules are negligible and each molecule's own volume is tiny compared with V/n.
One cubic metre of sea-level air (n = 2.55×10²⁵ molecules at 288.15 K) gives p = 1.013×10⁵ Pa.
""", equation=r"pV = n\,k_B T", ref="1.21")
nb.md("""
> ⚠️ **Common confusion — the symbol n switches meaning.** In §1.4 (P07, C04, D34) n was a *number density* [molecules
> per m³]. In (1.21) and in D11 below, n is a *count* of molecules [–], and the number density is n/V. Same letter,
> different quantity — watch the units.
""")
note("C39", "The constants", """
(CODATA, exact since 2019) k_B = 1.380 649×10⁻²³ J/K; Avogadro's number per kilomole A_o = 6.022 140 76×10²⁶ kmol⁻¹;
universal gas constant R_u = k_BA_o = 8314.46 J kmol⁻¹ K⁻¹; dry air M_w = 28.9644 kg/kmol (US Standard Atmosphere 1976),
so R = R_u/M_w = 287.06 J kg⁻¹ K⁻¹. (Kilomoles: primer P07.)
""")
D("D11", "From the molecular gas law to p = ρRT", ref="1.22",
  goal="Rewrite the molecule-counting gas law as a law for the continuum density, with one constant R per gas.",
  assumptions="""Perfect gas: non-interacting molecules (the start) · Kn ≪ 1, so nm/V is a continuum density (step 3) · one
average molecular weight for a mixture such as air (step 4).""",
  start=(r"pV = n\,k_B T", "n non-interacting molecules in a volume V at temperature T exert the pressure p (1.21); here n "
                           "counts molecules, as in the book, so n/V is the number density (which §1.4 itself called n — see the warning above)."),
  plan=["Divide by V.", "Create the density ρ = nm/V.", "Replace one molecule's mass by M_w/A_o.", "Name the constants."],
  uses=["(1.21) (C38, derived in D34)", "continuum density (C06)",
        "mole, kilomole, molecular weight and Avogadro's number (P07)", "constants k_B, A_o, R_u (C39)"],
  steps=[
      ("Divide by the volume", r"p = \dfrac{n}{V}\,k_B T",
       "Divide both sides by V ≠ 0; pressure is a local quantity and should not depend on how big a container we imagine.",
       "Pressure = molecules per volume × k_B × T."),
      ("Multiply and divide by the molecular mass", r"p = \dfrac{n\,m}{V}\,\dfrac{k_B}{m}\,T",
       "m/m = 1 changes nothing; we do it so that nm/V — mass per volume — appears in the formula.",
       "The same law, rearranged to expose mass."),
      ("Recognise the density", r"p = \rho\,\dfrac{k_B}{m}\,T",
       "nm/V is the total molecular mass per volume, the continuum density of C06 (valid for boxes in the continuum "
       "window, Kn ≪ 1).",
       "Density replaces the molecule count."),
      ("Write one molecule's mass per kilomole", r"\dfrac{k_B}{m} = \dfrac{k_B A_o}{M_w}",
       "One kilomole of molecules (A_o of them) has a mass of M_w kg, so m = M_w/A_o (P07); dividing by a fraction "
       "multiplies by its inverse.",
       "Mass of one molecule = molecular weight ÷ Avogadro's number."),
      ("Substitute into step 3", r"p = \rho\,\dfrac{k_B A_o}{M_w}\,T",
       "Replace k_B/m in step 3 by its equal from step 4; substituting equals keeps the equation true and removes the "
       "mass of a single molecule from the law.",
       "The law now contains only the density, measurable constants and T."),
      ("Name the universal gas constant", r"p = \rho\,\dfrac{R_u}{M_w}\,T,\qquad R_u \equiv k_B A_o",
       "k_BA_o = 1.380649×10⁻²³ J/K × 6.02214076×10²⁶ kmol⁻¹ = 8314.46 J kmol⁻¹ K⁻¹ is the same for every gas; naming it "
       "shortens the law.",
       "One universal constant divided by the gas's molecular weight."),
      ("Name the gas constant of this gas", r"p = \rho R T,\qquad R \equiv R_u/M_w",
       "For a given gas R_u/M_w is a single number (air: 8314.46/28.9644 = 287.06 J kg⁻¹ K⁻¹), so we give it one symbol.",
       "The perfect-gas law per unit mass."),
  ],
  result=(r"p = \rho R T", "pressure = density × gas constant × absolute temperature."),
  interpret="""Every atmospheric calculation in the book closes with this law. It fails when molecules attract each other
or fill a noticeable part of the volume — dense gases, gases near condensation (the van der Waals contrast in C41).""",
  check="""Units: kg m⁻³ · J kg⁻¹ K⁻¹ · K = J m⁻³ = Pa ✓. Number: ρ = 101 325/(287.06 × 288.15) = 1.225 kg/m³ ✓
(`perfect_gas_density`, and the from-scratch constants chain below). A mixture: dry air treated as one gas with
M_w = 28.96 ✓.""",
  traps="mixing mol and kmol (R_u = 8.314 J mol⁻¹ K⁻¹ vs 8314 J kmol⁻¹ K⁻¹, a factor 1000); confusing R_u with R; using °C "
        "for T.")
nb.worked_example("air density on three days", """
- Sea level, 15 °C: ρ = p/(RT) = 101 325/(287.06 × 288.15) = 1.225 kg/m³.
- A hot day, 35 °C: 101 325/(287.06 × 308.15) = 1.145 kg/m³ (6.5 % less).
- A high plateau, 84 kPa and 15 °C: 84 000/(287.06 × 288.15) = 1.016 kg/m³ — aircraft need longer runways there.
""")
nb.code("""
print(ch01.gas_constant(ch01.M_W_AIR))                           # R = R_u / M_w for dry air [J/(kg K)] → 287.058
p_days = np.array([101325.0, 101325.0, 84000.0])                 # pressures of the three days [Pa]
T_days = np.array([288.15, 308.15, 288.15])                      # temperatures [K]
rho_days = ch01.perfect_gas_density(p_days, T_days)              # Eq. (1.22) solved for ρ [kg/m³]
print(rho_days)                                                  # [1.225  1.1455 1.0155]
print(ch01.perfect_gas_state(p=101325.0, T=288.15))              # fill in the missing ρ: (p, ρ, T)
print(ch01.perfect_gas_pressure(1.225, 288.15))                  # and back to p = ρRT [Pa] → 101 327
""", explain="""
1. `gas_constant` divides R_u by the molecular weight (step 7 of D11).
2. `perfect_gas_density` is (1.22) solved for ρ, applied to three days at once.
3. `perfect_gas_state` returns whichever of p, ρ, T is missing; `perfect_gas_pressure` is (1.22) itself.
""")
nb.md("#### From scratch: the constants chain")
nb.check_agree("""
R_mine = ch01.K_B * ch01.N_A_KMOL / ch01.M_W_AIR                  # D11 steps 4–7: R = k_B A_o / M_w [J/(kg K)]
rho_mine = p_days / (R_mine * T_days)                             # ρ = p / (R T) [kg/m³]
assert np.allclose(rho_mine, ch01.perfect_gas_density(p_days, T_days))    # same as the library
n_days = p_days / (ch01.K_B * T_days)                             # the molecular route: n/V = p/(k_B T) [1/m³] (D34)
rho_molecular = n_days * ch01.molecular_mass(ch01.M_W_AIR)        # ρ = (n/V) m [kg/m³]
assert np.allclose(rho_molecular, rho_mine)                       # counting molecules gives the same density
print(rho_mine, rho_molecular)   # both routes: [1.225 1.1455 1.0155] kg/m³
""")
nb.figure("""
T_grid = np.linspace(200, 320, 200)                                        # temperatures [K]
fig, ax = plt.subplots(figsize=(7, 4))                                     # one panel
for p_iso, col in ((50e3, COLORS["blue"]), (70e3, COLORS["teal"]), (84e3, COLORS["orange"]), (101325.0, COLORS["accent"])):   # four pressures [Pa]
    ax.plot(T_grid, ch01.perfect_gas_density(p_iso, T_grid), color=col, label=f"p = {p_iso/1e3:.1f} kPa")   # isobar ρ(T)
ax.plot(T_days, rho_days, "o", color=COLORS["ink"], ms=8, label="the three days")                        # worked example
ax.set_xlabel("temperature T [K]"); ax.set_ylabel("density ρ [kg/m³]")                                   # labels with units
ax.set_title("At fixed pressure, warmer air is thinner: ρ = p/(RT)")                                     # the message
ax.legend(fontsize=8.5)                                                                                  # legend
plt.show()                                                                                               # display
""", see="""
Four falling curves, one per pressure (50, 70, 84 and 101.3 kPa), with the three worked-example days marked as black dots:
two on the sea-level curve, one on the 84 kPa curve.
""", read="""
Along one curve (fixed p) the density is inversely proportional to T. The vertical spacing between curves at one T shows
ρ ∝ p. The plateau day sits on the 84 kPa curve, below the two sea-level days.
""", change="""
Humid air has a smaller mean molecular weight (water vapour 18 < 28.96), hence a larger R and a lower density at the same
p and T — humid air is lighter than dry air.
""")
P("P41", "plotly 3-D surface", """
`go.Figure(go.Surface(x=…, y=…, z=…))` draws a surface you can rotate with the mouse; `go.Scatter3d` adds points or lines
on top. The published page keeps it interactive without Python.
""", code="""
import plotly.graph_objects as go                          # plotly's figure objects
xg, yg = np.meshgrid(np.linspace(-1, 1, 10), np.linspace(-1, 1, 10))   # a 10 × 10 grid of (x, y) points
fig = go.Figure(go.Surface(x=xg, y=yg, z=xg**2 + yg**2, showscale=False))   # the paraboloid z = x² + y²
fig.update_layout(height=320, margin=dict(l=0, r=0, t=10, b=0))            # compact
fig.show()                                                                  # drag to rotate
""")
nb.plotly("""
vg, Tg = np.meshgrid(np.linspace(0.5, 2.0, 40), np.linspace(200, 350, 40))        # a 40 × 40 grid of states (v, T)
pg = ch01.perfect_gas_pressure(1 / vg, Tg) / 1e3                                  # p = ρRT with ρ = 1/v [kPa]
fig = go.Figure(go.Surface(x=vg, y=Tg, z=pg, colorscale="Viridis", opacity=0.85, showscale=False, name="p(v, T)"))  # the surface
for kind in routes:                                                               # the three C25 routes as curves on it
    path = ch01.process_path(kind, (v1, 300.0), (2 * v1, 300.0), n=81)            # 81 samples per route
    fig.add_trace(go.Scatter3d(x=path["v"], y=path["T"], z=path["p"] / 1e3, mode="lines",   # a 3-D line on the surface
                               line=dict(color=colors[kind], width=6), name=routes[kind][:1]))   # route colour and label
fig.add_trace(go.Scatter3d(x=[v1], y=[300.0], z=[100.0], mode="markers", marker=dict(color="red", size=6), name="state 1"))  # the dot
fig.update_layout(height=520, margin=dict(l=0, r=0, t=40, b=0), title="Every state of air is a point on p = RT/v",   # size, title
                  scene=dict(xaxis_title="v [m³/kg]", yaxis_title="T [K]", zaxis_title="p [kPa]"))                 # axis labels
fig.show()                                                                        # rotate and hover
""", explain="""
1. `perfect_gas_pressure` evaluates p = ρRT on a grid of specific volumes and temperatures: the surface of all equilibrium
   states of air.
2. The three routes of C25 are drawn on the surface: a process is a curve on it, and "two properties fix the state"
   (C28) means that choosing (v, T) picks exactly one point.
""")
see_read_change(
    'A curved surface of pressure over the plane of specific volume and temperature, the three coloured routes of C25 lying on it, and a red dot for state 1.',
    'Height above the (v, T) floor is the pressure. Moving along the T direction at fixed v the surface rises in a straight line (p ∝ T); moving along v at fixed T it falls like 1/v. Every route between states stays on the surface, because every equilibrium state does.',
    'A van der Waals gas (C41) would give a slightly different surface that sags at moderate v, where molecular attraction lowers the pressure (very close to the volume b taken up by the molecules themselves it would instead shoot up).',
)
note("C41", "A perfect gas stores energy in temperature only", r"""
For a perfect gas the internal energy and the enthalpy depend on **temperature only**, and conversely only a perfect gas
has this property (stated; the proof combines the Gibbs relations with p = ρRT to show (∂e/∂v)_T = 0). A real gas with
molecular attractions does not: for CO₂ modelled as a van der Waals gas, e = C_vT − a/v, so squeezing it from 0.1 to
0.01 m³/kg at a fixed 300 K lowers e by 16.9 kJ/kg (code below).
""", equation=r"e = e(T), \qquad h = h(T)")
note("C48", "The expansion coefficient of a perfect gas", r"""
At fixed p, ρ = p/(RT) ∝ 1/T, so (stated) the expansion coefficient (1.20) is simply 1/T: 3.3×10⁻³ K⁻¹ for air at 300 K —
sixteen times that of water at 20 °C.
""", equation=r"\alpha = \frac{1}{T}", ref="1.28")
note("C62", "The isothermal atmosphere (book §1.10)", r"""
If a layer of atmosphere had one temperature T, hydrostatics (1.8) with ρ = p/(RT) would give dp/p = −(g/RT)dz, and
(stated) pressure falls exponentially with height. Why is one temperature a fair first guess? Over roughly the lowest
50 km the standard atmosphere's temperature stays within a band of about ±15 % around 250 K (the code below finds a largest
deviation of about 15 %), so read the isothermal model as an **order-of-magnitude** description, not a close fit. At 5 km
with T = 250 K our model gives p = 51.2 kPa, while the US Standard Atmosphere has 54.0 kPa — 5 % more, because the real
lower atmosphere is warmer than 250 K.
""", equation=r"p(z) = p_0\,e^{-gz/RT}")
note("C63", "The scale height (book §1.10)", r"""
The height over which the isothermal pressure falls by a factor e (one e-folding, P36) is the **scale height**. For
T = 250 K our `scale_height` gives 287.06 × 250/9.807 ≈ 7.32 km — a good one-number thickness of the atmosphere.
""", equation=r"H = \frac{RT}{g}")
nb.plotly("""
print(f"C63 scale height at 250 K: {ch01.scale_height(250.0):.0f} m")                            # H = RT/g
print(f"C62 isothermal p at 5 km, 250 K: {ch01.isothermal_pressure(5000.0, ch01.P_ATM, 250.0):.0f} Pa")  # p0 e^(−gz/RT)
print("C62 USSA-1976 (T, p, ρ) at 5 km geopotential height:", ch01.standard_atmosphere(5000.0))       # the standard atmosphere
T_50 = ch01.standard_atmosphere(np.linspace(0, 50e3, 501))[0]                                     # USSA temperatures, 0–50 km [K]
print(f"C62 largest deviation of the USSA temperature from 250 K below 50 km: {np.max(np.abs(T_50/250 - 1)):.1%}")  # ≈ 15 %
print(f"C48 α of air at 300 K: {ch01.perfect_gas_expansion_coefficient(300.0):.2e} 1/K")          # 1/T
de_vdw = (ch01.van_der_waals_internal_energy(300.0, 0.1, ch01.VDW_CO2["a"], 655.0)                # C41: e at v = 0.1 m³/kg …
          - ch01.van_der_waals_internal_energy(300.0, 0.01, ch01.VDW_CO2["a"], 655.0))            # … minus e at 0.01 m³/kg [J/kg]
print(f"C41 van der Waals CO₂ squeezed 0.1 → 0.01 m³/kg at 300 K: e drops by {de_vdw:.0f} J/kg")  # a real gas: e depends on v
z_iso = np.linspace(0, 40e3, 161)                                                                  # heights [m]
p_ussa = ch01.standard_atmosphere(z_iso)[1] / 1e3                                                 # USSA pressure [kPa]
def iso(T):                                                                                        # one slider position
    H = ch01.scale_height(T)                                                                       # scale height [m]
    return {"isothermal p(z)": (z_iso / 1e3, ch01.isothermal_pressure(z_iso, ch01.P_ATM, T) / 1e3),   # the model [km, kPa]
            "US Standard Atmosphere 1976": (z_iso / 1e3, p_ussa),                                   # the reference, fixed
            "p₀/e at z = H": ([H / 1e3], [ch01.P_ATM / np.e / 1e3])}                                # one e-folding up
fig = slider_figure(iso, "T", np.linspace(200, 300, 21), unit="K", xlabel="height z [km]", ylabel="pressure p [kPa]",  # 21 steps
                    title="One temperature gets the shape of p(z) roughly right", modes={"p₀/e at z = H": "markers"}, active=10)   # start at 250 K
fig.show()                                                                                         # drag T
""", explain="""
1. The printed lines are the numbers quoted in the C41, C48, C62 and C63 notes, computed now.
2. For each slider temperature the figure shows the isothermal model p₀e^(−gz/RT) (`isothermal_pressure`), the fixed US
   Standard Atmosphere profile, and a marker where the model has fallen to p₀/e — at z = H.
""")
see_read_change(
    "A falling exponential p(z) for the chosen temperature, the fixed US Standard Atmosphere profile, and a marker where the model's pressure has dropped to p₀/e.",
    "The marker's height is the scale height H = RT/g. Where the two curves overlap, one temperature describes the real atmosphere well; where they separate, the real temperature profile matters.",
    'At 200 K the model falls too quickly (H ≈ 5.9 km), at 300 K too slowly (H ≈ 8.8 km): the scale height is proportional to T.',
)
nb.md("""
**What would change if…** the gas were squeezed too quickly for heat to escape? Its temperature would rise, and p would
climb faster than ρ. **Next:** the isentropic law (C45).
""")
# ---------------------------------------------------------------------------------------------------- C45
core("C45", "The isentropic law: p/ρ^γ = const",
     "Squeeze a gas quickly — no time for heat to leak out and no friction. How does its pressure follow its density?")
nb.md("""
#### The problem in plain words
A bicycle pump's barrel heats up; a diesel engine ignites its fuel by compression alone; a rising air parcel cools as it
expands (§1.10); sound is a rapid squeeze. All follow the same power law.
""")
nb.md("""
#### The idea
| squeeze to 2× the density | isothermal (heat leaks out) | isentropic (no heat, no friction) |
|---|---|---|
| pressure ratio | 2 | 2^1.4 = 2.64 |
| temperature | unchanged | × 2^0.4 = 1.32 |

**The work done on the gas stays in it as internal energy, so the pressure climbs faster.** On log–log axes: slope 1
versus slope γ.
""")
note("C42", "R = C_p − C_v", r"""
For a perfect gas h = e + pv = e + RT, and e and h depend on T only (C41), so the partial derivatives in (1.14)–(1.15)
become ordinary ones (*gloss:* a total derivative d/dT, because nothing else varies). Differentiating h = e + RT with
respect to T gives (stated) the gap between the specific heats: 1004.7 − 717.6 = 287.1 J kg⁻¹ K⁻¹.
""", equation=r"R = C_p - C_v", ref="1.23")
note("C43", "The ratio of specific heats", """
Kinetic theory with rigid molecules gives γ = 5/3 for monatomic gases (He, Ar), 7/5 for linear molecules — all diatomic
ones (N₂, O₂, air) and also linear CO₂ — and 4/3 for non-linear molecules such as H₂O. Real CO₂ comes out lower (about 1.3),
because its molecular vibrations already store energy at room temperature. `ch01.GAMMA_BY_ATOMICITY` holds the rigid values.
""", equation=r"\gamma \equiv C_p/C_v", ref="1.24")
note("N16", "Air's values", """
Air at ordinary temperatures: γ = 1.40 and C_p ≈ 1005 J kg⁻¹ K⁻¹ (our constant 1004.7 follows from γ = 1.4 and R). Both
C_p and C_v rise slowly with temperature; we treat them as constants, which is what (1.25) needs.
""")
note("C44", "Adiabatic versus isentropic", """
**Adiabatic** = no heat crosses the boundary. **Isentropic** = adiabatic *and* frictionless, so s stays constant. Stirring an
insulated gas is adiabatic but not isentropic (Δs = 23.5 J kg⁻¹ K⁻¹ in the C34 example).
""")
P("P42", "separation of variables", r"""
For an equation like dy/y = k dx, put each variable on its own side and integrate both sides: ln y = kx + C, so
y = y₀e^{kx}. For dp/p = γ dρ/ρ the same move gives ln p = γ ln ρ + C.
""", code="""
print(np.log(2.0**1.4), 1.4 * np.log(2.0))   # ln(2^1.4) = 1.4 ln 2 = 0.9704: integrating dp/p = 1.4 dρ/ρ from ρ to 2ρ
""")
nb.md("""
> 🔁 Exponent rules — a^m a^n = a^(m+n), (a^m)^n = a^(mn), a^(−m) = 1/a^m, (a/b)^m = a^m/b^m — were primed in C06 (P43).
""")
D("D14", "The isentropic law for a perfect gas", ref="1.25",
  goal="Find how pressure follows density when a perfect gas is compressed or expanded with no heat and no friction.",
  assumptions="""Isentropic — adiabatic and frictionless, ds = 0 (steps 1, 3) · perfect gas (steps 2, 4) · constant C_p and
C_v, hence constant γ (step 10).""",
  start=(r"T\,ds = de + p\,dv = dh - v\,dp", "the Gibbs relations (1.18), valid for any process."),
  plan=["Set ds = 0 in both forms.", "Use e = e(T) and h = h(T) to bring in C_v and C_p.",
        "Divide the two equations so that T drops out.", "Switch from v to ρ and integrate."],
  uses=["Gibbs relations (C35, D10)", "e = e(T), h = h(T) for a perfect gas (C41)", "C_p, C_v (C30, C31)",
        "γ = C_p/C_v (C43)", "differentials (P34)", "separation of variables (P42)", "natural logarithm (P36)",
        "exponent rules (P43)"],
  steps=[
      ("Set ds = 0 in the first Gibbs form", r"0 = de + p\,dv",
       "Isentropic means s does not change, so T ds = 0 (T > 0). We start with the energy form because the work p dv "
       "appears in it.",
       "All the work done on the gas goes into internal energy."),
      ("Use e = e(T)", r"C_v\,dT = -p\,dv",
       "For a perfect gas e depends on T only (C41), so the partial derivative in (1.15) is an ordinary one and de = C_v dT; "
       "then move p dv to the other side.",
       "Compression (dv < 0) warms the gas."),
      ("Set ds = 0 in the second Gibbs form", r"0 = dh - v\,dp",
       "The same isentropic condition T ds = 0 applied to the enthalpy form of (1.18). We need a second equation, one that "
       "contains dp.",
       "At constant entropy, enthalpy rises by v dp when the pressure rises."),
      ("Use h = h(T)", r"C_p\,dT = v\,dp",
       "For a perfect gas h depends on T only (C41), so the partial derivative in (1.14) is an ordinary one and dh = C_p dT; "
       "then move v dp to the other side.",
       "Raising the pressure warms the gas."),
      ("Divide step 4 by step 2", r"\dfrac{C_p}{C_v} = \dfrac{v\,dp}{-\,p\,dv}",
       "Along a compression dT ≠ 0, so both sides of step 2 are nonzero and we may divide; dT cancels, which removes the "
       "temperature from the problem.",
       "The ratio of specific heats links the pressure change to the volume change."),
      ("Name the ratio γ", r"\gamma = -\dfrac{v\,dp}{p\,dv}",
       "γ ≡ C_p/C_v by definition (1.24); the minus sign is moved to the front of the fraction.",
       "γ measures how much faster p changes than v, in relative terms."),
      ("Separate the variables", r"\dfrac{dp}{p} = -\gamma\,\dfrac{dv}{v}",
       "Multiply both sides by −dv/v; now each side contains only one variable (P42).",
       "A 1 % decrease in volume raises the pressure by γ %."),
      ("Switch from specific volume to density", r"\dfrac{dv}{v} = -\dfrac{d\rho}{\rho}",
       "v = 1/ρ gives dv = −dρ/ρ² (the derivative of 1/ρ); dividing by v = 1/ρ gives −dρ/ρ. Watch the sign: this is the "
       "classic slip.",
       "A 1 % increase in density is a 1 % decrease in volume."),
      ("Substitute into step 7", r"\dfrac{dp}{p} = \gamma\,\dfrac{d\rho}{\rho}",
       "Replace dv/v in step 7 by −dρ/ρ from step 8; the two minus signs cancel.",
       "Relative pressure change = γ × relative density change."),
      ("Integrate both sides with constant γ", r"\ln p = \gamma\,\ln\rho + \text{const}",
       "∫dx/x = ln x (P36); γ can be taken outside the integral because C_p and C_v are constant (assumption).",
       "On log–log axes p against ρ is a straight line of slope γ."),
      ("Combine the logarithms and exponentiate", r"\dfrac{p}{\rho^\gamma} = \text{const}",
       "ln p − γ ln ρ = ln(p/ρ^γ) by the log and exponent rules (P43); if a logarithm is constant, so is its argument.",
       "Along an isentrope p/ρ^γ never changes."),
  ],
  result=(r"p/\rho^\gamma = \text{const}",
          "with no heat and no friction, a perfect gas whose C_p and C_v do not vary keeps p ∝ ρ^γ."),
  interpret="""Work done on the gas has nowhere to go but internal energy, so the gas heats up and its pressure rises faster
than in an isothermal squeeze (slope γ instead of 1). Combined with p = ρRT it gives the ratios (1.26), hence θ (D20), and
c = √(γRT) (C47). It fails with friction or heat exchange (stirring, slow compression in a conducting cylinder) and when
C_p varies strongly with T. In the explainer `heat_work_paths` below, the bold isentropic route lies exactly on the faint
isentrope.""",
  check="""Units: the constant carries Pa (m³/kg)^γ — fine, it is fixed by the starting state ✓. Limit γ → 1: p/ρ = const, the
isothermal law (RT constant) ✓. Number: a pump from 1 to 2 bar, ρ₂/ρ₁ = 2^(1/1.4) = 1.641 and T₂ = 288.15 × 2/1.641 = 351.3 K ✓;
the step-by-step integration of dp/dρ = γp/ρ below agrees with `isentropic_pressure` to 10⁻⁴.""",
  traps="the sign of dv/v = −dρ/ρ; assuming constant C_p and C_v silently; dividing step 2 by step 4 and getting 1/γ.")
nb.worked_example("a bicycle pump from 1 bar to 2 bar", """
Start at 288.15 K.
1. Density ratio from (1.25): ρ₂/ρ₁ = (p₂/p₁)^(1/γ) = 2^0.714 = 1.641.
2. Temperature from p = ρRT: T₂/T₁ = (p₂/p₁)(ρ₁/ρ₂) = 2/1.641 = 1.219 → T₂ = 351.3 K, 63 K hotter.
3. Isothermal comparison: the density ratio would be 2 and there would be no warming at all.
""")
nb.code("""
T_ratio, rho_ratio = ch01.isentropic_ratios(2.0)                  # Eq. (1.26) for p2/p1 = 2: (T2/T1, ρ2/ρ1)
print(T_ratio, rho_ratio, 288.15 * T_ratio)                      # 1.2190 1.6407 → T2 = 351.26 K
rho1 = ch01.perfect_gas_density(1e5, 288.15)                     # starting density at 1 bar, 15 °C [kg/m³]
print(ch01.isentropic_pressure(rho_ratio * rho1, 1e5, rho1))     # Eq. (1.25): the pressure at the new density → 2.0e5 Pa
print(ch01.cv_from_cp(ch01.CP_AIR), ch01.gamma_from_cp(ch01.CP_AIR))   # (1.23)–(1.24): C_v = 717.64, γ = 1.4
print(ch01.perfect_gas_sound_speed(np.array([250.0, 288.15, 300.0])))  # (1.27): c at three temperatures [m/s]
""", explain="""
1. `isentropic_ratios` applies (1.26) to a pressure ratio of 2 — the pump numbers.
2. `isentropic_pressure` is (1.25): starting from (ρ₁, p₁), it returns p at any density on the same isentrope.
3. `cv_from_cp` and `gamma_from_cp` use R = C_p − C_v and γ = C_p/C_v.
4. `perfect_gas_sound_speed` is c = √(γRT) of C47, below.
""")
nb.md("#### From scratch: follow the isentrope in small steps")
nb.check_agree("""
n_steps = 20_000                                                   # Euler steps (primer P30)
rho_grid = np.linspace(rho1, rho_ratio * rho1, n_steps + 1)        # densities from ρ1 to 1.641 ρ1 [kg/m³]
p_mine = np.empty(n_steps + 1); p_mine[0] = 1e5                    # start at 1 bar
for k in range(n_steps):                                           # step dp = γ p dρ/ρ (D14 step 9)
    p_mine[k + 1] = p_mine[k] + ch01.GAMMA_AIR * p_mine[k] / rho_grid[k] * (rho_grid[k + 1] - rho_grid[k])   # p_{k+1} = p_k + γ p_k Δρ/ρ_k
assert np.allclose(p_mine, ch01.isentropic_pressure(rho_grid, 1e5, rho1), rtol=1e-4)   # matches (1.25) everywhere
print(f"Euler end pressure {p_mine[-1]:.1f} Pa vs isentropic law {ch01.isentropic_pressure(rho_grid[-1], 1e5, rho1):.1f} Pa")   # hand integration vs (1.25)
""")
nb.figure("""
r = np.logspace(np.log10(0.5), np.log10(3.0), 100)                             # density ratios ρ/ρ0 [–]
fig, ax = plt.subplots(figsize=(6.5, 4.2))                                     # one log–log panel
ax.loglog(r, ch01.isentropic_pressure(r, 1.0, 1.0), color=COLORS["accent"], label="isentropic: p ∝ ρ^1.4 (no heat, no friction)")  # (1.25)
ax.loglog(r, r, color=COLORS["muted"], label="isothermal: p ∝ ρ (heat leaks freely)")          # p/p0 = ρ/ρ0
ax.loglog([1, rho_ratio], [1, 2], "o", color=COLORS["orange"], ms=8, label="pump: 1 → 2 bar")  # the worked example
ax.plot([1.5, 2.5, 2.5], [1.5**1.4, 1.5**1.4, 2.5**1.4], color=COLORS["accent"], lw=1)         # slope triangle
ax.text(2.55, 1.9**1.4 * 1.05, "slope 1.4", color=COLORS["accent"], fontsize=9)                # its label
ax.set_xlabel("density ratio ρ/ρ₀ [–]"); ax.set_ylabel("pressure ratio p/p₀ [–]")             # dimensionless axes
from matplotlib.ticker import FuncFormatter                                                   # control how tick labels are written
for axis in (ax.xaxis, ax.yaxis):                                                              # both log axes
    axis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))                              # plain numbers: 1, 2 instead of 10⁰, 2×10⁰
    axis.set_minor_formatter(FuncFormatter(lambda v, _: f"{v:g}" if v in (0.5, 2, 3, 5) else ""))  # label only a few minor ticks
ax.set_title("On log–log axes the exponent is the slope")                                      # the message
ax.legend(loc="upper left", fontsize=8.5)                                                      # legend
plt.show()                                                                                     # display
""", see="""
Two straight lines through the point (1, 1): a steeper purple isentrope and a grey isotherm, with the pump's start and end
states on the isentrope and a slope triangle.
""", read="""
The slope of each line is the exponent: 1.4 for a squeeze with no heat exchange, 1 when heat leaks away freely. Doubling
the pressure needs only a 1.64× denser gas on the isentrope, because the gas heats up.
""", change="""
Helium (γ = 5/3) gives an even steeper line — it heats more when pumped. A slow squeeze in a thin metal cylinder would move
the gas toward the grey line.
""")
note("C46", "The isentropic ratios", r"""
Combining (1.25) with p = ρRT gives (stated) the temperature and density along an isentrope directly from the pressure
ratio — the pump numbers above; θ in §1.10 is built from the first of them.
""", equation=r"\frac{T}{T_0} = \left(\frac{p}{p_0}\right)^{(\gamma-1)/\gamma}, \qquad \frac{\rho}{\rho_0} = \left(\frac{p}{p_0}\right)^{1/\gamma}",
     ref="1.26")
note("C47", "Sound speed of a perfect gas", r"""
Differentiating p = Kρ^γ at constant s gives (∂p/∂ρ)_s = γp/ρ = γRT, so (stated, closing C36) the speed of sound depends
on temperature only: 340.3 m/s at 288.15 K, 317 m/s at 250 K. **Trap:** Newton used the isothermal slope RT and got
√(RT) = 287.6 m/s — 15 % too slow.
""", equation=r"c = \sqrt{\gamma R T}", ref="1.27")
nb.explainer("heat_work_paths", heading="Same two states: what depends on the path?", why="""
Path dependence clicks only when you drag a route between the same two states and watch q and w change while Δe and Δs stay
fixed; the irreversible modes then show where T ds = δq breaks. It ties together C25, C35 and this block, and its
Derivation tab steps through D10 and D14.
""", tries=[
    "Run 'isothermal ×2', then 'cool then expand', and compare the q and w bars.",
    "Open the Derivation tab and step through D10 with the stirring mode switched on at step 7.",
    "Choose 'pump (isentropic)' and check that the bold curve lies on the faint isentrope.",
    "Free expansion: which bar stays at zero, and what does Δs do?",
])
nb.md("""
**What would change if…** a parcel of air were carried upward instead of pumped? Its pressure falls, so it expands and
cools along exactly this isentrope — while the air around it follows hydrostatics. Whether it ends up heavier or lighter
than its surroundings is the question of §1.10.
""")

# =====================================================================================================================
# A.10 §1.10 Stability of Stratified Fluid Media — C50, C51, C54, C55
# =====================================================================================================================
nb.section("1.10", "Stability of Stratified Fluid Media", intro="""
**What is this section about?** Whether a fluid whose density changes with height is stable: displace a small parcel and
see whether buoyancy pushes it back. For air this becomes a comparison of two lapse rates — written with **two opposite
sign conventions** in the literature, and both are used here — and it is captured in one line by the potential
temperature. (This book section also contains the isothermal atmosphere and the scale height, already stated as C62 and
C63 in §1.9.)

> **Lapse-rate conventions in every cell of this section.** We compute with Kundu's Γ ≡ dT/dz, print the meteorology
> value Γ_met = −dT/dz next to it, and always write the stability inequality in the form being used.
""")
# ---------------------------------------------------------------------------------------------------- C50
core("C50", "The displaced parcel: buoyancy as a restoring or runaway force",
     "Nudge a small blob of fluid upward and let go. What force acts on it, and does it come back?")
nb.md("""
#### The problem in plain words
Morning fog stays trapped in a valley; a chimney plume spreads out flat under an inversion (a layer where the temperature rises with height); the ocean's thermocline (the layer, typically a few hundred metres down, where temperature and density change fastest with depth) blocks
mixing between warm surface water and the cold deep; on a hot afternoon air parcels shoot up into towering cumulus. The
same test decides all of them.
""")
nb.md("""
#### The idea
```
z_o + ζ ──  parcel density  ρ(z_o) + (dρ_a/dz)·ζ        environment  ρ(z_o) + (dρ/dz)·ζ
   ↑ ζ       the parcel changes along ITS OWN path        the surroundings keep THEIR profile
z_o ─────── parcel = environment (at rest)
parcel heavier than its new surroundings  → pushed back down  (stable)
parcel lighter                            → pushed further up (unstable)
```
The tools are already in place: Newton's second law (P09, C06), the first-order Taylor expansion (P26, C20) and buoyancy =
weight of the displaced fluid (D37, C20).
""")
P("P45", "square root of a negative number", r"""
There is no real number whose square is negative, so we name one: $i$ with $i^2 = -1$, and then $\sqrt{-a} = i\sqrt a$.
Two facts turn exponentials into the curves we will see. **Euler's formula** $e^{i\theta} = \cos\theta + i\sin\theta$, so
$\tfrac12(e^{iNt} + e^{-iNt}) = \cos Nt$ — an oscillation. The **hyperbolic cosine** $\cosh x = \tfrac12(e^{x} + e^{-x})$ —
equal to 1 at x = 0, then growing like ½eˣ: a runaway. Its partner, the **hyperbolic sine** $\sinh x = \tfrac12(e^{x} - e^{-x})$, starts at 0 and also grows like ½eˣ. Replacing a real frequency by an imaginary one turns cos into cosh:
$\cos(i\sigma t) = \cosh(\sigma t)$.
""", code="""
theta = 0.7                                                    # any angle [rad]
print(np.exp(1j * theta), np.cos(theta) + 1j * np.sin(theta))  # Euler's formula (1j is Python's i): both (0.7648+0.6442j)
print(np.cosh(2.0), (np.exp(2.0) + np.exp(-2.0)) / 2)          # cosh 2 = (e² + e⁻²)/2 = 3.762
print(np.cos(1j * 2.0).real)                                   # cos(2i) = cosh 2 = 3.762: an imaginary frequency means growth
""")
P("P44", "linear second-order ODE", r"""
$\zeta'' + N^2\zeta = 0$ with a constant $N^2$. Try $\zeta = e^{\lambda t}$: every derivative brings down a factor λ, so
$\lambda^2 + N^2 = 0$ and $\lambda = \pm\sqrt{-N^2}$. If $N^2 > 0$ the roots are $\pm iN$ and (P45) the solutions are
cos Nt and sin Nt — oscillation. If $N^2 < 0$ the roots are real, $\pm\sigma$ with $\sigma = \sqrt{-N^2}$, and the
solutions are cosh σt and sinh σt — runaway. If $N^2 = 0$ both roots are 0 (a *double* root): the equation is just ζ″ = 0,
whose solutions are 1 and t, so ζ = A + Bt — a straight line.
""", code="""
N_demo, t_demo = 0.01, 100.0                                           # a frequency [1/s] and a time [s]
zeta_dd = -5 * N_demo**2 * np.cos(N_demo * t_demo)                     # second derivative of ζ = 5 cos(N t)
print(zeta_dd + N_demo**2 * 5 * np.cos(N_demo * t_demo))               # ζ'' + N² ζ → 0.0: it solves the equation
""")
D("D18", "The parcel equation of §1.10 and N²", ref="1.29",
  goal="""Find the equation of motion of a small blob of fluid pushed a little way up or down in a stratified fluid at
rest, and the number that decides whether it comes back.""",
  assumptions="""Small displacement ζ (steps 5–8) · the parcel moves without friction and adiabatically, with no mixing
(step 5: it follows its own isentropic density) · its pressure instantly equals that of its surroundings (step 2: buoyancy
from the environment) · the background is static (step 6) · released from rest (steps 13–14).""",
  start=(r"\rho_p V\,\dfrac{d^2\zeta}{dt^2} = \sum F_z",
         "Newton's second law for a parcel of volume V and density ρ_p displaced by ζ from its rest height z_o."),
  plan=["Forces: weight and buoyancy.", "Densities of the parcel and of its new surroundings, to first order in ζ.",
        "Keep only first-order terms.", "Name the coefficient N² and solve."],
  uses=["Newton's second law (P09)", "buoyancy (D37, C20)", "first-order Taylor expansion (P26)",
        "square root of a negative number, Euler's formula and cosh (P45)", "linear second-order ODE (P44)"],
  steps=[
      ("Name the parcel's two vertical forces", r"\sum F_z = -\rho_p V g + F_b",
       "Its weight points down; the surrounding fluid's pressure gives a net force F_b. Nothing else acts on a frictionless "
       "parcel. We list the forces before using Newton's law.",
       "Weight down, buoyancy up."),
      ("Use the buoyancy of D37 at the new height", r"F_b = \rho(z_o + \zeta)\,V g",
       "The parcel's pressure matches its surroundings, so the pressure around it is the environment's hydrostatic pressure "
       "at z_o + ζ; D37 says the net push is the weight of displaced *environment* fluid there.",
       "The parcel is buoyed up by the fluid it pushes aside at its new height."),
      ("Insert both forces into Newton's law", r"\rho_p V\,\dfrac{d^2\zeta}{dt^2} = -\rho_p V g + \rho(z_o+\zeta)\,V g",
       "The start says mass × acceleration equals the sum of the forces; steps 1 and 2 give that sum term by term.",
       "The parcel accelerates according to the difference between buoyancy and weight."),
      ("Divide by the parcel's mass ρ_pV", r"\dfrac{d^2\zeta}{dt^2} = -\,g\,\dfrac{\rho_p - \rho(z_o+\zeta)}{\rho_p}",
       "ρ_pV ≠ 0, so we may divide both sides by it; the g terms are collected as −g(ρ_p − ρ(z_o+ζ))/ρ_p, which leaves the "
       "acceleration alone on the left.",
       "A parcel heavier than its surroundings accelerates down."),
      ("Expand the parcel's density along its own path", r"\rho_p = \rho(z_o) + \dfrac{d\rho_a}{dz}\,\zeta",
       "First-order Taylor expansion (P26). The parcel started with the environment's density ρ(z_o) and changes "
       "adiabatically as it meets new pressure; dρ_a/dz is that isentropic rate of change.",
       "The parcel's density changes at its own rate."),
      ("Expand the environment's density at the new height", r"\rho(z_o+\zeta) = \rho(z_o) + \dfrac{d\rho}{dz}\,\zeta",
       "First-order Taylor expansion of the static background profile; dρ/dz is the environment's gradient.",
       "The surroundings have their own, different rate."),
      ("Subtract the two", r"\rho_p - \rho(z_o+\zeta) = \Big(\dfrac{d\rho_a}{dz} - \dfrac{d\rho}{dz}\Big)\zeta",
       "ρ(z_o) cancels because parcel and environment matched at rest; what remains is first order in ζ.",
       "The density difference grows in proportion to the displacement."),
      ("Keep only the leading term in the denominator", r"\rho_p \approx \rho(z_o)",
       "The numerator is already proportional to ζ; the correction (dρ_a/dz)ζ in the denominator would only add terms of "
       "order ζ², which we drop for a small displacement.",
       "For the inertia, the parcel's density is simply its starting density."),
      ("Substitute steps 7 and 8 into step 4",
       r"\dfrac{d^2\zeta}{dt^2} = -\,\dfrac{g}{\rho(z_o)}\Big(\dfrac{d\rho_a}{dz} - \dfrac{d\rho}{dz}\Big)\zeta",
       "Replace the density difference in the numerator by step 7 and the denominator by step 8; nothing else changes.",
       "The acceleration is proportional to the displacement."),
      ("Swap the order in the bracket",
       r"\dfrac{d^2\zeta}{dt^2} - \dfrac{g}{\rho(z_o)}\Big(\dfrac{d\rho}{dz} - \dfrac{d\rho_a}{dz}\Big)\zeta = 0",
       "−(a − b) = (b − a); then move the right-hand side to the left. This is the form the book writes.",
       "The parcel equation of §1.10."),
      ("Name the coefficient N²", r"N^2 \equiv -\dfrac{g}{\rho(z_o)}\Big(\dfrac{d\rho}{dz} - \dfrac{d\rho_a}{dz}\Big)",
       "Definition (1.29), chosen so that the equation reads ζ″ + N²ζ = 0 — the oscillator form of P44.",
       "One number collects gravity and the two density gradients."),
      ("Try an exponential solution", r"\zeta = e^{\lambda t}\ \Rightarrow\ \lambda^2 + N^2 = 0",
       "For a linear equation with constant coefficients, e^{λt} turns derivatives into powers of λ (P44); e^{λt} ≠ 0 "
       "divides out.",
       "The motion is set by the roots λ = ±√(−N²)."),
      ("Apply release from rest when N² ≠ 0",
       r"\zeta = \zeta_0\cos Nt\ \ (N^2>0),\qquad \zeta = \zeta_0\cosh\!\big(\sqrt{-N^2}\,t\big)\ \ (N^2<0)",
       "With two different roots ±λ the general solution is ae^{λt} + be^{−λt}. Starting at rest, ζ′(0) = λ(a − b) = 0, so "
       "a = b = ζ₀/2 (from ζ(0) = ζ₀): the even combination ½ζ₀(e^{λt} + e^{−λt}). For imaginary λ = iN that is ζ₀ cos Nt "
       "(Euler's formula, P45); for real λ = √(−N²) it is ζ₀ cosh √(−N²)t (P45).",
       "Stable columns make the parcel oscillate; unstable ones make it run away."),
      ("Treat the double root N² = 0", r"\zeta = A + Bt,\quad \zeta'(0) = B = 0\ \Rightarrow\ \zeta = \zeta_0",
       "With N² = 0 both roots are 0 and e^{λt} gives only the constant; the equation is ζ″ = 0, whose general solution is "
       "the straight line A + Bt (P44). Release from rest means ζ′(0) = B = 0, and ζ(0) = A = ζ₀.",
       "A neutral column leaves the parcel where it was put."),
  ],
  result=(r"\dfrac{d^2\zeta}{dt^2} + N^2\zeta = 0,\qquad N^2 = -\dfrac{g}{\rho(z_o)}\Big(\dfrac{d\rho}{dz} - \dfrac{d\rho_a}{dz}\Big)",
          "a displaced parcel is a harmonic oscillator with squared frequency N², positive when the environment's density "
          "falls with height faster than the parcel's own."),
  interpret="""Stability is a race between two density gradients: the environment's and the parcel's own adiabatic one. It
is the book's first stability analysis — the same "linearise, then look at the sign of the growth rate" logic returns in
Ch. 11. N is the highest frequency internal waves can have (Ch. 7 §7.8). The model has no friction, so a stable parcel rings
forever; real fluids damp it. It fails for large displacements (the Taylor steps) and when the parcel mixes with its
surroundings.""",
  check="""Units: (m s⁻² / kg m⁻³) × kg m⁻⁴ = s⁻² ✓. Limits: a uniform incompressible fluid with an incompressible parcel,
dρ/dz = dρ_a/dz = 0 → N² = 0, neutral ✓; an environment that is exactly isentropic, dρ/dz = dρ_a/dz → neutral ✓. Number:
ocean thermocline, dρ/dz = −0.01 kg m⁻⁴, dρ_a/dz = 0, ρ = 1025 kg/m³ → N² = 9.57×10⁻⁵ s⁻², period 642 s ✓ (the code below:
`brunt_vaisala_sq`, `parcel_displacement`, an Euler–Cromer loop and the unlinearised `parcel_ode_from_gradients` agree).""",
  traps="the sign of buoyancy; using the environment's gradient for the parcel; keeping O(ζ²) terms inconsistently; "
        "replacing ρ_p by ρ(z_o) *before* subtracting (that would lose the whole effect).")
nb.worked_example("a parcel in the ocean thermocline", """
ρ = 1025 kg/m³; the density falls upward by 1 kg/m³ per 100 m, so dρ/dz = −0.01 kg m⁻⁴; treat the parcel as
incompressible, dρ_a/dz = 0.
1. N² = −(g/ρ)(dρ/dz − dρ_a/dz) = −(9.81/1025)(−0.01 − 0) = 9.57×10⁻⁵ s⁻².
2. N = 9.78×10⁻³ s⁻¹; period 2π/N = 642 s ≈ 10.7 min.
3. Released 5 m above its rest height: ζ(t) = 5 cos(0.00978 t) m.
""")
nb.code("""
N2_oc = ch01.brunt_vaisala_sq(1025.0, -0.01, 0.0)                     # Eq. (1.29): thermocline, incompressible parcel [1/s²]
t = np.linspace(0, 1800, 601)                                          # half an hour, every 3 s [s]
zeta = ch01.parcel_displacement(t, 5.0, N2_oc)                         # D18 step 13: ζ = 5 cos(N t) [m]
print(f"N² = {N2_oc:.4g} 1/s², {ch01.stability_timescale(N2_oc)}")    # 9.567e-05 and ('period', ≈ 642 s)
drho_a_demo = -0.005                                                   # a parcel that gets denser as it sinks [kg/m⁴] (illustrative)
N2_demo = ch01.brunt_vaisala_sq(1025.0, -0.01, drho_a_demo)            # its N² from (1.29) [1/s²]
for z0 in (5.0, 300.0, 1000.0):                                        # small, larger and very large releases [m]
    lin = ch01.parcel_displacement(t, z0, N2_demo)                     # linear solution (steps 13–14)
    nl = ch01.parcel_ode_from_gradients(t, z0, 1025.0, -0.01, drho_a_demo)   # step 4 before linearising, solve_ivp (P31)
    print(f"ζ0 = {z0:6.0f} m: largest linear − unlinearised gap = {np.abs(lin - nl).max():.2e} m "
          f"({np.abs(lin - nl).max() / z0:.1e} of the amplitude)")    # grows with the release height
""", explain="""
1. `brunt_vaisala_sq` is (1.29) with the two density gradients; `stability_timescale` gives the period.
2. `parcel_displacement` is the solution of steps 13–14 (cos, constant or cosh, chosen by the sign of N²).
3. `parcel_ode_from_gradients` solves Newton's law of step 4 *before* the small-displacement steps 5–8, keeping the parcel's
   changing density ρ_p in the denominator. (With an incompressible parcel, dρ_a/dz = 0, that denominator never changes and
   the two equations would be identical, so here the parcel is given a density gradient of its own.)
4. The gap relative to the amplitude grows with the release height — from about 10⁻⁵ for 5 m to a few ×10⁻³ for 1 km:
   linearisation (step 8) is excellent for small nudges and slowly degrades for large ones.
""")
nb.md("#### From scratch: an Euler–Cromer loop")
nb.check_agree("""
dt_p = 1.0                                          # time step [s]
z_mine = np.empty(1801); w = 0.0                    # displacement history [m] and vertical velocity [m/s]
z_mine[0] = 5.0                                     # released 5 m up, from rest
for k in range(1800):                               # 1800 steps = 30 min
    w += -N2_oc * z_mine[k] * dt_p                  # velocity first: dw/dt = −N² ζ (the parcel equation)
    z_mine[k + 1] = z_mine[k] + w * dt_p            # then position with the NEW velocity (Euler–Cromer, P30)
assert np.allclose(z_mine[::3], zeta, atol=0.05)    # within 5 cm (1 % of the amplitude) of the exact cosine
print(f"after 30 min: loop {z_mine[-1]:.3f} m, exact {zeta[-1]:.3f} m")   # loop vs cosine
""")
nb.animation("""
nfr_p = 48 if not FAST else 24                                                  # frames
t_an = np.linspace(0, 600, nfr_p)                                               # ten minutes [s]
cases = [(1e-4, "stable, N² = +1e-4 s⁻²", COLORS["teal"]), (0.0, "neutral, N² = 0", COLORS["amber"]),   # (N², label, colour) for …
         (-1e-4, "unstable, N² = −1e-4 s⁻²", COLORS["rose"])]                   # the three columns
top = 200.0                                                                     # height of the drawn column top [m]
tracks = [np.minimum(ch01.parcel_displacement(t_an, 20.0, N2c), top) for N2c, _, _ in cases]   # ζ(t), capped at the top [m]
fig = plt.figure(figsize=(7.6, 5.0))                                            # three columns above, ζ(t) below
gs = fig.add_gridspec(2, 3, height_ratios=[1.3, 1])   # 2 rows × 3 columns layout
dots_p = []   # the three parcel markers
for j, (N2c, lab, col) in enumerate(cases):                                     # draw each column once
    axc = fig.add_subplot(gs[0, j])   # column j in the top row
    shade = np.linspace(1, 0, 50)[:, None] if N2c > 0 else (np.zeros((50, 1)) + 0.5 if N2c == 0 else np.linspace(0, 1, 50)[:, None])   # density shading: stable, neutral, unstable
    axc.imshow(shade, extent=(-1, 1, -100, top), aspect="auto", cmap="Blues", vmin=-0.5, vmax=1.5, origin="lower")  # dark = dense
    axc.axhline(0, color=COLORS["muted"], lw=0.8, ls=":")                       # rest height z_o
    (dp_,) = axc.plot([0], [20.0], "o", color=col, ms=14)                       # the parcel
    dots_p.append(dp_)   # keep it for the update function
    axc.set_xticks([]); axc.set_ylim(-100, top); axc.set_title(lab, fontsize=9)   # fixed height axis and label
    axc.set_ylabel("ζ [m]" if j == 0 else "")   # y label on the first column only
axz = fig.add_subplot(gs[1, :])                                                 # the ζ(t) traces
lines_p = []   # the three growing traces
for (N2c, lab, col), trk in zip(cases, tracks):   # one trace per column
    axz.plot(t_an, trk, color=col, alpha=0.2, lw=1)                             # ghost of the whole solution
    (ln,) = axz.plot([], [], color=col, lw=2.5)                                 # the part played so far
    lines_p.append(ln)   # keep it for the update function
axz.set_xlim(0, 600); axz.set_ylim(-30, top + 10)   # fixed axes
axz.set_xlabel("time t [s]"); axz.set_ylabel("displacement ζ [m]")   # axis labels with units
def update(i):                                                                  # frame i: time t_an[i]
    for dp_, ln, trk in zip(dots_p, lines_p, tracks):   # each column in turn
        dp_.set_data([0], [trk[i]])                                             # move each parcel
        ln.set_data(t_an[:i + 1], trk[:i + 1])                                  # extend each trace
    return dots_p + lines_p   # the artists that changed
show_animation(animate(update, frames=nfr_p, fig=fig, interval=80), player="video")   # smooth MP4
""", explain="""
1. `parcel_displacement` gives ζ(t) for the three signs of N², each released 20 m above its rest height; the runaway
   parcel is stopped at the drawn column top with `np.minimum`.
2. Each column is shaded by density (dark = dense): heavy fluid below light fluid in the stable column, uniform in the
   neutral one, heavy over light in the unstable one.
3. The bottom panel draws each trace as time advances, over a faint ghost of the full solution.
""")
nb.md("""
**What you see.** One parcel bobbing up and down (period 10.5 min), one resting where it was released, one accelerating
away until it hits the top of the drawing.

**How to read it.** The three curves are a cosine, a flat line and a cosh; nothing but the sign of N² differs between the
columns.

**What would change if…** N² were doubled: the period would shrink by √2, from 10.5 to 7.4 minutes, and the runaway parcel
would escape faster by the same factor.

**What would change if…** the parcel felt friction or exchanged heat with its surroundings? The oscillation would die out —
the model has neither, so it rings forever. **Next:** everything hinges on the number N² — let us read it (C51).
""")
# ---------------------------------------------------------------------------------------------------- C51
core("C51", "The Brunt–Väisälä frequency N²: one number that classifies a column",
     "What do the sign and the size of N² tell you about a fluid column before you push anything?")
nb.md("""
#### The problem in plain words
Internal waves in the ocean, lee waves behind mountains, the Richardson number of Ch. 11, the convection switch in climate
models — all read N². It is the natural frequency of a stratified fluid.
""")
nb.md("""
#### The idea
| N² | motion after release | time scale |
|---|---|---|
| > 0 | oscillation, cos Nt | period 2π/N |
| = 0 | stays (or drifts at its initial speed) | none |
| < 0 | runaway, cosh √(−N²) t | e-folding time 1/√(−N²) |
""")
nb.md("""
> 🔁 An imaginary frequency means growth: cos(iσt) = cosh σt, and cosh x = ½(eˣ + e⁻ˣ) (primer P45 in the C50 block).
""")
nb.md(r"""
#### The maths, step by step
1. From D18:
   $$N^2 = -\frac{g}{\rho(z_o)}\left(\frac{d\rho}{dz} - \frac{d\rho_a}{dz}\right) \qquad (1.29) \quad [\mathrm{1/s^2}]$$
   with g [m/s²], ρ [kg/m³] and the density gradients in kg/m⁴.
2. N² > 0 ⇔ dρ/dz < dρ_a/dz: the environment's density falls with height **faster** than a displaced parcel's does.
3. N = √N² [rad/s], and the period is T_N = 2π/N [s].
4. N² < 0: growth rate σ = √(−N²) [1/s] and e-folding time 1/σ. After 3 e-foldings the displacement has grown by
   cosh 3 ≈ ½e³ ≈ 10× — the ½ in cosh halves what a pure exponential (e³ ≈ 20×) would give.
""")
note("C52", "Stable, neutral, unstable", """
**Stable** (N² > 0): a displaced parcel returns and oscillates. **Unstable** (N² < 0): it accelerates away. **Neutral**
(N² = 0): it stays where it is put — reached either when density does not vary at all (dρ/dz = dρ_a/dz = 0) or when the
environment's density falls exactly as fast as a parcel's (dρ/dz = dρ_a/dz), an atmosphere of uniform entropy.
""")
nb.worked_example("reading three values of N²", """
- N² = 10⁻⁴ s⁻² → N = 0.01 s⁻¹ → period 628 s = 10.5 min.
- N² = 4×10⁻⁴ s⁻² → N = 0.02 s⁻¹ → period 314 s.
- N² = −10⁻⁴ s⁻² → σ = 0.01 s⁻¹ → one e-folding every 100 s; after 5 minutes (3 e-foldings) the parcel is cosh(3) = 10.1
  times as far from its rest height as it started.
""")
P("P46", "np.where and np.select", """
`np.where(cond, a, b)` picks a where the condition holds and b elsewhere, element by element; `np.select([c1, c2], [a1, a2],
default)` handles several regimes — that is how labels like "stable" are assigned to whole arrays.
""", code="""
x = np.array([1e-4, 0.0, -1e-4])                                          # three values of N² [1/s²]
print(np.where(x > 0, "stable", "not stable"))                           # two regimes
print(np.select([x > 0, x < 0], ["stable", "unstable"], "neutral"))      # three regimes
""")
nb.code("""
N2_set = np.array([4e-4, 1e-4, 0.0, -1e-4])                   # four columns [1/s²]
print(ch01.classify_stability(N2_set))                        # C52 labels → ['stable' 'stable' 'neutral' 'unstable']
for x in N2_set:                                              # the time scale of each
    kind, secs = ch01.stability_timescale(x)                  # ('period', s), ('none', inf) or ('efold', s)
    print(f"N² = {x:+.0e} 1/s²: {kind} {secs:.1f} s")   # kind and time scale
""", explain="""
1. `classify_stability` labels each N² (with `np.select`-style logic, P46).
2. `stability_timescale` returns the period 2π/N for stable columns, the e-folding time 1/√(−N²) for unstable ones and
   no time scale for a neutral column.
""")
nb.figure("""
t_f = np.linspace(0, 900, 901)                                             # 15 minutes [s]
fig, ax = plt.subplots(figsize=(7.2, 4))                                   # one panel
for N2c, col, ls in ((4e-4, COLORS["teal"], "-"), (1e-4, COLORS["teal"], "--"), (0.0, COLORS["amber"], "-"),   # five columns:
                     (-2.5e-5, COLORS["rose"], "--"), (-1e-4, COLORS["rose"], "-")):                            # N², colour, line style
    ax.plot(t_f, np.clip(ch01.parcel_displacement(t_f, 1.0, N2c), -5, 5), color=col, ls=ls, label=f"N² = {N2c:+.1e} s⁻²")  # ζ/ζ0
    kind, secs = ch01.stability_timescale(N2c)                             # its time scale
    if kind == "period":                                                   # oscillating column
        ax.axvline(secs, color=col, ls=ls, lw=0.8)                         # a dashed vertical at each period
    elif kind == "efold":                                                  # runaway column
        ax.plot([secs], [np.cosh(1.0)], "o", color=col)                    # a dot one e-folding time after release
ax.set_ylim(-5.2, 5.2)                                                     # fixed range; runaways leave the frame
ax.set_xlabel("time t [s]"); ax.set_ylabel("ζ / ζ₀ [–] (clipped at ±5)")   # axis labels
ax.set_title("The sign of N² picks the shape, its size the time scale")   # the message
ax.legend(loc="lower left", fontsize=8, ncol=2)                            # legend in two columns
plt.show()                                                                 # display
""", see="""
Two cosines of different periods (teal), a flat line (amber) and two curves that leave the frame upward (rose). Vertical
lines mark the two periods (314 s and 628 s); dots mark one e-folding time on the growing curves (100 s and 200 s).
""", read="""
The sign of N² decides the shape of the motion; its size sets the time scale — the period 2π/N or the e-folding time
1/√(−N²).
""", change="""
Halve every |N²| and every time scale grows by √2: the periods move right to 444 s and 889 s, and the runaways take
141 s and 283 s per e-folding.
""")
P("P47", "live widgets", """
`live(fn, a=(lo, hi, step))` builds ipywidgets sliders that re-run `fn` each time you release a slider — but only while
Python is running (Colab or a local Jupyter). On the published web page the cell is frozen and replaced by a note; the
static figure above carries the idea there.
""")
nb.live("""
from fluidpy.core.interact import live                                             # ipywidgets sliders (kernel only)
def parcel_plot(N2=1e-4, zeta0=20.0):                                              # called again whenever a slider moves
    t_l = np.linspace(0, 1800, 400)                                                # half an hour [s]
    rho0 = 1025.0                                                                  # seawater density [kg/m³]
    lin = ch01.parcel_displacement(t_l, zeta0, N2)                                 # linear solution (D18 steps 13–14)
    drho_a = -0.005                                                                # the parcel's own density gradient [kg/m⁴]
    drho_env = drho_a - N2 * rho0 / ch01.G0                                        # environment gradient that gives this N² (1.29)
    nl = ch01.parcel_ode_from_gradients(t_l, zeta0, rho0, drho_env, drho_a, zeta_max=2000.0)   # unlinearised, stopped at 2 km
    fig, ax = plt.subplots(figsize=(6.5, 3.2))                                     # one panel
    ax.plot(t_l, lin, color=COLORS["accent"], label="linear ζ₀ cos/cosh")   # linear solution
    ax.plot(t_l, nl, "--", color=COLORS["orange"], label="unlinearised (stops at 2 km)")   # unlinearised solution
    ax.set_ylim(-1.2 * zeta0, 2000); ax.set_xlabel("t [s]"); ax.set_ylabel("ζ [m]"); ax.legend(fontsize=8)   # axes and legend
    plt.show()                                                                     # draw the figure for these slider values
_ = live(parcel_plot, N2=(-2e-4, 4e-4, 1e-5), zeta0=(1.0, 200.0, 1.0))            # two sliders (_ hides the return value)
""", explain="""
1. The callback draws the linear solution and the unlinearised Newton law for the chosen N² and release height. The
   parcel has its own density gradient dρ_a/dz = −0.005 kg m⁻⁴, and the environment gradient is chosen so that (1.29)
   gives the slider's N²; the unlinearised law therefore keeps the parcel's changing density in its denominator.
2. On this 2 km axis the two curves lie on top of each other: as the C50 printout showed, their gap stays a small fraction
   of the displacement (about 10⁻³ at a few hundred metres). The only visible difference is for N² < 0, where the runaway
   unlinearised run is stopped at 2 km while the linear cosh keeps growing.
""")
note("C60", "Seawater", r"""
Seawater density depends on temperature, pressure **and salinity** S (grams of salt per kilogram of seawater, about 35 g/kg
on average). The **potential density** ρ_θ is the density a parcel would have if it were brought isentropically (no heat,
no friction) to a fixed reference pressure — defined properly in C58, below in this section. A displaced parcel keeps its
salt, so for seawater ρ_θ is taken at constant S, and oceanographers refer densities to sea-level pressure. Our linear model: 10 °C and 35 g/kg → 1027.0 kg/m³; 20 °C → 1025.3 kg/m³.
""", equation=r"\rho = \rho(T, p, S)")
note("N21", "How fast a sinking ocean parcel gets denser", r"""
A parcel moved vertically in the ocean is squeezed or relaxed by the pressure of its surroundings, at constant entropy and
salinity, so (stated) its density changes at a rate set by the speed of sound: with c = 1500 m/s and ρ = 1025 kg/m³,
dρ_a/dz = −4.47×10⁻³ kg m⁻⁴.
""", equation=r"\frac{d\rho_a}{dz} \cong -\frac{\rho g}{c^2}")
note("C61", "The ocean's stability criterion", r"""
With N21, (1.29) becomes the ocean's stability test: the column is stable when the compressibility-corrected density
gradient below is negative. **We read (1.35) as "has the same sign as" dρ_θ/dz** — the two sides agree only up to a positive
factor near the reference pressure. Thermocline with dρ/dz = −0.01 kg m⁻⁴: −0.01 + 0.00447 = −5.53×10⁻³ kg m⁻⁴ < 0 → stable,
but N² drops from 9.57×10⁻⁵ to 5.29×10⁻⁵ s⁻² (the period grows from 10.7 to 14.4 min), because the sinking parcel is also
squeezed denser.
""", equation=r"\frac{d\rho_\theta}{dz} \;\sim\; \frac{d\rho}{dz} + \frac{\rho g}{c^2}", ref="1.35")
nb.code("""
print("C60 linear seawater ρ at 10 and 20 °C, 35 g/kg:", ch01.seawater_density_linear(np.array([283.15, 293.15]), 35.0))   # warmer water is lighter
drho_a = ch01.isentropic_density_gradient(1025.0, 1500.0)                         # N21: −ρ g / c² [kg/m⁴]
print(f"N21 parcel density gradient: {drho_a:.3e} kg/m⁴")   # the parcel's own gradient
print(f"C61 corrected gradient: {ch01.ocean_potential_density_gradient(-0.01, 1025.0, 1500.0):.3e} kg/m⁴")   # (1.35)
N2_c = ch01.brunt_vaisala_sq(1025.0, -0.01, drho_a)                              # (1.29) with a compressible parcel
print(f"N² with a compressible parcel: {N2_c:.3e} 1/s², period {ch01.stability_timescale(N2_c)[1]/60:.1f} min")   # weaker stability, longer period
""", explain="""
1. `seawater_density_linear` is a linear equation of state in T and S around 10 °C and 35 g/kg.
2. `isentropic_density_gradient` is −ρg/c² (N21); `ocean_potential_density_gradient` adds it to the observed gradient (C61).
3. Putting the parcel's own gradient into `brunt_vaisala_sq` shows the stability is weaker than the incompressible estimate.
""")
nb.md("""
**What would change if…** the fluid were air, whose parcels expand a lot as they rise? Then dρ_a/dz is large, and it is
easier to compare temperatures than densities. **Next:** the adiabatic lapse rate (C54).
""")
# ---------------------------------------------------------------------------------------------------- C54
core("C54", "The adiabatic lapse rate — and the two ways its sign is written",
     "Air gets colder as you go up. Does that make the atmosphere unstable?")
nb.md("""
#### The problem in plain words
Mountain tops are cold, yet the air over them is usually stable; smog sits under inversions; thunderstorms build on hot
afternoons. In climate science the **lapse-rate feedback** — how the vertical temperature profile changes as the planet
warms — is one of the main feedbacks, and those papers write the lapse rate as a *positive* number for cooling with height.
You need both sign habits.
""")
nb.md("""
#### The idea
```
z ↑      environment T(z): cools at 6.5 K/km
  |   \\    a parcel lifted without heating: cools at 9.8 K/km (it expands as the pressure drops)
  |    \\ \\
  |     \\  \\          at the same height the parcel is COLDER than its surroundings
  |______\\___\\___ T   → denser → sinks back → STABLE, although T falls with height
```
**Stability compares two cooling rates: the environment's and a rising parcel's own.**
""")
note("C49", "A static atmosphere from one profile", r"""
In a fluid at rest, hydrostatics (1.8) and the equation of state (1.12) link p, ρ and T, so one profile fixes the other two:
give T(z), integrate dp/dz = −pg/(RT(z)), then ρ = p/(RT). For T(z) = 288.15 K − 6.5 K/km × z, 5 km up: p = 54.0 kPa and
ρ = 0.736 kg/m³ — the standard atmosphere's values (code below).
""")
note("C53", "The lapse rate in two sign conventions", r"""
The **lapse rate** is the environment's vertical temperature gradient. **Two sign conventions are in daily use; both are
shown everywhere in this section:**

| Convention | Definition | Standard troposphere | Dry adiabatic Γ_a | Stable when |
|---|---|---|---|---|
| Kundu (this book, fluidpy code) | Γ ≡ dT/dz | −6.5 K/km | −9.8 K/km | Γ > Γ_a, e.g. −6.5 > −9.8 |
| Meteorology (climate literature) | Γ ≡ −dT/dz | +6.5 K/km | +9.8 K/km | Γ < Γ_a, e.g. 6.5 < 9.8 |

**Converting: negate the number and flip the inequality** — Γ_met = −Γ_Kundu. The physics never changes.
""")
P("P48", "inequalities under a sign change", """
Multiplying both sides of an inequality by a negative number reverses it: a > b ⇔ −a < −b. That is exactly the move between
the two lapse-rate conventions.
""", code="""
a, b = -6.5, -9.8                  # Kundu lapse rates [K/km]: environment and dry adiabat
print(a > b, -a < -b)              # True True: "−6.5 > −9.8" and "6.5 < 9.8" say the same thing
""")
nb.code("""
Ga = ch01.adiabatic_lapse_rate()                                   # Eq. (1.30), perfect gas: Γa = −g/C_p (Kundu) [K/m]
Ga_met = ch01.lapse_rate_convention(Ga, "meteorology")             # the same number in the meteorology convention [K/m]
print(f"Γa = {Ga*1e3:.3f} K/km (Kundu)  =  {Ga_met*1e3:+.3f} K/km (meteorology)")   # same physics, two signs
for dTdz in (-6.5e-3, -12.0e-3, 0.0, Ga):                          # environment gradients dT/dz [K/m]
    k = ch01.lapse_rate_stability(dTdz)                             # verdict + inequality in Kundu's form
    m = ch01.lapse_rate_stability(dTdz, convention="meteorology")   # the same test written the meteorology way
    print(f"{k.verdict:9s} | Kundu: {k.text:38s} | meteorology: {m.text}")   # verdict and both inequalities
    assert k.code == m.code                                         # the verdict never depends on the convention
""", explain="""
1. `adiabatic_lapse_rate()` returns Γa = −g/C_p in Kundu's convention (negative for air).
2. `lapse_rate_convention` only flips the sign for the meteorology convention.
3. `lapse_rate_stability` computes the verdict once, in Kundu's convention, and writes the criterion with numbers in the
   requested convention; its fields `.verdict`, `.text` and `.code` (+1 stable, 0 neutral, −1 unstable) are printed.
4. The assert proves that the verdict is the same in both conventions — only the bookkeeping differs.
""")
P("P49", "chain rule", r"""
If f depends on x and y and both change, $df = (\partial f/\partial x)_y\,dx + (\partial f/\partial y)_x\,dy$. Divide by dz
to get the rate along a path.
""", code="""
f2 = lambda x, y: x**2 * y                       # a function of two variables
dx_, dy_ = 1e-6, 2e-6                            # small changes of x and y at (1, 2)
print(f2(1 + dx_, 2 + dy_) - f2(1, 2))           # the actual change ≈ 6.0e-06
print(2 * 1 * 2 * dx_ + 1**2 * dy_)              # (∂f/∂x) dx + (∂f/∂y) dy = 4e-6 + 2e-6 = 6.0e-06
""")
P("P50", "Gibbs free energy", r"""
g ≡ h − Ts [J/kg]. We use it only as a device: combined with the Gibbs relation dh = T ds + v dp, its differential is
simply dg = v dp − s dT — a state function whose changes involve only dp and dT.
""", code="""
T_, s_, v_, ds_, dp_, dT_ = sp.symbols("T s v ds dp dT")      # symbols for a small change
dh_ = T_ * ds_ + v_ * dp_                                     # Gibbs: dh = T ds + v dp
dg_ = dh_ - (T_ * ds_ + s_ * dT_)                             # dg = dh − d(Ts) with the product rule (P38)
print(sp.expand(dg_))                                         # → dp·v − dT·s: the T ds terms cancel
""")
P("P51", "exact differentials and Maxwell relations", r"""
If dg = A dp + B dT for a smooth state function g, then A = (∂g/∂p)_T and B = (∂g/∂T)_p, and because mixed second
derivatives are equal, (∂A/∂T)_p = (∂B/∂p)_T. Applied to dg = v dp − s dT this gives a **Maxwell relation**,
(∂v/∂T)_p = −(∂s/∂p)_T: an entropy derivative that cannot be measured directly equals one that can.
""", code="""
xs_, ys_ = sp.symbols("x y")                                  # two variables
f_mix = xs_**2 * sp.sin(ys_)                                  # a smooth function
print(sp.diff(f_mix, xs_, ys_) == sp.diff(f_mix, ys_, xs_))   # ∂²f/∂x∂y equals ∂²f/∂y∂x → True
""")
D("D19", "The adiabatic lapse rate without perfect-gas relations", ref="1.30",
  goal="""Find the rate at which the temperature of a parcel changes as it moves isentropically up through a fluid in
hydrostatic balance — for *any* fluid, not just a perfect gas — and then read it in both lapse-rate sign conventions.""",
  assumptions="""Isentropic parcel, ds = 0 (step 1) · a single-component fluid whose state is fixed by T and p (step 2) ·
smooth state functions, so mixed partials commute (step 7) · parcel pressure = environment pressure, and the environment
is hydrostatic (step 14) · evaluated where the parcel still has the environment's density (step 14) · perfect gas only in
the last step (step 16).""",
  start=(r"T\,ds = dh - v\,dp \ \ (1.18) \text{ with } ds = 0, \qquad \dfrac{dp}{dz} = -\rho g \ \ (1.8)",
         "the parcel's entropy stays constant, and the pressure it feels is the hydrostatic pressure of its surroundings."),
  plan=["Write dh for the isentropic parcel.",
        "Expand dh in T and p, which brings in C_p and an unknown (∂h/∂p)_T.",
        "Find (∂h/∂p)_T with the Gibbs free energy, a Maxwell relation and α.",
        "Divide by dz and use hydrostatics."],
  uses=["Gibbs relations (C35, D10)", "C_p (1.14) (C30)", "α (1.20) (C37)", "two properties fix the state (C28)",
        "hydrostatic law (C20, D05)", "chain rule (P49)", "partial derivative with a variable held fixed (P39)",
        "Gibbs free energy (P50)", "exact differentials and Maxwell relations (P51)",
        "inequalities under a sign change (P48, for the result)"],
  steps=[
      ("Set ds = 0 in the second Gibbs form", r"dh = v\,dp",
       "The parcel moves without heat or friction, so its entropy is constant (T ds = 0). We use the enthalpy form because "
       "it contains dp, which hydrostatics will supply.",
       "The parcel's enthalpy changes only through the pressure it meets."),
      ("Expand dh in T and p", r"dh = \Big(\dfrac{\partial h}{\partial T}\Big)_p dT + \Big(\dfrac{\partial h}{\partial p}\Big)_T dp",
       "Two properties fix the state (C28), so h = h(T, p); the chain rule for a function of two variables (P49) gives its "
       "differential. We want dT to appear.",
       "h changes because T changes and because p changes."),
      ("Recognise C_p", r"dh = C_p\,dT + \Big(\dfrac{\partial h}{\partial p}\Big)_T dp",
       "(1.14) defines C_p = (∂h/∂T)_p. One coefficient is now a known property; the other still needs work.",
       "The first part is the specific heat times the temperature change."),
      ("Get (∂h/∂p)_T from Gibbs at fixed T", r"\Big(\dfrac{\partial h}{\partial p}\Big)_T = T\Big(\dfrac{\partial s}{\partial p}\Big)_T + v",
       "Gibbs dh = T ds + v dp holds for any change (D10); take a change at constant T and divide by dp (P39). This trades "
       "the unknown for an entropy derivative.",
       "How h depends on p at fixed T involves how s does."),
      ("Introduce the Gibbs free energy", r"dg = v\,dp - s\,dT,\qquad g \equiv h - Ts",
       "dg = dh − T ds − s dT (product rule on Ts, P38) and dh = T ds + v dp; the T ds terms cancel (P50). We choose g because "
       "its natural variables are exactly p and T.",
       "A helper state function whose changes involve only dp and dT."),
      ("Read off the partial derivatives of g", r"\Big(\dfrac{\partial g}{\partial p}\Big)_T = v,\qquad \Big(\dfrac{\partial g}{\partial T}\Big)_p = -s",
       "Compare dg = v dp − s dT with the chain rule dg = (∂g/∂p)_T dp + (∂g/∂T)_p dT term by term (P51).",
       "v and −s are the slopes of g."),
      ("Equate the mixed second derivatives", r"\Big(\dfrac{\partial v}{\partial T}\Big)_p = -\Big(\dfrac{\partial s}{\partial p}\Big)_T",
       "g is a smooth state function, so ∂²g/∂T∂p = ∂²g/∂p∂T (equality of mixed partials, P51); differentiating step 6's "
       "first result in T and the second in p gives this **Maxwell relation**. It converts an entropy derivative into a "
       "measurable one.",
       "How entropy changes with pressure is set by thermal expansion."),
      ("Express (∂v/∂T)_p with α", r"\Big(\dfrac{\partial v}{\partial T}\Big)_p = v\,\alpha",
       "v = 1/ρ gives (∂v/∂T)_p = −(1/ρ²)(∂ρ/∂T)_p; by (1.20) (∂ρ/∂T)_p = −ρα, so the result is α/ρ = vα.",
       "Heating at fixed pressure swells unit mass by vα per kelvin."),
      ("Insert into the Maxwell relation", r"\Big(\dfrac{\partial s}{\partial p}\Big)_T = -\,v\,\alpha",
       "Steps 7 and 8 together. The sign matters: squeezing at fixed T lowers the entropy when α > 0.",
       "Entropy falls with pressure at the rate vα."),
      ("Complete (∂h/∂p)_T", r"\Big(\dfrac{\partial h}{\partial p}\Big)_T = v\,(1 - \alpha T)",
       "Substitute step 9 into step 4: T(−vα) + v = v(1 − αT). The unknown coefficient is now written with measurable "
       "properties only.",
       "The pressure part of dh, in measurable properties."),
      ("Equate the two expressions for dh", r"C_p\,dT + v(1 - \alpha T)\,dp = v\,dp",
       "Step 3 (with step 10) and step 1 describe the same change of the same parcel, so their right-hand sides are equal.",
       "Two bookkeepings of one enthalpy change."),
      ("Subtract v dp from both sides", r"C_p\,dT = v\,\alpha\,T\,dp",
       "v(1 − αT)dp − v dp = −vαT dp; moving it across makes it positive. Temperature change is now tied to pressure "
       "change alone.",
       "A parcel warms when compressed and cools when the pressure drops (for α > 0)."),
      ("Divide by dz along the parcel's path", r"C_p\,\dfrac{dT_a}{dz} = v\,\alpha\,T\,\dfrac{dp}{dz}",
       "Both sides are changes along the vertical path; dividing by dz turns them into rates per metre. The subscript a "
       "marks the adiabatic parcel.",
       "Temperature rate per metre ↔ pressure rate per metre."),
      ("Use hydrostatics with v = 1/ρ", r"C_p\,\dfrac{dT_a}{dz} = -\,\alpha\,T\,g",
       "The parcel's pressure is the environment's, dp/dz = −ρg (1.8); at the release height the parcel's v equals the "
       "environment's 1/ρ, so vρ = 1.",
       "The pressure drop with height sets the cooling."),
      ("Divide by C_p", r"\dfrac{dT_a}{dz} \equiv \Gamma_a = -\,\dfrac{g\,\alpha\,T}{C_p}",
       "C_p > 0, so dividing by it is allowed and changes no sign. This is (1.30), derived without any perfect-gas relation.",
       "The adiabatic lapse rate of any fluid."),
      ("Specialise to a perfect gas", r"\Gamma_a = -\,\dfrac{g}{C_p} = -9.76\ \text{K/km}",
       "α = 1/T for a perfect gas (1.28, C48), so αT = 1; with g = 9.807 m/s² and C_p = 1004.7 J kg⁻¹ K⁻¹ the number follows.",
       "Dry air cools 9.76 K for every kilometre it rises."),
  ],
  result=(r"\Gamma_a = \dfrac{dT_a}{dz} = -\dfrac{g\,\alpha\,T}{C_p}; \qquad \text{dry air: } \Gamma_a = -g/C_p = -9.76 \text{ K/km}",
          "**Kundu (Γ ≡ dT/dz):** Γ_a = −9.76 K/km, and the environment is stable when **Γ > Γ_a** (e.g. −6.5 > −9.76). "
          "**Meteorology (Γ ≡ −dT/dz):** Γ_a,met = +gαT/C_p = +9.76 K/km, and it is stable when **Γ < Γ_a** (6.5 < 9.76). "
          "Converting multiplies both sides by −1, which reverses the inequality (P48); the verdict is identical."),
  interpret="""Γ_a is the slope a well-mixed (isentropic) layer settles to — the reference against which every observed
profile is judged. A column cooling more slowly than 9.76 K/km (dT/dz > −9.76 K/km, or Γ_met < 9.76 K/km) is stable even
though it cools with height. It fails for moist air (latent heat makes the saturated rate about 4–7 K/km), for large
excursions where the parcel's temperature departs from the environment's (then the local rate is −(g/C_p)T_parcel/T_env),
and for mixtures whose composition changes (seawater salinity). In the explainer `parcel_stability` (C55 block), the dashed
adiabat's slope is exactly this number.""",
  check="""Units: m s⁻² · K⁻¹ · K / (J kg⁻¹ K⁻¹) = K/m ✓. Limits: α = 0 (a fluid that does not expand) → Γ_a = 0, no
adiabatic temperature change ✓; perfect gas → −g/C_p ✓. Numbers: air −9.76 K/km; water at 283 K (α = 1.5×10⁻⁴ K⁻¹,
C_p = 4190 J kg⁻¹ K⁻¹) → −0.10 K/km ✓. The sympy cell below verifies steps 7–16 for an arbitrary Gibbs free energy, and
the from-scratch cell differentiates a lifted parcel's temperature and finds −g/C_p at the release height ✓.""",
  traps="the sign of the Maxwell relation; using dh = C_p dT (true only for a perfect gas); using the parcel's v with the "
        "environment's ρ away from the release point; **negating Γ without flipping the inequality** when switching to "
        "the meteorology convention.",
  check_src="""
import sympy as sp                                             # symbolic algebra (primer P40)
T_g, p_g, g_g, R_g, cp_g = sp.symbols("T p g R c_p", positive=True)   # temperature, pressure, gravity, gas constant, C_p
G = sp.Function("G")(T_g, p_g)                                 # ANY smooth Gibbs free energy per unit mass g(T, p)
v_g = sp.diff(G, p_g)                                          # step 6: v = (∂g/∂p)_T
s_g = -sp.diff(G, T_g)                                         # step 6: s = −(∂g/∂T)_p
h_g = G + T_g * s_g                                            # g = h − Ts rearranged: h = g + Ts
Cp_g = sp.diff(h_g, T_g)                                       # (1.14): C_p = (∂h/∂T)_p
alpha_g = sp.diff(v_g, T_g) / v_g                              # step 8: α = (∂v/∂T)_p / v
# Slope of an isentrope. Along it s stays constant, so by the chain rule (P49)
#   ds = (∂s/∂T)_p dT + (∂s/∂p)_T dp = 0   ⇒   dT/dp = −(∂s/∂p)_T / (∂s/∂T)_p   (the "implicit-function rule")
dTdp_s = -sp.diff(s_g, p_g) / sp.diff(s_g, T_g)                # that ratio, with s written through the arbitrary G
print(sp.simplify(dTdp_s - v_g * alpha_g * T_g / Cp_g))        # 0 → step 12 (C_p dT = v α T dp) holds for any fluid
print(sp.simplify(dTdp_s * (-g_g / v_g) + g_g * alpha_g * T_g / Cp_g))   # 0 → (1.30) after dp/dz = −g/v (steps 13–15)
G_perfect = cp_g * (T_g - T_g * sp.log(T_g)) + R_g * T_g * sp.log(p_g)  # a perfect gas with constant C_p
# .subs(G, G_perfect) replaces the unknown function G by the concrete perfect-gas formula; the derivatives of G are
# then still written as unevaluated "Derivative(...)" objects, and .doit() carries them out on the concrete formula.
print(sp.simplify((-g_g * alpha_g * T_g / Cp_g).subs(G, G_perfect).doit()))   # −g/c_p → step 16
""")
note("N17", "The number, in both conventions", r"""
Earth's dry air: Γ_a = −g/C_p = −9.81/1004.7 = **−9.76 K/km in Kundu's convention = +9.76 K/km in the meteorology
convention** (usually quoted as 9.8 K/km). Water is hardly affected by adiabatic compression: α = 1.5×10⁻⁴ K⁻¹ and
C_p = 4190 J kg⁻¹ K⁻¹ at 283 K give Γ_a = −0.10 K/km.
""")
nb.worked_example("is the standard troposphere stable?", """
Environment dT/dz = −6.5 K/km; Γ_a = −9.76 K/km.
- **Kundu:** −6.5 > −9.76 ✓ stable. **Meteorology:** Γ = +6.5, Γ_a = +9.76, and 6.5 < 9.76 ✓ stable. The margin is 3.26 K/km
  in both.
- Lift a parcel 1 km from 288.15 K: it cools to 278.39 K, while the environment there is 281.65 K — the parcel is 3.26 K
  colder, hence denser, so it sinks back.
- N² = (g/T)(dT/dz − Γ_a) = (9.81/288.15)(3.26×10⁻³) = 1.11×10⁻⁴ s⁻² → period 9.9 min.
- Over hot ground with dT/dz = −12 K/km: −12 < −9.76 (meteorology: 12 > 9.76) → unstable; at 300 K a displacement e-folds
  every 117 s.
""")
nb.code("""
N2_trop = ch01.brunt_vaisala_sq_from_lapse(288.15, -6.5e-3)          # N² = (g/T)(dT/dz + g/C_p) at the ground [1/s²]
print(f"standard troposphere: N² = {N2_trop:.3e} 1/s², {ch01.stability_timescale(N2_trop)}")      # period ≈ 596 s
print(f"parcel 1 km up: {ch01.parcel_temperature(288.15, 1000.0):.2f} K; environment {288.15 - 6.5:.2f} K")   # parcel colder than environment → sinks
N2_hot = ch01.brunt_vaisala_sq_from_lapse(300.0, -12e-3)             # superadiabatic layer over hot ground
print(f"hot ground: N² = {N2_hot:.3e} 1/s², {ch01.stability_timescale(N2_hot)}")                  # e-folds in ≈ 117 s
print(f"water: Γa = {ch01.adiabatic_lapse_rate(T=283.15, cp=4190.0, alpha=1.5e-4)*1e3:.3f} K/km")   # (1.30) general form
z5 = np.array([0.0, 5000.0])                                          # the ground and 5 km [m]
p5, rho5, T5 = ch01.atmosphere_from_temperature(z5, lambda zz: 288.15 - 6.5e-3 * zz)   # C49: p, ρ, T from T(z) alone
print(f"C49 at 5 km: p = {p5[-1]/1e3:.1f} kPa, ρ = {rho5[-1]:.3f} kg/m³; USSA-1976: {ch01.standard_atmosphere(5000.0)[1]/1e3:.1f} kPa")   # matches the standard atmosphere
""", explain="""
1. `brunt_vaisala_sq_from_lapse` is N² in lapse-rate form (derived as D36 in the next block); `stability_timescale` turns it
   into a period or an e-folding time.
2. `parcel_temperature` follows the dry adiabat from the ground, T₀ + Γ_a z.
3. `adiabatic_lapse_rate` with α and C_p given uses the general form (1.30) — for water.
4. `atmosphere_from_temperature` integrates hydrostatics with p = ρRT for a given T(z): one profile fixes the other two (C49).
""")
nb.md("""
#### From scratch: lift a parcel and measure its slope
Take the standard atmosphere's pressure profile, lift a dry parcel from the ground along the isentrope (1.26), and
differentiate its temperature numerically.
""")
nb.check_agree("""
z_lift = np.linspace(0, 200, 201)                                     # heights every metre up to 200 m [m]
T_env, p_env, _ = ch01.standard_atmosphere(z_lift)                    # the environment the parcel rises through (USSA-1976)
T_parcel = 288.15 * (p_env / p_env[0])**((ch01.GAMMA_AIR - 1) / ch01.GAMMA_AIR)   # (1.26): its own temperature at each height
dTdz_parcel = np.gradient(T_parcel, z_lift, edge_order=2)             # its slope dT/dz along the path (P22) [K/m]
assert np.isclose(dTdz_parcel[0], ch01.adiabatic_lapse_rate(), rtol=1e-3)   # = Γa = −g/C_p at the release height (Kundu)
assert np.isclose(-dTdz_parcel[0], ch01.lapse_rate_convention(ch01.adiabatic_lapse_rate(), "meteorology"), rtol=1e-3)  # = +g/C_p
print(f"parcel slope at release: {dTdz_parcel[0]*1e3:.3f} K/km (Kundu) = {-dTdz_parcel[0]*1e3:+.3f} K/km (meteorology)")   # both conventions
""")
nb.md("""
At the release height parcel and environment coincide, so the slope is exactly (1.30). Higher up the parcel is colder than
the air around it and its slope becomes −(g/C_p)(T_parcel/T_env) — (1.30) is the *local* rate where the two still match,
the point step 14 of D19 makes.
""")
nb.figure("""
z_km = np.linspace(0, 3, 61)                                                 # heights [km]
fig, ax = plt.subplots(figsize=(7.6, 4.4))                                   # one panel: temperature across, height up
rows = []                                                                    # text rows for the side table
for dTdz_km, col, name in ((-6.5, COLORS["ink"], "standard"), (-12.0, COLORS["rose"], "hot ground"), (5.0, COLORS["teal"], "inversion")):  # 3 profiles
    ax.plot(288.15 + dTdz_km * z_km, z_km, color=col, lw=2.5, label=f"{name}: dT/dz = {dTdz_km:+.1f} K/km".replace("-", "−"))  # T(z)
    k = ch01.lapse_rate_stability(dTdz_km * 1e-3, prefix=False)              # Kundu: verdict and bare inequality
    m = ch01.lapse_rate_stability(dTdz_km * 1e-3, convention="meteorology", prefix=False)   # meteorology form
    rows.append(f"{name:10s} {k.verdict:8s}  {k.text:22s}  ⇔  {m.text}")    # one row of the table
ax.plot(ch01.parcel_temperature(288.15, z_km * 1e3), z_km, "--", color=COLORS["accent"], lw=2, label="dry adiabat (parcel)")  # Γa line
ax.text(0.02, 0.03, "Kundu (dT/dz vs Γa)                ⇔  meteorology (Γ vs Γa)\\n" + "\\n".join(rows),   # header + rows
        transform=ax.transAxes, fontsize=7.5, family="monospace", va="bottom",   # placed in axes coordinates, fixed-width font
        bbox=dict(facecolor="white", edgecolor=COLORS["grid"]))              # the convention table inside the axes
ax.set_xlabel("temperature T [K]"); ax.set_ylabel("height z [km]")          # axis labels with units
ax.set_ylim(0, 3); ax.set_xlim(195, 310)                                    # room on the left for the table
ax.set_title("Stable where the environment cools more slowly than the dashed adiabat")   # the message
ax.legend(loc="upper right", fontsize=8)                                      # legend
savefig(fig, "ch01", "c54_lapse_rates")                                       # keep a copy for review
plt.show()                                                                    # display
""", see="""
Three environment profiles starting at 288.15 K — black (cooling 6.5 K/km), rose (cooling 12 K/km) and teal (warming 5 K/km,
an inversion) — and one dashed purple dry adiabat. The box lists each profile's verdict with its inequality written both
ways.
""", read="""
Where an environment line leans *less* than the dashed adiabat (cools more slowly), a lifted parcel ends up colder than its
surroundings → stable; where it leans *more*, the parcel ends up warmer → unstable. In the box, each Kundu inequality and its
meteorology twin always give the same verdict.
""", change="""
Make the environment cool at exactly 9.76 K/km and its line lies on top of the adiabat: neutral — a parcel is at home at
every height, and both inequalities become equalities.
""")
nb.plotly("""
z3 = np.linspace(0, 3e3, 61)                                                   # heights [m]
grads = np.linspace(-15e-3, 10e-3, 26)                                         # environment dT/dz from −15 to +10 K/km [K/m]
def column(dTdz):                                                              # one slider position
    return {"environment T(z)": (288.15 + dTdz * z3, z3 / 1e3),   # environment line [K, km]
            "dry adiabat": (ch01.parcel_temperature(288.15, z3), z3 / 1e3)}   # the parcel's line, fixed
fig = slider_figure(column, "dT/dz", grads * 1e3, unit="K/km", xlabel="temperature T [K]", ylabel="height z [km]",   # 26 slider positions
                    title="", active=int(np.argmin(np.abs(grads + 6.5e-3))))   # start at −6.5 K/km
for i, x in enumerate(grads):                                                  # per-position title with both conventions
    k = ch01.lapse_rate_stability(x, prefix=False)                             # Kundu form
    m = ch01.lapse_rate_stability(x, convention="meteorology", prefix=False)   # meteorology form
    kind, secs = ch01.stability_timescale(ch01.brunt_vaisala_sq_from_lapse(288.15, x))   # period or e-folding time
    extra = {"period": f"period {secs/60:.1f} min", "efold": f"e-folds in {secs:.0f} s", "none": "no restoring force"}[kind]   # time-scale text
    st = fig.layout.sliders[0].steps[i]                                        # the slider step (method "update")
    st.args = [st.args[0], {"title.text": f"dT/dz: {k.text}  ⇔  Γ_met: {m.text}: <b>{k.verdict}</b>, {extra}"}]   # this position's title
fig.layout.title.text = fig.layout.sliders[0].steps[fig.layout.sliders[0].active].args[1]["title.text"]   # starting title
fig.show()                                                                     # drag the lapse rate
""", explain="""
1. For each environment gradient from −15 to +10 K/km the figure shows the environment line and the fixed dry adiabat.
2. `lapse_rate_stability` writes the criterion in Kundu's form and in the meteorology form; the two texts are joined with ⇔
   in each slider position's title (the step's second dictionary updates the layout, as in the C20 slider).
3. `brunt_vaisala_sq_from_lapse` and `stability_timescale` add the period or the e-folding time.
""")
see_read_change(
    'An environment temperature line and the fixed dashed dry adiabat, both starting at 288.15 K at the ground; the title states the criterion in both conventions, the verdict and the time scale.',
    'If the environment line leans less than the adiabat (cools more slowly), a lifted parcel is colder than its surroundings: stable, and the title gives an oscillation period. If it leans more, the title gives an e-folding time.',
    'Put the slider on −10 K/km, just steeper than Γa: the verdict turns unstable, the Kundu inequality reads −10.0 < −9.8 and the meteorology one 10.0 > 9.8 — both flip together.',
)
nb.md("""
**What would change if…** you read a paper saying "the lapse rate increases under warming"? In its convention Γ = −dT/dz
grows, i.e. the air cools *faster* with height and moves *toward* Γ_a — it becomes less stable. In our code that is dT/dz
becoming more negative. **Next:** a temperature label that removes the parcel's own cooling altogether (C55).
""")
# ---------------------------------------------------------------------------------------------------- C55
core("C55", "Potential temperature θ: a temperature label that survives lifting",
     "Can we label a parcel of air with a temperature that does not change when it rises or sinks without heating?")
nb.md("""
#### The problem in plain words
Is air at 500 hPa and −23 °C "colder" than surface air at 15 °C? Bring it down and it warms by compression. Weather maps,
isentropic analysis and all of dry atmospheric dynamics (Ch. 13) use a temperature with that effect removed.
""")
nb.md("""
#### The idea
Bring every parcel adiabatically to a reference pressure p_o and read its thermometer there.

**Which p_o?** The book takes p_o = p(0), the sea-level pressure, close to 100 kPa. Meteorology — and fluidpy's
`potential_temperature` (`P_REF`) — fixes it at exactly 1000 hPa = 100 kPa. The two differ by about 1 % at sea level
(1013 vs 1000 hPa), which is why the code below finds θ ≈ 287.1 K rather than 288.15 K at the ground. What matters is
that p_o is **one fixed constant** for every parcel; fluidpy keeps it (`p_ref`) separate from the statics pressure `p0`,
the pressure at z = 0 of a particular column, which can vary from day to day.

```
500 hPa    T = 250 K  ──(compress adiabatically down to 1000 hPa)──►  θ = 305 K
1000 hPa   T = 288 K  ──(already there)────────────────────────────►  θ = 288 K
θ increases upward  ⇒  a lifted parcel is always cooler (in θ) than the air above it  ⇒  stable
```
The tools: the exponent rules (P43, C06 block) and the isentropic ratio (1.26) (C46, C45 block).
""")
D("D20", "Potential temperature", ref="1.31",
  goal="""Define a temperature label that a parcel keeps when it moves up or down without heating: the temperature it would
have if brought adiabatically to a standard pressure p_o.""",
  assumptions="""Perfect gas with constant γ (start) · the imagined trip to p_o is adiabatic and reversible (step 1) · p_o is
one fixed reference pressure for all parcels (step 1): sea-level pressure ≈ 100 kPa in the book, exactly 1000 hPa in
meteorology and in fluidpy.""",
  start=(r"\dfrac{T}{T_0} = \Big(\dfrac{p}{p_0}\Big)^{(\gamma-1)/\gamma}",
         "along one isentrope of a perfect gas, the temperature ratio of any two states is a power of their pressure ratio (1.26)."),
  plan=["Choose the two states as 'now' and 'at the reference pressure'.", "Solve for the reference temperature.",
        "Rewrite the exponent with R/C_p."],
  uses=["isentropic ratios (C46, C45)", "exponent rules (P43)", "R = C_p − C_v and γ = C_p/C_v (C42, C43)"],
  steps=[
      ("Take the reference state (θ, p_o) and the current state (T, p)", r"\dfrac{T}{\theta} = \Big(\dfrac{p}{p_o}\Big)^{(\gamma-1)/\gamma}",
       "(1.26) links any two states on one isentrope; by definition the parcel taken adiabatically to p_o has the temperature "
       "θ, so rename T₀ → θ and p₀ → p_o.",
       "The parcel now and the parcel at the reference pressure lie on the same isentrope."),
      ("Multiply both sides by θ", r"T = \theta\,\Big(\dfrac{p}{p_o}\Big)^{(\gamma-1)/\gamma}",
       "θ ≠ 0, so multiplying both sides by it keeps the equation true; this is the book's form (1.31).",
       "Actual temperature = potential temperature scaled by a pressure factor."),
      ("Solve for θ", r"\theta = T\,\Big(\dfrac{p_o}{p}\Big)^{(\gamma-1)/\gamma}",
       "Divide both sides by the power and use 1/(a/b)^k = (b/a)^k (P43). This is the form we compute.",
       "Correct the thermometer reading for the pressure."),
      ("Rewrite the exponent with C_v and C_p", r"\dfrac{\gamma-1}{\gamma} = 1 - \dfrac{C_v}{C_p}",
       "Split the fraction: (γ − 1)/γ = 1 − 1/γ, and 1/γ = C_v/C_p by (1.24).",
       "The exponent is one minus the inverse of the specific-heat ratio."),
      ("Use R = C_p − C_v", r"\theta = T\,\Big(\dfrac{p_o}{p}\Big)^{R/C_p}",
       "1 − C_v/C_p = (C_p − C_v)/C_p = R/C_p by (1.23); for air R/C_p = 287.06/1004.7 = 2/7.",
       "The form used in meteorology, with exponent 0.286."),
  ],
  result=(r"T = \theta\,(p/p_o)^{(\gamma-1)/\gamma} \ \ (1.31), \qquad \theta = T\,(p_o/p)^{R/C_p}",
          "potential temperature is the temperature a parcel would have at the reference pressure p_o (1000 hPa in fluidpy)."),
  interpret="""θ is the conserved temperature of dry adiabatic motion: two parcels can be compared regardless of their
heights. Its vertical gradient is the stability test (C57, D36). It changes when heat is added (radiation, condensation) —
which is how diabatic heating is measured.""",
  check="""The exponent is dimensionless and so is p_o/p ✓. At p = p_o, θ = T ✓. Adiabatic motion keeps θ constant, because
it only moves the parcel along its isentrope ✓. Number: 500 hPa, 250 K → θ = 250 × 2^(2/7) = 304.8 K ✓
(`potential_temperature(250.0, 5e4)` and the from-scratch formula below agree).""",
  traps="p_o is a fixed constant — do not use the local pressure p, or a surface pressure that varies from place to "
        "place, in its place (fluidpy keeps `p_ref` separate from the statics `p0` for this reason); inverting the ratio "
        "(p/p_o instead of p_o/p) in the computed form.")
nb.worked_example("two parcels", """
θ = T(p_o/p)^(R/C_p) with R/C_p = 287.06/1004.7 = 0.2857 = 2/7.
- A parcel at 500 hPa and 250 K: θ = 250 × 2^0.2857 = 250 × 1.2190 = 304.8 K.
- Surface air at 1000 hPa and 288 K: θ = 288 K.
- The upper parcel is 38 K colder in T but 17 K warmer in θ — the column is stable.
""")
nb.code("""
print(f"θ(250 K, 500 hPa) = {ch01.potential_temperature(250.0, 5.0e4):.2f} K")        # Eq. (1.31) → 304.75 K
z_tr = np.linspace(0, 11e3, 111)                                                        # the troposphere, every 100 m [m]
T_tr, p_tr, rho_tr = ch01.standard_atmosphere(z_tr)                                     # USSA-1976 (geopotential height)
theta_tr = ch01.potential_temperature(T_tr, p_tr)                                       # θ at every height [K] (p_o = 1000 hPa)
dth = np.gradient(theta_tr, z_tr, edge_order=2)                                          # dθ/dz [K/m] (P22)
N2_theta = ch01.brunt_vaisala_sq_from_theta(theta_tr, dth)                               # N² = (g/θ) dθ/dz (D36 result)
N2_lapse = ch01.brunt_vaisala_sq_from_lapse(T_tr, -6.5e-3)                               # N² = (g/T)(dT/dz + g/C_p)
assert np.allclose(N2_theta, N2_lapse, rtol=1e-3)                                        # two forms of the same N²
print(f"θ rises from {theta_tr[0]:.2f} K at the ground to {theta_tr[-1]:.1f} K at 11 km")   # stable: θ increases upward
print(f"N² at the ground: from θ {N2_theta[0]:.4e}, from the lapse rate {N2_lapse[0]:.4e} 1/s²")   # the two forms agree
print(f"C56 dθ/dz at the ground: {ch01.potential_temperature_gradient(T_tr[0], -6.5e-3, p=p_tr[0])*1e3:.2f} K/km")   # (1.32)
""", explain="""
1. `potential_temperature` is (1.31) with p_o = 1000 hPa — note θ at the ground is about 287.1 K, not 288.15 K, because
   the standard ground pressure 1013 hPa is slightly above the reference pressure (see "Which p_o?" above).
2. θ increases steadily through the standard troposphere: it is stable everywhere.
3. N² from the θ gradient (`brunt_vaisala_sq_from_theta`) and from the lapse rate (`brunt_vaisala_sq_from_lapse`) agree — the
   result of D36 below, checked on real profiles.
4. `potential_temperature_gradient` is (1.32): dθ/dz = (θ/T)(dT/dz + g/C_p).
""")
nb.md("#### From scratch: θ by the formula")
nb.check_agree("""
theta_mine = T_tr * (ch01.P_REF / p_tr)**((ch01.GAMMA_AIR - 1) / ch01.GAMMA_AIR)   # θ = T (p_o/p)^((γ−1)/γ), D20 step 3
assert np.allclose(theta_mine, ch01.potential_temperature(T_tr, p_tr))              # = the library
assert np.allclose(ch01.temperature_from_potential(theta_mine, p_tr), T_tr)         # and back again (D20 step 2)
print(f"θ at 5 km: {theta_mine[50]:.2f} K (mine) = {ch01.potential_temperature(T_tr[50], p_tr[50]):.2f} K (library)")   # one sample value
""")
note("C56", "θ's gradient and the two lapse rates", r"""
Taking logarithms of (1.31), differentiating in z, and using dp/dz = −ρg, p = ρRT and α = 1/T gives (stated) the link
between θ's gradient and the two lapse rates, written here in both conventions. Standard troposphere at the ground:
(θ/T)(−6.5 + 9.76) K/km = 3.25 K/km.
""", equation=r"\frac{T}{\theta}\frac{d\theta}{dz} = \frac{dT}{dz} + \frac{g}{C_p} = \Gamma - \Gamma_a\ \text{(Kundu)} = \Gamma_{a,\rm met} - \Gamma_{\rm met}\ \text{(meteorology)}",
     ref="1.32")
note("N19", "The log-derivative of (1.31)", """
The logarithmic derivative of (1.31) is step 5 of the next derivation; the same move gives the atmosphere's θ profiles in
Ch. 13 §13.2.
""")
P("P52", "logarithmic differentiation", r"""
d(ln f)/dz = (1/f) df/dz, and logarithms turn products and powers into sums: ln(ab^k) = ln a + k ln b. So the relative rate
of change of a product of powers is the weighted sum of the relative rates of its factors.
""", code="""
zq, hq = 2.0, 1e-6                                                 # a point and a small step
print((np.log((zq + hq)**3) - np.log(zq**3)) / hq, 3 / zq)          # d ln(z³)/dz = 3/z = 1.5 at z = 2
""")
D("D36", "N² from the potential-temperature gradient (links Eqs. 1.29 and 1.32; our addition)", ref="",
  goal="""Show that for a perfect-gas atmosphere the density-based N² of (1.29) is simply g/θ times the vertical gradient of
θ — so stability means "θ increases upward" — and connect it to the two lapse rates.""",
  assumptions="""Perfect gas with constant γ (steps 2–3) · the parcel's pressure equals the environment's (step 3) · at the
parcel's rest height ρ_a = ρ (step 1).""",
  start=(r"N^2 = -\dfrac{g}{\rho}\Big(\dfrac{d\rho}{dz} - \dfrac{d\rho_a}{dz}\Big)",
         "the parcel result (1.29) of D18, evaluated at the height where parcel and environment share the density ρ."),
  plan=["Write both density gradients as relative (logarithmic) rates.",
        "Environment from p = ρRT; parcel from the isentropic ratio.",
        "Subtract: the pressure terms combine into θ.", "Use (1.32) for the lapse-rate form."],
  uses=["(1.29) (C51, D18)", "logarithmic differentiation (P52)", "p = ρRT (C40)", "isentropic ratios (1.26) (C46)",
        "θ (1.31) (D20)", "(1.32) (C56)", "γ and R/C_p (C42, C43)"],
  steps=[
      ("Take ρ inside the bracket", r"N^2 = -g\Big(\dfrac1\rho\dfrac{d\rho}{dz} - \dfrac1\rho\dfrac{d\rho_a}{dz}\Big)",
       "Distribute 1/ρ over both terms; at the rest height the parcel's density equals ρ, so each term is a relative "
       "(logarithmic) rate. Relative rates are what logarithms make easy.",
       "N² compares relative density changes."),
      ("Log-differentiate the environment's p = ρRT", r"\dfrac1\rho\dfrac{d\rho}{dz} = \dfrac1p\dfrac{dp}{dz} - \dfrac1T\dfrac{dT}{dz}",
       "ρ = p/(RT), so ln ρ = ln p − ln R − ln T; differentiate in z with R constant (P52).",
       "The environment's density falls because the pressure falls, and rises where it gets colder."),
      ("Log-differentiate the parcel's isentropic density", r"\dfrac1{\rho_a}\dfrac{d\rho_a}{dz} = \dfrac1\gamma\,\dfrac1p\dfrac{dp}{dz}",
       "Along its isentrope ρ_a ∝ p^(1/γ) (1.26), and the parcel's pressure is the environment's p(z); so ln ρ_a = (1/γ) ln p "
       "+ const.",
       "The parcel's density follows the pressure only, and more weakly."),
      ("Subtract step 3 from step 2",
       r"\dfrac1\rho\dfrac{d\rho}{dz} - \dfrac1\rho\dfrac{d\rho_a}{dz} = \Big(1 - \dfrac1\gamma\Big)\dfrac1p\dfrac{dp}{dz} - \dfrac1T\dfrac{dT}{dz}",
       "Collect the two pressure terms; ρ_a = ρ at the rest height, so the second term on the left is step 3.",
       "What decides stability: a pressure part and a temperature part."),
      ("Log-differentiate θ = T(p_o/p)^((γ−1)/γ)",
       r"\dfrac1\theta\dfrac{d\theta}{dz} = \dfrac1T\dfrac{dT}{dz} - \dfrac{\gamma-1}{\gamma}\,\dfrac1p\dfrac{dp}{dz}",
       "ln θ = ln T + ((γ − 1)/γ)(ln p_o − ln p) (P52); p_o is a constant. This is the log-derivative of (1.31) (N19).",
       "θ rises with T and with falling pressure."),
      ("Recognise the bracket of step 4",
       r"\Big(1 - \dfrac1\gamma\Big)\dfrac1p\dfrac{dp}{dz} - \dfrac1T\dfrac{dT}{dz} = -\,\dfrac1\theta\dfrac{d\theta}{dz}",
       "1 − 1/γ = (γ − 1)/γ, so the right side of step 4 is exactly minus the right side of step 5.",
       "The density race is the θ gradient in disguise."),
      ("Substitute into step 1", r"N^2 = \dfrac{g}{\theta}\,\dfrac{d\theta}{dz}",
       "Steps 4 and 6 replace the bracket by −(1/θ)dθ/dz; −g × (−(1/θ)dθ/dz) — the two minus signs cancel.",
       "N² is g over θ times θ's gradient."),
      ("Use (1.32) for the lapse-rate form", r"N^2 = \dfrac{g}{T}\Big(\dfrac{dT}{dz} + \dfrac{g}{C_p}\Big) = \dfrac{g}{T}\,(\Gamma - \Gamma_a)",
       "(1.32) gives (1/θ)dθ/dz = (1/T)(dT/dz + g/C_p) (C56); with Γ_a = −g/C_p (D19) the bracket is Γ − Γ_a in Kundu's "
       "convention.",
       "N² is proportional to how much more slowly the environment cools than a parcel."),
  ],
  result=(r"N^2 = \dfrac{g}{\theta}\dfrac{d\theta}{dz} = \dfrac{g}{T}\,(\Gamma - \Gamma_a)\ \text{(Kundu)} = \dfrac{g}{T}\,(\Gamma_{a,\rm met} - \Gamma_{\rm met})\ \text{(meteorology)}",
          "a perfect-gas atmosphere is stable exactly where θ increases upward; in lapse rates, where the environment cools "
          "more slowly than 9.76 K/km."),
  interpret="""One profile — θ(z) — tells you everything about dry static stability; this is how Ch. 13 and weather analysis
read soundings, and why isentropic surfaces are the natural coordinates of dry dynamics. It fails for moist air (use the
equivalent potential temperature) and for the ocean, where salinity requires potential density (C61). In the explainer
`parcel_stability` below, the N² computed from θ equals the N² readout computed from the lapse rate.""",
  check="""Units: (m s⁻²/K) × K/m = s⁻² ✓. Neutral: dθ/dz = 0 ⇔ Γ = Γ_a ✓. An isothermal layer (dT/dz = 0) at 250 K:
N² = g²/(C_pT) = 3.83×10⁻⁴ s⁻², period 5.4 min ✓. Standard troposphere at 288.15 K: 1.11×10⁻⁴ s⁻², period 9.9 min ✓ — the
code above computed N² from the θ profile of USSA-1976 with `np.gradient` and from the lapse rate, and they agree to 10⁻³.""",
  traps="using the environment's density gradient for the parcel; forgetting (γ − 1)/γ = R/C_p when comparing with (1.32); "
        "substituting dp/dz = −ρg too early (it is not needed until step 8).")
note("C57", "Stability from θ", r"""
So a dry atmosphere is **stable where θ increases with height, neutral where it stays constant, unstable where it decreases** — the gradient of
θ, not of T, decides. In a laboratory tank the difference hardly matters: over 1 m the adiabatic cooling is only 9.8×10⁻³ K,
so T ≈ θ.
""", equation=r"N^2 = \frac{g}{\theta}\frac{d\theta}{dz}")
note("N18", "A typical lower atmosphere", """
A typical lower atmosphere (our synthetic profile, not digitised from the book's figure): a **mixed layer** near the
ground cooling at almost exactly Γ_a (near-neutral, stirred by turbulence), an **inversion** where T *rises* with height
(very stable: smog and fog stay below it), and a **stable layer** above that cools more slowly than Γ_a. In θ: vertical in
the mixed layer, a sharp increase through the inversion, a steady increase above. The slider below builds it with the mixed
layer set exactly to Γ_a.
""")
nb.plotly("""
z_bl = np.linspace(0, 2000, 201)                                                  # heights up to 2 km, every 10 m [m]
inv = np.linspace(-5e-3, 20e-3, 26)                                               # inversion-layer dT/dz, −5 … +20 K/km [K/m]
def boundary_layer(x_km):                                                         # one slider position (K/km)
    col = ch01.synthetic_boundary_layer_column(z_bl, inversion_dT_dz=x_km * 1e-3,   # T, p, ρ, θ, N² layer by layer …
                                               mixed_dT_dz=ch01.adiabatic_lapse_rate())   # … with the mixed layer exactly neutral
    stab = col["N2"] > 1e-12                                                      # stable where N² > 0
    return {"θ(z)": (col["theta"], z_bl / 1e3),                                   # potential temperature [K] against z [km]
            "θ where stable (N² > 0)": (np.where(stab, col["theta"], np.nan), z_bl / 1e3),   # teal markers
            "T(z)": (col["T"], z_bl / 1e3),                                       # actual temperature [K]
            "θ where neutral or unstable (N² ≤ 0)": (np.where(~stab, col["theta"], np.nan), z_bl / 1e3)}   # rose markers
fig = slider_figure(boundary_layer, "inversion dT/dz", inv * 1e3, unit="K/km", xlabel="temperature T, θ [K]",   # 26 positions
                    ylabel="height z [km]", title="", active=15,                  # start at +10 K/km
                    modes={"θ where stable (N² > 0)": "markers", "θ where neutral or unstable (N² ≤ 0)": "markers"})   # masked θ as dots
palette = {"θ(z)": COLORS["accent"], "T(z)": COLORS["ink"], "θ where stable (N² > 0)": COLORS["teal"],   # colour per trace name
           "θ where neutral or unstable (N² ≤ 0)": COLORS["rose"]}               # the colours used throughout this section
for tr in fig.data:                                                               # recolour every precomputed trace
    tr.line.color = palette[tr.name]; tr.marker.color = palette[tr.name]; tr.marker.size = 4   # line, marker colour and size
for i, x in enumerate(inv):                                                       # per-position title: N² in the inversion
    N2_inv = ch01.brunt_vaisala_sq_from_lapse(288.0, x)                            # at a typical temperature there
    kind, secs = ch01.stability_timescale(N2_inv)                                  # period or e-folding time
    k = ch01.lapse_rate_stability(x, prefix=False); m = ch01.lapse_rate_stability(x, convention="meteorology", prefix=False)  # both forms
    tail = f"period {secs/60:.1f} min" if kind == "period" else (f"e-folds in {secs:.0f} s" if kind == "efold" else "neutral")  # time scale text
    st = fig.layout.sliders[0].steps[i]                                            # the slider step (method "update")
    st.args = [st.args[0], {"title.text": f"inversion layer: {k.text} ⇔ {m.text} ({k.verdict}), N² = {N2_inv:.1e} s⁻², {tail}"}]  # its title
fig.layout.title.text = fig.layout.sliders[0].steps[15].args[1]["title.text"]    # starting title
fig.show()                                                                        # drag the inversion strength
""", explain="""
1. `synthetic_boundary_layer_column` builds T, p, ρ, θ and N² layer by layer (hydrostatics and p = ρRT) for a mixed layer
   at exactly Γ_a, an inversion of adjustable strength from 800 to 1000 m, and a stable layer above.
2. θ is drawn twice more as markers, masked with `np.where` (P46): teal where N² > 0, rose where N² ≤ 0.
3. Each slider position's title gives the inversion layer's criterion in both conventions and its N².
""")
see_read_change(
    'Temperature T (black) and potential temperature θ (purple) against height up to 2 km; θ is overlaid with teal dots where the air is stable and rose dots where it is neutral or unstable.',
    'In the mixed layer θ is vertical (neutral, rose); through the inversion θ jumps to the right (very stable, teal); above it θ keeps increasing slowly (stable). Stability is read from the tilt of θ, not of T: the upper layer is stable although T falls with height.',
    "Drag the inversion strength below −9.8 K/km: that layer's θ starts to decrease with height, its dots turn rose, and the title reports an e-folding time instead of a period.",
)
note("C58", "Potential density", r"""
The same trick for density: the **potential density** ρ_θ is the density a parcel would have after an isentropic trip to
p_o (stated). Air with ρ = 0.70 kg/m³ at 500 hPa has ρ_θ = 1.148 kg/m³ (code below).
""", equation=r"\rho(z) = \rho_\theta(z)\,\big(p(z)/p_o\big)^{1/\gamma}", ref="1.33")
note("N20", "A constant product", """
Multiplying (1.31) and (1.33) and using p = ρRT, the exponents add up to 1 and θρ_θ = p_o/R is the same for every parcel — one
line of algebra, not used again (checked in the code below).
""")
note("C59", "Potential density must decrease upward", """
Taking logarithms of θρ_θ = const gives (1.34), −(1/ρ_θ)dρ_θ/dz = (1/θ)dθ/dz: for stability the potential density must
**decrease** with height — the form the ocean uses (C61), and Ch. 13 §13.2 uses for vertical density variation.
""")
nb.code("""
print(f"C58 ρθ for ρ = 0.70 kg/m³ at 500 hPa: {ch01.potential_density(0.7, 5.0e4):.4f} kg/m³")   # Eq. (1.33)
rho_th = ch01.potential_density(rho_tr, p_tr)                                                     # ρθ along the troposphere
assert np.allclose(theta_tr * rho_th, ch01.P_REF / ch01.R_AIR)                                    # N20: θ ρθ = p_o / R everywhere
print(f"N20: θ ρθ = {theta_tr[0]*rho_th[0]:.3f} = p_o/R = {ch01.P_REF/ch01.R_AIR:.3f} kg K/m³")   # the constant product
print("C59: ρθ decreases upward:", bool(np.all(np.diff(rho_th) < 0)))                             # stable troposphere
""", explain="""
1. `potential_density` is (1.33) solved for ρ_θ.
2. The product θρ_θ is the same constant p_o/R at every height of the standard troposphere (N20).
3. ρ_θ decreases monotonically upward, the density form of the stability test (C59).
""")
nb.explainer("parcel_stability", heading="Push a parcel up: does it come back?", why="""
Stability links three things at once — the environment's profile, the parcel's own adiabat and ζ(t) — and the sign
convention adds a fourth. Reshaping the profile while the parcel moves, and flipping Kundu ↔ meteorology while nothing
physical changes, makes the criterion stick. Its Derivation tab walks through D18, D19, D20 and D36.
""", tries=[
    "Start from 'standard atmosphere', press play, then switch the convention chip and read the status line in both forms.",
    "Drag dT/dz below −9.8 K/km and watch the cosine turn into a cosh.",
    "Turn on the nocturnal inversion and release the parcel inside it.",
    "Ocean mode: toggle compressibility and compare the periods (Check question 4).",
])
nb.md("""
**What would change if…** the air were moist and clouds formed? Condensation releases heat, so a saturated parcel cools more
slowly than Γ_a (about 4–7 K/km instead of 9.8, depending on temperature), and a column that is stable for dry air can be unstable for cloudy air —
beyond this chapter. **Next:** §1.11 asks what dimensional reasoning alone can tell us.
""")

# =====================================================================================================================
# A.11 §1.11 Dimensional Analysis — C64, C67, C69
# =====================================================================================================================
nb.section("1.11", "Dimensional Analysis", intro="""
**What is this section about?** Because nature does not care about our units, every correct equation has matching
dimensions, and every law can be written with dimensionless groups — exactly n − r of them. That turns a seven-variable
experiment into a four-number one.
""")
# ---------------------------------------------------------------------------------------------------- C64
core("C64", "Dimensional homogeneity: a correct law cannot care about units",
     "Could a true law of nature give a different answer in feet than in metres?")
nb.md("""
#### The problem in plain words
A spacecraft was lost because one team worked in pound-force seconds and another in newton seconds. Checking that every term
of an equation has the same dimensions catches such mistakes — and every test in this project does it with `pint`. The same
principle lets us predict the form of laws before solving anything.
""")
nb.md(r"""
#### The idea
**Units** are our choice (m, ft); **dimensions** are what is measured (a length). Write a dimension as powers of the base
dimensions, $[q] = \mathrm{M}^a\,\mathrm{L}^b\,\mathrm{T}^c\,\Theta^d$ — an exponent vector (a, b, c, d). Adding two terms is
only meaningful when their vectors are equal.

| quantity | M | L | T | Θ |
|---|---|---|---|---|
| p | 1 | −1 | −2 | 0 |
| ρ | 1 | −3 | 0 | 0 |
| g | 0 | 1 | −2 | 0 |
| z | 0 | 1 | 0 | 0 |

Reminders: R01 (§1.2) — a unit is a product of powers of base units; R02 — prefixes only rescale numbers, they never change
a dimension.
""")
nb.md(r"""
#### The maths, step by step
1. $[\rho g z]$ = (1, −3, 0, 0) + (0, 1, −2, 0) + (0, 1, 0, 0) = (1, −1, −2, 0) — exponents **add** when quantities multiply
   (P43).
2. $[p_0]$ = (1, −1, −2, 0): the same, so $p = p_0 - \rho g z$ can be right.
3. A wrong formula $p = p_0 - \rho g z^2$ has (1, 0, −2, 0) ≠ [p]: rejected without any experiment.
4. Divide the correct law by $p_0$: $p/p_0 = 1 - \rho g z/p_0$ — every term is now a pure number, the same in any system of
   units. Two conclusions: the terms of a law must match (**dimensional homogeneity**), and every law can be written in
   **dimensionless form**.
""")
nb.worked_example("10 m of water in SI and in cgs", """
- SI: ρgz/p₀ = 1000 × 9.81 × 10 / 101 325 = 0.968.
- cgs: ρ = 1 g/cm³, g = 981 cm/s², z = 1000 cm, p₀ = 1 013 250 dyn/cm² → 1 × 981 × 1000 / 1 013 250 = 0.968.

The group is identical, although the individual numbers changed by factors of 1000, 100 and 10.
""")
nb.code("""
import pint                                                                         # the units library behind Q_ (primer P02)
print(ch01.dimension_vector("kg/m**3") + ch01.dimension_vector("m/s**2") + ch01.dimension_vector("m"))  # [ρgz] → [ 1 -1 -2  0]
print(ch01.dimension_vector("Pa"))                                                  # [p] → [ 1 -1 -2  0]: the same
print((Q_(1000, "kg/m**3") * Q_(9.81, "m/s**2") * Q_(10, "m")).to("Pa"))           # ρgz converts to pascals: 98100 Pa
try:                                                                                # the wrong formula ρ g z² …
    (Q_(1000, "kg/m**3") * Q_(9.81, "m/s**2") * Q_(10, "m")**2).to("Pa")   # kg/s² cannot become Pa
except pint.DimensionalityError as err:                                             # … is refused by pint
    print("rejected:", str(err)[:70])   # pint's message
V = {"rho": "kg/m**3", "g": "m/s**2", "z": "m", "p0": "Pa"}                         # the variables and their SI units
vals = {"rho": 1000.0, "g": 9.81, "z": 10.0, "p0": 101325.0}                        # their values in SI
group = {"rho": 1, "g": 1, "z": 1, "p0": -1}                                        # the group ρ g z / p0 as exponents
for system in ("SI", "cgs", "imperial"):                                            # the same physical values in three unit systems
    print(system, round(ch01.group_value(group, ch01.rescale_units(vals, V, system)), 4))   # 0.9682 every time
""", explain="""
1. `dimension_vector` returns the (M, L, T, Θ) exponents of a unit; adding the vectors multiplies the quantities.
2. `pint` converts ρgz to pascals but raises `DimensionalityError` for ρgz², whose dimensions are not a pressure.
3. `rescale_units` expresses the same physical values in another unit system, and `group_value` multiplies the powers of a
   group: the group ρgz/p₀ has the same value in SI, cgs and imperial units.
""")
nb.figure("""
basis = ["M", "L", "T", "Θ"]                                                       # base dimensions
vec = lambda u: ch01.dimension_vector(u)                                           # exponent vector of a unit
terms_ok = {"p": vec("Pa"), "p₀": vec("Pa"), "ρgz": vec("kg/m**3") + vec("m/s**2") + vec("m")}          # correct law
terms_bad = {"p": vec("Pa"), "p₀": vec("Pa"), "ρgz²": vec("kg/m**3") + vec("m/s**2") + 2 * vec("m")}    # wrong law: z²
cols = [COLORS["orange"], COLORS["blue"], COLORS["teal"], COLORS["amber"]]         # M orange, L blue, T teal, Θ amber
fig, axes = plt.subplots(1, 2, figsize=(9, 3.4), sharey=True)                      # left correct, right wrong
for ax, terms, title in ((axes[0], terms_ok, "p = p₀ − ρgz  ✓"), (axes[1], terms_bad, "p = p₀ − ρgz²  ✗")):   # the two panels
    ref = terms["p"]                                                                # the dimension every term must have
    for j, (name, e) in enumerate(terms.items()):                                   # one group of four bars per term
        for b in range(4):                                                          # one bar per base dimension
            bad = e[b] != ref[b]                                                    # this exponent does not match
            x_b = j + (b - 1.5) * 0.2                                               # bar position
            ax.bar(x_b, e[b], width=0.18, color=cols[b], label=basis[b] if j == 0 else None)   # the bar (height = exponent)
            if e[b] == 0:                                                           # a zero exponent has no height …
                ax.plot([x_b - 0.08, x_b + 0.08], [0, 0], color=cols[b], lw=4)      # … so draw a short thick stub instead
            if bad:                                                                 # this exponent differs from p's
                ax.bar(x_b, ref[b], width=0.18, fill=False, edgecolor=COLORS["rose"], lw=2, ls="--")   # where it should be
                ax.annotate(f"{basis[b]}: {e[b]} ≠ {ref[b]}".replace("-", "−"), (x_b, ref[b]), (x_b - 0.9, -2.6),   # label it
                            color=COLORS["rose"], fontsize=10, arrowprops=dict(arrowstyle="->", color=COLORS["rose"]))
    ax.set_xticks(range(3)); ax.set_xticklabels(list(terms)); ax.axhline(0, color=COLORS["muted"], lw=0.8)   # term names, zero line
    ax.set_title(title)                                                             # which formula
axes[0].set_ylim(-3.0, 1.5)                                                         # room for the annotation
axes[0].set_ylabel("exponent of the base dimension [–]")                            # shared y label
axes[0].legend(ncol=4, fontsize=8, loc="lower left")                                # which colour is which dimension
plt.show()                                                                          # display
""", see="""
Two panels of grouped bars — the exponents of M (orange), L (blue), T (teal) and Θ (amber) for each term of a formula; an
exponent of 0 (Θ everywhere, since no temperature appears) is drawn as a short flat stub. On the left every term has the same
pattern (1, −1, −2, 0). On the right the ρgz² term has L exponent 0 where p has −1: its L bar is only a stub, the dashed
rose outline shows where it should reach, and the label reads "L: 0 ≠ −1".
""", read="""
Equal bars for every term mean the formula is dimensionally homogeneous and may be right; one mismatched bar is enough to
reject a formula without any experiment.
""", change="""
Add a term ½ρu² (the dynamic pressure) to the left panel: its bars would match p's too — Bernoulli's equation (Ch. 4) passes
this test.
""")
nb.md("""
**What would change if…** a problem had seven variables? Checking one equation is easy; finding *all* the unit-free
combinations needs bookkeeping. **Next:** put the exponent vectors side by side as the columns of a matrix (C67).
""")
# ---------------------------------------------------------------------------------------------------- C67
core("C67", "The dimensional matrix: dimensions as columns of numbers",
     "How do we write down what dimensions each variable has, so that linear algebra can reason about them?")
nb.md("""
#### The problem in plain words
An engineer wants the pressure drop along a pipe for any fluid, diameter, roughness and speed. Seven variables at ten values
each would be 10⁷ experiments. Dimensional analysis promises far fewer — and it starts by tabulating dimensions.
""")
nb.md("""
#### The idea
Stack the exponent vectors of C64 as columns: the rows are the base dimensions, the columns the variables. Everything
dimensional analysis needs is now in one small integer matrix.
""")
note("C65", "Step 1 — choose the variables", r"""
The most important step: one solution variable (Δp) plus everything it can depend on — geometry (Δx, d, ε), flow (U) and
material (ρ, μ). Leave one out and no amount of algebra recovers it (dropping μ is tested in the C69 collapse figure); add
unnecessary ones and the result gets weaker.
""", equation=r"f(\Delta p, \Delta x, d, \varepsilon, U, \rho, \mu) = 0", ref="1.38")
note("C66", "Step 2 — dimensions", """
[q] means "the dimension of q", written with M, L, T and Θ: [U] = L/T, [p] = M/(LT²), [C_p] = L²/(ΘT²). Temperature often
appears only as k_BΘ, RΘ or C_pΘ, whose dimension is L²/T², so three base dimensions suffice for most flows. (The kilomole in
R_u is treated as a pure number, as the book's Example 1.2 does; `dimension_vector` says so with a warning.)
""")
P("P53", "matrices, determinants and minors", r"""
A matrix is a table of numbers; the **determinant** of a square one is a single number that is zero exactly when its columns
are dependent. A **minor** is the determinant of a square piece cut out of a bigger matrix. 3×3 determinants by **cofactor
expansion** along the first row: $a_{11}(a_{22}a_{33} - a_{23}a_{32}) - a_{12}(a_{21}a_{33} - a_{23}a_{31}) +
a_{13}(a_{21}a_{32} - a_{22}a_{31})$.
""", code="""
print(np.linalg.det(np.array([[2, 0, 0], [0, 3, 0], [0, 0, 4]])))   # a diagonal matrix: det = 2 × 3 × 4 = 24.0
""")
P("P54", "linear independence and rank", """
Vectors are independent if none is a combination of the others. The **rank** of a matrix is the number of independent
columns (which equals the number of independent rows) — the size of its largest nonzero minor.
""", code="""
print(np.linalg.matrix_rank(np.array([[1, 2], [2, 4]])))   # the second column is 2 × the first → rank 1
""")
nb.md(r"""
#### The maths, step by step
1. Columns in the order of (1.38): Δp, Δx, d, ε, U, ρ, μ.
2. [Δp] = M L⁻¹ T⁻² → column (1, −1, −2); [Δx] = [d] = [ε] = L → (0, 1, 0); [U] = L T⁻¹ → (0, 1, −1); [ρ] = M L⁻³ →
   (1, −3, 0); [μ] = M L⁻¹ T⁻¹ → (1, −1, −1).
3. The matrix (no Θ row is needed):
$$\begin{array}{c|ccccccc} & \Delta p & \Delta x & d & \varepsilon & U & \rho & \mu\\\hline \mathrm M & 1&0&0&0&0&1&1\\ \mathrm L & -1&1&1&1&1&-3&-1\\ \mathrm T & -2&0&0&0&-1&0&-1\end{array} \qquad (1.39)$$
4. Its rank r is at most 3, because it has three rows.
""")
note("C68", "Step 3 — rank by minors", """
r is the size of the largest square sub-matrix with a nonzero determinant. The minor of the first three columns (Δp, Δx, d)
is 0 — Δx and d are both pure lengths, so those columns are dependent — but the minor of the last three (U, ρ, μ) is −1 ≠ 0,
so r = 3 (worked by hand below). A rank below the number of rows happens when one row is a combination of the others, for
example in statics problems where mass enters only through force.
""")
nb.worked_example("the (U, ρ, μ) minor by hand", """
Rows M, L, T of the columns U, ρ, μ: [[0, 1, 1], [1, −3, −1], [−1, 0, −1]].

Cofactor expansion along row 1: 0·((−3)(−1) − (−1)(0)) − 1·((1)(−1) − (−1)(−1)) + 1·((1)(0) − (−3)(−1))
= 0 − 1·(−1 − 1) + (0 − 3) = 2 − 3 = **−1** ≠ 0 → r = 3.

The first-three-columns minor [[1, 0, 0], [−1, 1, 1], [−2, 0, 0]]: expansion gives 1·(1·0 − 1·0) − 0 + 0 = 0.
""")
P("P55", "itertools.combinations", """
`combinations(range(7), 3)` lists every way to choose 3 of 7 column indices, without repeats and ignoring order — 35 of them.
""", code="""
import itertools                                          # standard library: combinatorics
print(len(list(itertools.combinations(range(7), 3))))     # 7·6·5/(3·2·1) = 35 ways to pick 3 columns
""")
nb.code("""
A, names, rows = ch01.dimensional_matrix(ch01.PIPE)                  # Eq. (1.39): integer matrix, column names, row names
print(rows, names)                                                   # ['M', 'L', 'T'] ['dp', 'dx', 'd', 'eps', 'U', 'rho', 'mu']
print(A)                                                             # the 3 × 7 matrix
r, ri, ci = ch01.rank_by_minors(A)                                   # rank and the first nonzero r × r minor found
print(f"rank {r}; first nonzero minor: rows {ri}, columns {[names[c] for c in ci]}")   # r = 3 and its witness
print("det of (U, ρ, μ):", ch01.minor_determinant(A, (0, 1, 2), (4, 5, 6)))    # −1, exact integer arithmetic
print("det of (Δp, Δx, d):", ch01.minor_determinant(A, (0, 1, 2), (0, 1, 2)))  # 0
""", explain="""
1. `dimensional_matrix` builds (1.39) from the dictionary of variables and units (`ch01.PIPE`), dropping the empty Θ row.
2. `rank_by_minors` searches the minors in order and returns the rank with the first nonzero witness — here the columns
   (Δp, Δx, U); any nonzero 3×3 minor proves r = 3.
3. `minor_determinant` evaluates a chosen minor exactly by cofactor expansion: −1 for (U, ρ, μ), 0 for (Δp, Δx, d).
""")
P("P56", "np.linalg.det and np.linalg.matrix_rank", """
numpy's floating-point versions: `det` returns, for example, −0.9999999999999998 for an integer −1 (round-off), so round it
before comparing; `matrix_rank` counts the independent columns but, working in floating point, ignores directions smaller than a tiny tolerance so that round-off is not mistaken for independence.
""", code="""
print(np.linalg.det(A[:, 4:7]), round(np.linalg.det(A[:, 4:7])))   # float round-off, then the exact integer −1
""")
nb.md("#### From scratch: every 3×3 minor by cofactor expansion")
nb.check_agree("""
def det3(M):                                                        # determinant of a 3×3 matrix by cofactor expansion (P53)
    return (M[0][0] * (M[1][1] * M[2][2] - M[1][2] * M[2][1])   # a11 (a22 a33 − a23 a32)
            - M[0][1] * (M[1][0] * M[2][2] - M[1][2] * M[2][0])   # − a12 (a21 a33 − a23 a31)
            + M[0][2] * (M[1][0] * M[2][1] - M[1][1] * M[2][0]))   # + a13 (a21 a32 − a22 a31)
all_c = list(itertools.combinations(range(7), 3))                   # the 35 choices of three columns
nonzero = [c for c in all_c if det3(A[:, c]) != 0]                  # the minors that are not zero
assert len(nonzero) > 0                                             # at least one → rank 3
assert 3 == ch01.rank_by_minors(A)[0] == np.linalg.matrix_rank(A)   # hand, library and numpy agree on r = 3
assert all(det3(A[:, c]) == round(np.linalg.det(A[:, c])) for c in all_c)   # every minor matches numpy
print(f"{len(nonzero)} of {len(all_c)} minors are nonzero; rank = 3")   # how many minors survive
""")
nb.figure("""
symbols = ["Δp", "Δx", "d", "ε", "U", "ρ", "μ"]                                      # column labels
fig, ax = plt.subplots(figsize=(7.4, 3.2))                                          # one heatmap
ax.imshow(A, cmap="RdBu", vmin=-3, vmax=3)                                          # imshow colours each matrix cell
for i in range(A.shape[0]):                                                         # print the integer in every cell
    for j in range(A.shape[1]):   # every column
        ax.text(j, i, f"{A[i, j]:d}".replace("-", "−"), ha="center", va="center", fontsize=12,   # the integer, with a true minus sign
                color="white" if abs(A[i, j]) >= 2 else COLORS["ink"])             # white text on the dark cells
ax.set_xticks(range(7)); ax.set_xticklabels(symbols); ax.set_yticks(range(3)); ax.set_yticklabels(rows)   # label columns and rows
ax.grid(False)                                                                      # no grid lines over the cells
ax.add_patch(plt.Rectangle((3.5, -0.5), 3, 3, fill=False, ec=COLORS["accent"], lw=3))   # the (U, ρ, μ) minor
ax.text(5, 2.75, "det = −1", color=COLORS["accent"], ha="center", va="top", fontsize=10)   # its determinant
ax.add_patch(plt.Rectangle((-0.5, -0.5), 3, 3, fill=False, ec=COLORS["muted"], lw=3, ls="--"))   # the (Δp, Δx, d) minor
ax.text(1, 2.75, "det = 0", color=COLORS["muted"], ha="center", va="top", fontsize=10)   # its determinant
ax.set_ylim(3.0, -0.6)                                                              # leave room for the labels below
ax.set_title(f"Eq. (1.39): one nonzero 3×3 minor proves r = 3 ({len(nonzero)} of 35 are nonzero)")   # the message
savefig(fig, "ch01", "c67_dimensional_matrix")                                      # keep a copy for review
plt.show()                                                                          # display
""", see="""
The 3 × 7 integer matrix (1.39) as coloured cells, with the (U, ρ, μ) block outlined in purple (determinant −1) and the
(Δp, Δx, d) block outlined in grey (determinant 0).
""", read="""
One nonzero 3×3 block is enough to prove r = 3. The three identical pure-length columns (Δx, d, ε) are what make many minors
vanish: any block containing two of them has two equal columns.
""", change="""
Add the temperature T as an eighth variable and a Θ row appears, so the rank could reach 4 — Example 1.2 (C73) has r = 4.
""")
nb.md("""
**What would change if…** we asked which combinations of exponents make a product of the variables dimensionless? That is
solving A·k = 0 — the null space of this matrix, and the Π theorem (C69).
""")
# ---------------------------------------------------------------------------------------------------- C69
core("C69", "Buckingham's Π theorem: a law depends on n − r numbers",
     "Why can seven pipe-flow variables be squeezed into four dimensionless numbers — and why exactly four?")
nb.md("""
#### The problem in plain words
Wind-tunnel models of aircraft, scale models of harbours, rotating-tank experiments for the atmosphere: a small experiment
stands in for a big system when the right dimensionless numbers match. The Reynolds, Froude and Rossby numbers (Ch. 4 §4.11,
Ch. 13) all come from this theorem.
""")
nb.md(r"""
#### The idea
A product of powers $q_1^{k_1}\cdots q_n^{k_n}$ is dimensionless exactly when its exponent vector **k** solves A·**k** = 0.
The solutions form a subspace — the **null space** — whose dimension is n − r. Each basis vector is one Π group.

```
7 variables  →  matrix of rank 3  →  null space of dimension 7 − 3 = 4  →  Π₁, Π₂, Π₃, Π₄
```
""")
note("N22", "The starting relation", """
The starting point is a relation among the n variables — written in code as the input dictionary of `pi_groups`; Ch. 4 §4.11
starts from the same form.
""", equation=r"f(q_1, q_2, \dots, q_n) = 0", ref="1.36")
note("N23", "Step 4 — count the groups", "The number of independent groups is n − r. Pipe: 7 − 3 = 4.")
P("P58", "null space and rank–nullity theorem", """
The null space of A is the set of all vectors k with A·k = 0 (sums and multiples of solutions are again solutions).
**Rank–nullity:** for a matrix with n columns, rank + (dimension of the null space) = n.
""", code="""
print(sp.Matrix([[1, 1, 0], [0, 0, 1]]).nullspace())   # rank 2, 3 columns → one null-space vector: (−1, 1, 0)
""")
P("P59", "scaling the base units", """
Change the unit of length by a factor λ_L (m → cm: λ_L = 100), of mass by λ_M, of time by λ_T. A quantity with dimension
M^a L^b T^c gets its **number** multiplied by λ_M^a λ_L^b λ_T^c; the physical quantity itself is unchanged.
""", code="""
print(1000 * 1000**1 * 100**-3)    # 1000 kg/m³ in g/cm³: λ_M = 1000 (kg → g), λ_L = 100 (m → cm), dimension M L⁻³ → 1.0
""")
P("P61", "sympy Matrix and nullspace", """
`sp.Matrix(...)` keeps integers exact; `.rank()` and `.nullspace()` return exact results — fractions, never round-off.
""", code="""
print(sp.Matrix([[1, 0, 0, 0, 0, 1, 1], [-1, 1, 1, 1, 1, -3, -1], [-2, 0, 0, 0, -1, 0, -1]]).rank())   # (1.39) → 3
""")
D("D28", "Buckingham's Π theorem", ref="1.37",
  goal="""Show that any correct relation among n dimensional variables can be rewritten as a relation among exactly n − r
independent dimensionless groups, with r the rank of the matrix A of dimensions — and see why the groups are not unique.""",
  assumptions="""The variable list is complete (step 13: nothing else can enter f) · f is dimensionally homogeneous, i.e. true
in every unit system (steps 12–13) · all dimensions are products of powers of the base dimensions (step 2).""",
  start=(r"f(q_1, q_2, \dots, q_n) = 0 \ \ (1.36), \text{ with dimensional matrix } A \text{ of rank } r",
         "a complete, dimensionally homogeneous relation among n variables whose dimensions we have tabulated (rows = base "
         "dimensions, columns = variables)."),
  plan=["A product of powers is dimensionless exactly when its exponent vector solves A·k = 0.",
        "Rank–nullity counts the independent solutions: n − r.",
        "Choosing the units cleverly shows that the relation can depend only on those groups."],
  uses=["dimensional homogeneity (C64)", "the dimensional matrix and its rank (C67, C68)", "exponent rules (P43)",
        "linear independence and rank (P54)", "null space and rank–nullity (P58)", "scaling the base units (P59)",
        "matrices and minors (P53)"],
  steps=[
      ("Form a general product of powers", r"\Pi = q_1^{k_1}\,q_2^{k_2}\cdots q_n^{k_n}",
       "We *look for* unit-free combinations of this single-product (monomial) form, because a product of powers of the "
       "variables has a dimension that is again a product of powers of M, L, T, Θ, which we can set to 1. Steps 9–13 show "
       "that nothing more general is ever needed. The exponents k = (k₁, …, k_n) are what we must find.",
       "A group is some variables multiplied together, each raised to a power."),
      ("Write its dimension with the matrix entries", r"[\Pi] = \prod_d d^{\,\sum_j A_{dj} k_j}",
       "[q_j] = M^{A_Mj} L^{A_Lj} T^{A_Tj} Θ^{A_Θj}; raising to k_j multiplies the exponents and multiplying quantities adds "
       "them (P43).",
       "The exponent of each base dimension is a weighted sum of the k's."),
      ("Require every exponent to vanish", r"A\,\mathbf k = \mathbf 0",
       "Dimensionless means M⁰L⁰T⁰Θ⁰; each row d of A gives one equation Σ_j A_dj k_j = 0.",
       "A group is dimensionless exactly when its exponent vector solves this linear system."),
      ("Name the solution set", r"\{\mathbf k : A\mathbf k = \mathbf 0\} = \operatorname{null}(A)",
       "Definition of the null space (P58); sums and multiples of solutions are solutions, so it is a vector space.",
       "All dimensionless groups live in the null space of A."),
      ("Count its dimension", r"\dim\operatorname{null}(A) = n - r",
       "The rank–nullity theorem for a matrix with n columns: rank + nullity = n (P58).",
       "There are exactly n − r independent exponent vectors."),
      ("Pick a basis and call the groups Π₁ … Π_{n−r}", r"\mathbf k^{(1)},\dots,\mathbf k^{(n-r)} \;\to\; \Pi_1,\dots,\Pi_{n-r}",
       "Any basis of the null space will do; its vectors are linearly independent (P54), so no group is a product of powers "
       "of the others (see C71, the note just after this derivation). Different bases give different, equivalent sets — the groups are not unique.",
       "n − r independent dimensionless numbers."),
      ("Show every group is built from these", r"\mathbf k = \sum_i c_i\mathbf k^{(i)} \;\Rightarrow\; \Pi = \prod_i \Pi_i^{\,c_i}",
       "A basis spans the space; adding exponent vectors multiplies the corresponding powers (P43).",
       "Any other dimensionless combination is a product of powers of Π₁ … Π_{n−r}."),
      ("Choose r variables with independent columns", r"\det A_{\rm rep} \neq 0",
       "Rank r means some r×r minor is nonzero (C68); its r variables are the repeating set. If A has more rows than r, drop "
       "the dependent rows first.",
       "Pick r variables that between them carry all the independent dimensions."),
      ("Build one group for each other variable",
       r"A_{\rm rep}\,\mathbf a_j = -A_{\cdot j}\ \Rightarrow\ \Pi_j = q_j \prod_{i\in\text{rep}} q_i^{\,a_{ji}}\quad (j \notin \text{rep})",
       "The exponent vector of Π_j has 1 in place j, the unknowns a_j in the repeating places and 0 elsewhere, so A·k = 0 "
       "(step 3) becomes the square system A_rep a_j = −A_{·j} (the exponent algebra worked with numbers in C70 below). It has exactly one solution because "
       "det A_rep ≠ 0 (step 8). This gives n − r groups; each contains its own q_j and no other non-repeating variable, so "
       "their exponent vectors are independent (P54), and n − r independent vectors in a space of dimension n − r (step 5) "
       "form a basis of null(A).",
       "Every non-repeating variable, made dimensionless with the repeating ones, is a group — and together they are a basis."),
      ("Rescale the base units", r"q_j \to q_j\prod_d \lambda_d^{A_{dj}}",
       "Changing the unit of dimension d by a factor λ_d multiplies the *number* of every quantity by λ_d to its exponent "
       "(P59); the physics is unchanged.",
       "New units, new numbers, same world."),
      ("Choose the scales that make the repeating variables equal to 1", r"\sum_d A_{dj}\ln\lambda_d = -\ln q_j \quad (j \in \text{rep})",
       "Taking logs turns step 10 into r linear equations for the ln λ_d, whose matrix is the (transposed) repeating block — "
       "nonsingular by step 8 — so a solution exists.",
       "We can always pick units in which the r repeating variables are exactly one."),
      ("Read the other variables in those units", r"q_j' = \Pi_j \quad (j \notin \text{rep})",
       "Π_j = q_j × (powers of the repeating variables) from step 9 is dimensionless, so its number is the same in every "
       "unit system (C64); in the new units the repeating factors are 1, so q_j′ equals Π_j.",
       "In these units each remaining variable *is* its group."),
      ("Use homogeneity to drop the ones", r"f(1,\dots,1,\Pi_1,\dots,\Pi_{n-r}) = 0 \;\Rightarrow\; \phi(\Pi_1,\dots,\Pi_{n-r}) = 0",
       "A correct law holds in every unit system (C64), in particular in these; the r constant ones are absorbed into a new "
       "function φ. Completeness of the variable list guarantees that nothing else hides in f.",
       "The law depends only on the n − r groups."),
  ],
  result=(r"\phi(\Pi_1, \Pi_2, \dots, \Pi_{n-r}) = 0 \quad\text{or}\quad \Pi_1 = \varphi(\Pi_2, \dots, \Pi_{n-r})",
          "n variables with r independent dimensions can always be traded for n − r independent dimensionless groups."),
  interpret="""Dimensional analysis is linear algebra: groups are null-space vectors, their number is rank–nullity, and a
different repeating set is a change of basis. It is why experiments are plotted against Re, Fr and Ro, and why a small model
can stand for a big system (Ch. 4 §4.11). It fails — silently — if a relevant variable is missing (drop μ from the pipe: 3
groups, no Reynolds number, and laminar data will not collapse) or if the "law" mixes quantities that are only numerically
related in one unit system. In the explainer `buckingham_pi_machine` below, switch the units to cgs and every Π value stays
put.""",
  check="""Pipe: n = 7, r = 3 → 4 groups (1.40) ✓; Π₁ = Δp/ρU² = 10 in SI and in cgs ✓ (code below); `groups_independent`
rejects a fifth group ✓; Example 1.4 (n = 4, r = 3) → one group, hence a constant ✓. The sympy cell checks the null space
exactly.""",
  traps="believing the groups are unique; forgetting that the relation must be complete and homogeneous; choosing a "
        "singular repeating set (two pure lengths such as d and ε cannot cancel time).",
  check_src="""
import sympy as sp                                            # exact linear algebra (primer P61)
A_s = sp.Matrix([[1, 0, 0, 0, 0, 1, 1],                       # M row of (1.39): columns dp dx d eps U rho mu
                 [-1, 1, 1, 1, 1, -3, -1],                    # L row
                 [-2, 0, 0, 0, -1, 0, -1]])                   # T row
ns = A_s.nullspace()                                          # a basis of all dimensionless exponent vectors (step 4)
print(A_s.rank(), len(ns), A_s.shape[1] - A_s.rank())         # 3 4 4 → rank–nullity (step 5)
print(all((A_s * k).is_zero_matrix for k in ns))              # True → every basis vector is dimensionless (step 3)
k1 = sp.Matrix([1, 0, 0, 0, -2, -1, 0])                       # exponents of Π₁ = Δp U⁻² ρ⁻¹
print((A_s * k1).is_zero_matrix)                              # True → Π₁ is dimensionless
print(sp.Matrix.hstack(*ns, k1).rank() == len(ns))            # True → Π₁ is a combination of the basis (step 7)
rep = [4, 2, 5]                                               # step 8: repeating columns U, d, rho (their minor is nonzero)
A_rep = A_s[:, rep]                                           # the 3×3 repeating block
print(A_rep.det() != 0)                                       # True → step 9's systems have unique solutions
groups_k = []                                                 # exponent vectors built as in step 9
for j in (0, 1, 3, 6):                                        # the non-repeating columns dp, dx, eps, mu
    a_j = A_rep.LUsolve(-A_s[:, j])                           # solve A_rep a_j = −A_{·j} exactly
    k = sp.zeros(7, 1); k[j] = 1                              # 1 in the variable's own place …
    for i, col in enumerate(rep):                             # … and the solved exponents in the repeating places
        k[col] = a_j[i]                                        # exponent of repeating variable col
    groups_k.append(k)                                        # keep this group's exponent vector
print(all((A_s * k).is_zero_matrix for k in groups_k), sp.Matrix.hstack(*groups_k).rank())   # True 4 → a basis of null(A)
""")
note("C70", "Step 5 — build the groups (exponent algebra)", r"""
Choose r = 3 repeating variables whose minor is nonzero (U, d, ρ). Each group is one other variable times powers of them.
For Δp: Π₁ = Δp U^a d^b ρ^c; matching M: c + 1 = 0; L: a + b − 3c − 1 = 0; T: −a − 2 = 0 gives a = −2, b = 0, c = −1
(stated; solved in code below). The same with Δx, ε and μ gives:
""", equation=r"\Pi_1 = \frac{\Delta p}{\rho U^2},\quad \Pi_2 = \frac{\Delta x}{d},\quad \Pi_3 = \frac{\varepsilon}{d},\quad \Pi_4 = \frac{\mu}{\rho U d}")
note("N24", "Groups by inspection", """
The same groups can be found **by inspection**, cancelling M, L and T one ratio at a time — the way Navier–Stokes is made
dimensionless in Ch. 4 §4.11.
""")
note("C71", "Combining groups", """
Products and powers of groups are groups too: Δp d²ρ/μ² = Π₁/Π₄² and ε/Δx = Π₃/Π₂. Any such set may replace the original,
but **only n − r of them are independent** — a fifth is always a combination of four (`groups_independent` → False, below).
""")
note("C72", "Step 6 — the dimensionless law", "Π₄ is the inverse of the Reynolds number, 1/Re.",
     equation=r"\frac{\Delta p}{\rho U^2} = \varphi\!\left(\frac{\Delta x}{d}, \frac{\varepsilon}{d}, \frac{\mu}{\rho U d}\right)",
     ref="1.40")
note("N25", "Step 7 — add physics", """
Far from the inlet the pressure drop is proportional to the length, so Δp/ρU² = (Δx/d)·φ₂(ε/d, Re): three numbers become two.
For slow laminar flow in a smooth pipe the full solution (Ch. 8) gives φ₂ = 32/Re — the line the collapse figure below
follows.
""")
nb.worked_example("Π₁ in two unit systems", """
- SI: Δp = 100 Pa, ρ = 1000 kg/m³, U = 0.1 m/s → Π₁ = 100/(1000 × 0.01) = 10.
- cgs: Δp = 1000 dyn/cm², ρ = 1 g/cm³, U = 10 cm/s → Π₁ = 1000/(1 × 100) = 10.
- With d = 1 cm and μ = 10⁻³ Pa s: Π₄ = 10⁻³/(1000 × 0.1 × 0.01) = 10⁻³, i.e. Re = 1000.
""")
P("P60", "fractions.Fraction", """
`Fraction(-1, 2)` is an exact rational number; fluidpy returns exponents this way so that ½ never becomes 0.49999.
""", code="""
from fractions import Fraction              # exact rational numbers (standard library)
print(Fraction(1, 3) + Fraction(1, 6))      # → 1/2 exactly
""")
nb.code("""
groups = ch01.pi_groups(ch01.PIPE, solution="dp", repeating=("U", "d", "rho"))   # (1.40): n − r = 4 groups, Δp in the first
print([ch01.group_latex(g) for g in groups])                                     # their LaTeX forms
print(ch01.solve_exponents("dp", ("U", "d", "rho"), ch01.PIPE))                  # C70: exponents of Π₁ as exact fractions
pipe_vals = {"dp": 100.0, "dx": 1.0, "d": 0.01, "eps": 1e-5, "U": 0.1, "rho": 1000.0, "mu": 1e-3}   # one pipe state in SI
for system in ("SI", "cgs", "imperial"):                                         # the same state in three unit systems
    nums = ch01.rescale_units(pipe_vals, ch01.PIPE, system)                      # new numbers for the same physics
    print(system, [round(ch01.group_value(g, nums), 6) for g in groups])        # [10.0, 100.0, 0.001, 0.001] each time
extra = {"dp": 1, "d": 2, "rho": 1, "mu": -2}                                    # C71: Δp d² ρ / μ² = Π₁/Π₄²
print("four groups independent:", ch01.groups_independent(groups))              # True
print("with a fifth group:", ch01.groups_independent(groups + [extra]))         # False: only n − r are independent
""", explain="""
1. `pi_groups` finds a basis of the null space with the chosen repeating variables; `group_latex` prints each group.
2. `solve_exponents` solves the 3×3 exponent system of C70 exactly.
3. `rescale_units` + `group_value`: every group keeps its value in SI, cgs and imperial units.
4. `groups_independent` checks the rank of the groups' exponent matrix: four are independent, a fifth never is.
""")
P("P57", "np.linalg.solve", """
`np.linalg.solve(B, b)` finds the vector x with B·x = b for a square, nonsingular matrix B.
""", code="""
print(np.linalg.solve(np.array([[2.0, 0.0], [0.0, 4.0]]), np.array([2.0, 8.0])))   # 2x = 2, 4y = 8 → [1. 2.]
""")
nb.md("#### From scratch: the exponent system of Π₁ with numpy")
nb.check_agree("""
B = A[:, [names.index(k) for k in ("U", "d", "rho")]]         # columns of the repeating variables (a nonsingular minor)
b = -A[:, names.index("dp")]                                  # move the Δp column (exponent 1) to the right-hand side
x = np.linalg.solve(B.astype(float), b.astype(float))         # exponents a, b, c of U, d, ρ (P57)
assert np.allclose(x, [-2, 0, -1])                            # Π₁ = Δp U⁻² d⁰ ρ⁻¹
lib = ch01.solve_exponents("dp", ("U", "d", "rho"), ch01.PIPE)   # the library's exact exponents
assert np.allclose(x, [float(lib[k]) for k in ("U", "d", "rho")])   # same as the library's exact fractions
print("exponents of U, d, ρ in Π₁:", x)   # [−2, 0, −1]
""")
nb.figure("""
rng_p = np.random.default_rng(7)                                                  # seeded pipe states (P10)
lo_hi = {"dp": (1, 1e5), "dx": (0.1, 100), "d": (1e-3, 1), "eps": (1e-6, 1e-3), "U": (1e-3, 10), "rho": (1, 1e3), "mu": (1e-5, 1)}  # SI ranges
states = [{k: 10**rng_p.uniform(np.log10(a), np.log10(b_)) for k, (a, b_) in lo_hi.items()} for _ in range(60)]  # log-uniform
fig, (axg, axr) = plt.subplots(1, 2, figsize=(9, 3.8))                            # left: groups, right: raw variables
for system, col in (("cgs", COLORS["teal"]), ("imperial", COLORS["orange"])):     # two other unit systems
    xs_g, ys_g, xs_r, ys_r = [], [], [], []                                       # collected points
    for s in states:                                                              # every random pipe state
        new = ch01.rescale_units(s, ch01.PIPE, system)                            # the same state in other units
        xs_g += [ch01.group_value(g, s) for g in groups]; ys_g += [ch01.group_value(g, new) for g in groups]   # Π in SI vs other
        xs_r += [s["dp"], s["mu"]]; ys_r += [new["dp"], new["mu"]]                # two raw variables for contrast
    axg.loglog(xs_g, ys_g, "o", ms=3, color=col, label=system)                    # groups
    axr.loglog(xs_r, ys_r, "o", ms=3, color=col, label=system)                    # raw numbers
for ax, title in ((axg, "Π groups: every point on the 1:1 line"), (axr, "raw Δp and μ: shifted by the unit factors")):   # both panels
    lim = [1e-12, 1e12]; ax.plot(lim, lim, color=COLORS["muted"], lw=1); ax.set_xlim(lim); ax.set_ylim(lim)   # 1:1 line, same range
    ax.set_xlabel("value in SI units"); ax.set_title(title, fontsize=10); ax.legend(fontsize=8)              # labels
axg.set_ylabel("value in the other unit system")                                  # shared meaning of y
plt.show()                                                                        # display
""", see="""
Left: the four Π groups of 60 random pipe states, recomputed in cgs (teal) and imperial units (orange), all exactly on the
1:1 line. Right: the raw numbers for Δp and μ of the same states, lying on lines parallel to the diagonal but shifted away
from it.
""", read="""
A group's value does not depend on the units; a dimensional variable's number does (by the factor λ_M^a λ_L^b λ_T^c of P59).
That invariance is what step 12 of D28 uses.
""", change="""
Plot a combination that is *not* dimensionless, such as Δp/(ρU): its points would leave the diagonal by the length-unit
factor, like the raw variables on the right.
""")
nb.figure("""
fig, (ax_raw, ax_pi) = plt.subplots(1, 2, figsize=(9, 3.8))                        # left raw data, right Π coordinates
dx_p = 1.0                                                                          # pipe length [m]
fam = 0                                                                             # family counter, to stagger the Re values
for name, mk in (("water", "o"), ("glycerine", "s"), ("air", "^")):                # three fluids
    fp = ch01.FLUIDS[name]                                                          # ρ [kg/m³], μ [Pa s] from fluidpy's property table
    for d_p, col in zip((2e-3, 5e-3, 1e-2, 2e-2), (COLORS["accent"], COLORS["teal"], COLORS["orange"], COLORS["blue"])):   # four diameters [m]
        Re = np.logspace(0.2, np.log10(2000) - 0.24, 8) * 10**(0.04 * fam - 0.2)     # laminar Reynolds numbers 1 … 2000, staggered
        fam += 1                                                                    # next family gets slightly larger Re
        U_p = Re * fp["mu"] / (fp["rho"] * d_p)                                     # the speeds that give them [m/s]
        dp_p = ch01.poiseuille_pressure_drop(fp["mu"], U_p, dx_p, d_p)              # laminar Δp (result derived in Ch. 8; data generator only)
        ax_raw.loglog(U_p, dp_p / dx_p, mk, ms=4, color=col)                        # raw: pressure gradient vs speed
        ax_pi.loglog(Re, dp_p / (fp["rho"] * U_p**2) * d_p / dx_p, mk, ms=4, color=col)   # Π₁ · d/Δx vs Re
Re_line = np.logspace(0, np.log10(2000), 50)   # Reynolds numbers for the reference line
ax_pi.loglog(Re_line, 32 / Re_line, color=COLORS["muted"], lw=1.5, label="32/Re")  # the laminar law
ax_raw.set_xlabel("mean speed U [m/s]"); ax_raw.set_ylabel("Δp/Δx [Pa/m]")   # raw axes with units
ax_raw.set_title("raw data: 12 separate families", fontsize=10)   # left message
ax_pi.set_xlabel("Re = ρUd/μ = 1/Π₄ [–]"); ax_pi.set_ylabel("Π₁ · d/Δx [–]")   # dimensionless axes
ax_pi.set_title("Π coordinates: one line", fontsize=10)   # right message
from matplotlib.lines import Line2D                                                 # simple legend handles
handles = [Line2D([], [], ls="", marker=m, color=COLORS["ink"], label=n) for n, m in (("water", "o"), ("glycerine", "s"), ("air", "^"))]   # marker = fluid
handles += [Line2D([], [], ls="", marker="o", color=c, label=f"d = {d*1e3:g} mm")     # colour = diameter
            for d, c in zip((2e-3, 5e-3, 1e-2, 2e-2), (COLORS["accent"], COLORS["teal"], COLORS["orange"], COLORS["blue"]))]
ax_raw.legend(handles=handles, fontsize=7.5, ncol=2, loc="upper left")              # marker = fluid, colour = diameter
ax_pi.legend(fontsize=8)                                                            # the 32/Re line
savefig(fig, "ch01", "c69_collapse")                                                # keep a copy for review
plt.show()                                                                          # display
""", see="""
Left: pressure gradient against speed for water (circles), glycerine (squares) and air (triangles) in four pipe diameters
(colours) — twelve families spread over many decades. Right: the same 96 points in Π coordinates, all lying on one grey line.
""", read="""
Once plotted as Π₁·d/Δx against Re, the data depend on a single number — why engineers plot friction factors against the
Reynolds number. Without μ the only groups left would be Δx/d and ε/d, which are the same for every point in one pipe, so
nothing could collapse the spread of Π₁: leaving out a variable is fatal (C65).
""", change="""
Turbulent data from rough pipes would spread out again by ε/d — then Π₃ matters, and one line becomes a family of curves
(the Moody chart).
""")
note("C73", "Example 1.2 — the scale height from dimensions alone (stated)", """
Variables H, T₀, M_w, g, R_u (n = 5); M, L, T and Θ all appear and r = 4, so there is one group, HgM_w/(R_uT₀) — and a lone
group must be a constant: H = const · R_uT₀/(gM_w). The isothermal atmosphere (C63) fixes the constant to 1; our numbers
at 250 K give 7.32 km both ways (code below).
""")
note("C74", "Example 1.3 — Pythagoras", """
Even Pythagoras follows: a right triangle's area is C²φ(β) for hypotenuse C and angle β, and the two smaller similar
triangles (*gloss:* same angles, sides in one ratio) add up to it, giving A² + B² = C² (φ(β) = ¼ sin 2β as a check) — a
curiosity, not used later.
""")
note("C75", "Example 1.4 — the energy of a blast wave (stated)", """
G. I. Taylor estimated an explosion's energy from photographs: E, ρ, the front radius D and the time t give r = 3 and one
group, so E = KρD⁵/t² and D ∝ t^(2/5). K comes only from the full similarity solution (*gloss:* a solution whose shape
stays the same as it grows): K ≈ 0.86 for a **free spherical** blast in air (γ = 1.4). Our number: a front of radius 100 m
after 25 ms in air of 1.2 kg/m³ gives E = 1.9×10¹³ K J ≈ 4.6K kilotons of TNT. A blast at the ground spreads as a hemisphere,
which behaves like half of a free sphere of energy 2E — so the same photograph then means half that energy.
""")
note("C76", "Example 1.5 — why the sky is blue (stated)", r"""
Scattered intensity S from a small particle of volume V at distance d, for light of wavelength λ (*gloss:* the colour of
light; its intensity grows with the square of the wave amplitude). Only L and the intensity dimension appear (the T row is
−3 times the M row, so r = 2); energy spreading forces d⁻², the dipole physics forces V², and then the groups force λ⁻⁴.
Blue light (450 nm) scatters (700/450)⁴ = 5.9 times more than red (700 nm).
""", equation=r"\frac{S}{I} = \frac{V^2}{d^2\lambda^4}\,\varphi_3(n_s)")
nb.code("""
import warnings                                                                      # standard library
with warnings.catch_warnings():                                                      # fluidpy warns that the kmol is treated as a
    warnings.simplefilter("ignore", UserWarning)                                     #   pure number (C66) — expected here
    r_sh = ch01.rank_by_minors(ch01.dimensional_matrix(ch01.SCALE_HEIGHT)[0])[0]     # C73: rank of the scale-height matrix
    g_sh = [ch01.group_latex(g) for g in ch01.pi_groups(ch01.SCALE_HEIGHT)]          # its single group
print("C73: r =", r_sh, "| group:", g_sh)   # rank 4, one group
print(f"C73: R_u T/(g M_w) = {ch01.R_U*250/(ch01.G0*ch01.M_W_AIR):.1f} m; scale_height(250 K) = {ch01.scale_height(250.0):.1f} m")   # dimensional estimate = scale height
print(f"C74: φ(30°) = {ch01.pythagoras_phi(np.radians(30)):.4f}; ¼ sin 60° = {0.25*np.sin(np.radians(60)):.4f}")   # area function check
E_over_K = ch01.blast_energy(100.0, 0.025, 1.2)                                      # C75: E/K for D = 100 m at t = 25 ms [J]
print(f"C75: E = {E_over_K:.2e} K J = {E_over_K/4.184e12:.2f} K kt of TNT; free-sphere K = {ch01.TAYLOR_K_GAMMA14}")   # energy per unit K
ratio = ch01.rayleigh_scattering_ratio(1, 1, 450e-9) / ch01.rayleigh_scattering_ratio(1, 1, 700e-9)   # C76: blue/red
print(f"C76: blue/red scattering = {ratio:.3f}; rank of the Rayleigh matrix = {ch01.rank_by_minors(ch01.dimensional_matrix(ch01.RAYLEIGH)[0])[0]}")   # λ⁻⁴ law
""", explain="""
1. `dimensional_matrix` and `rank_by_minors` for Example 1.2 give r = 4 and one group; the dimensional estimate R_uT/(gM_w)
   equals `scale_height` exactly.
2. `pythagoras_phi` is the area function of Example 1.3.
3. `blast_energy` evaluates E = KρD⁵/t² with K = 1 (so the output is "per unit K"); `TAYLOR_K_GAMMA14` is Taylor's constant
   for a free spherical blast.
4. `rayleigh_scattering_ratio` gives the λ⁻⁴ law's blue-to-red ratio, and the Rayleigh matrix has rank 2.
""")
nb.explainer("buckingham_pi_machine", heading="Why can 7 pipe variables shrink to 4 numbers?", why="""
Toggling variables, swapping the repeating set and switching units changes the matrix, its rank and the groups at once — you
experiment with the method instead of watching it done once. Its Derivation tab walks through D28.
""", tries=[
    "Choose the repeating set (d, ε, ρ) and read why the machine refuses.",
    "Remove μ: how many groups are left, and which physics is lost?",
    "Switch to the blast preset: what power of t does the radius follow?",
    "Change SI → imperial and watch every Π value stay put.",
])
nb.md("""
**What would change if…** the problem included rotation (the Earth's rotation rate Ω) as an eighth variable? One more
column, the same rank, one more group — the Rossby number of Ch. 13.
""")
nb.pointer("""
The chapter's 30 exercises are in the book and are not reproduced here. Of the derivations the book leaves to Exercises 1.10,
1.11, 1.13 and 1.14, D14 (isentropic law), D18 (parcel equation) and D19 (adiabatic lapse rate) are written out above; the
result of Exercise 1.10 (e = e(T) for a perfect gas) is stated in the C40 block. `S01`
""")
nb.pointer("""
The book's reading list is at the end of the printed chapter; the public data used here (CODATA, USSA-1976, IAPWS, Jennings,
Taylor) are cited in reference/ch01/SOURCES.md. `S02`
""")

# =====================================================================================================================
# summary
# =====================================================================================================================
nb.summary(
    clicked=[
        "A density 'at a point' is the plateau of box averages between molecular noise and the flow's own variation (C06).",
        "Viscosity is momentum diffusion: τ = μ du/dy, and ν = μ/ρ sets how fast motion spreads (C12).",
        "In still fluid pressure falls with height at exactly ρg per metre; buoyancy is the resulting face-pressure difference (C20).",
        "Heat plus work changes internal energy; q and w depend on the route, Δe does not (C25).",
        "Entropy is the route-independent sum of dq/T, and the Gibbs relations tie state functions together for any process (C35).",
        "The speed of sound is √((∂p/∂ρ)_s); an incompressible fluid has c → ∞ (C36).",
        "p = ρRT is the molecular gas law per unit mass, with R = R_u/M_w (C40).",
        "With no heat and no friction a perfect gas follows p ∝ ρ^γ (C45).",
        "A displaced parcel obeys ζ″ + N²ζ = 0 (C50).",
        "N² > 0 oscillates, N² = 0 stays, N² < 0 runs away (C51).",
        "Air is stable when it cools more slowly than 9.8 K/km: dT/dz > −9.8 K/km (Kundu) ⇔ Γ_met < +9.8 K/km (meteorology) (C54).",
        "Potential temperature removes adiabatic cooling; stability is θ increasing upward, N² = (g/θ)dθ/dz (C55).",
        "Every term of a correct law has the same dimensions (C64).",
        "The dimensions of all variables form an integer matrix whose rank counts the independent dimensions (C67).",
        "Dimensionless groups are the null space of that matrix: there are n − r of them (C69).",
    ],
    feeds_forward=[
        "Ch. 2–3: fields and fluid particles on the continuum",
        "Ch. 4: the tensor viscosity law, the energy equation, Boussinesq (α), dimensionless Navier–Stokes (Π)",
        "Ch. 7 and 13: N² for internal waves; θ and hydrostatics for geophysical fluid dynamics",
        "Ch. 8: the Couette start-up and Poiseuille results used here as checks",
        "Ch. 15: the speed of sound and the isentropic relations",
    ],
    left_out=[
        "the 30 exercises (in the book)",
        "moist thermodynamics: saturated adiabats and equivalent potential temperature (beyond Ch. 1)",
        "the proof of the speed-of-sound formula (1.19) (Ch. 15 §15.2)",
    ],
)

path = nb.save()
print(f"wrote {path.relative_to(pathlib.Path(__file__).resolve().parents[1]).as_posix()} "
      f"({len(nb.cells)} cells, {len(nb.cores)} CORE blocks, {len(nb.derivations)} derivations, "
      f"{len(nb.primers)} primers, {len(nb.explainers)} explainers)")
