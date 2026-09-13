---
name: python-viz
description: Figures, animations and Python-side interactive plots for fluidpy notebooks - the house matplotlib style, when to use a static figure vs a matplotlib animation (video or frame player) vs a plotly slider/animation figure (works on the published page) vs ipywidgets (kernel only), size budgets, and the exact helper calls in fluidpy.core.style, fluidpy.core.anim and fluidpy.core.interact. Load when writing any figure, animation or interactive cell.
---

# python-viz — pictures that teach, in every host

## 1. Choose the right kind
| The reader needs to… | Use | Works on the web page? |
|---|---|---|
| see one relationship | static matplotlib figure + What you see / How to read it / What would change if | yes |
| see how a curve changes with one parameter | `slider_figure` (plotly; every slider position precomputed) | **yes** |
| watch something evolve in time (surface, particles, profile diffusing) | `animate` + `show_animation(player="video")` | yes |
| step through discrete stages (time steps of a scheme, iterations) | `show_animation(player="frames")` or `animate_figure` (plotly) | yes |
| explore 2+ parameters freely | `live(...)` ipywidgets, **plus** a slider_figure for the page | kernel only (page shows a note) |
| manipulate a phenomenon with a guided story | an HTML explainer (skill `interactive-viz`) — 4–5 per chapter | yes |

## 2. House style (`fluidpy.core.style`)
- The setup cell calls `setup_notebook()` → `use_style()` (palette shared with the explainers: purple `#6c5ce7`, teal,
  orange, blue, rose), plotly renderer (Colab: `colab`; elsewhere `plotly_mimetype+notebook_connected` so the exported
  page loads plotly from the CDN), numpy print options, and returns `FAST`.
- Axis labels always carry units: `ax.set_xlabel("wavenumber $k$ [m$^{-1}$]")`. Titles say the message, not the variable
  ("Short waves are slower in deep water"). Legends name parameter values ("$H$ = 10 m").
- Non-dimensional plots say so: "$c/\sqrt{gH}$ [–]". Mark regimes/limits with light dashed lines and a text label.
- Save figures you want to reuse: `savefig(fig, "ch07", "dispersion")` → `outputs/ch07/dispersion.png` (git-ignored).
- One idea per figure; two panels only when comparing (`fig, (a, b) = plt.subplots(1, 2, figsize=(9, 3.6))`).

## 3. Animations (`fluidpy.core.anim`)
```python
from fluidpy.core.anim import animate, show_animation
fig, ax = plt.subplots(figsize=(6.5, 3.2))              # keep ≤ 7×4 in
(surface,) = ax.plot(x, eta(x, 0))                      # draw the first frame once
dots = ax.scatter(xp, zp, s=8)                          # particles
def update(i):                                          # called for every frame i
    t = i * dt
    surface.set_ydata(eta(x, t))                        # move the free surface
    dots.set_offsets(np.c_[xp + dx(t), zp + dz(t)])     # move each particle on its orbit
    return surface, dots
show_animation(animate(update, frames=60 if not FAST else 30, fig=fig, interval=50))   # MP4 (ffmpeg) or frame player
```
- `player="video"` (default when ffmpeg exists: imageio-ffmpeg locally, system ffmpeg on Colab) is ~10× smaller;
  `player="frames"` gives ◀ ▮▮ ▶ and single-step buttons — use it when stepping matters.
- ≤ 120 frames, dpi 80, fixed axis limits (never autoscale per frame), one animation idea per cell.
- Always follow with What you see / How to read it / What would change if.

## 4. Plotly figures that survive publishing (`fluidpy.core.interact`)
```python
from fluidpy.core.interact import slider_figure, animate_figure
fig = slider_figure(lambda H: {"deep limit": (k, np.sqrt(g / k) + 0*k),
                               "finite depth": (k, ch07.phase_speed(k, H))},
                    "H", np.linspace(0.5, 50, 30), unit="m",
                    xlabel="k [1/m]", ylabel="c [m/s]", title="Depth only matters for long waves")
fig.show()
```
- ≤ 40 slider steps × ≤ 4 traces × ≤ 400 points (~200 kB). Fix axis ranges (the helper does, from the global data).
- `animate_figure(frame_fn, times, …)` gives Play/Pause + a time slider without a kernel.
- 3-D surfaces: `go.Surface` with a dropdown of precomputed cases, `height=520`.

## 5. Live widgets (kernel only)
`nb.live(src)` writes the "Live cell" note and tags the cell `live-only`; inside use `fluidpy.core.interact.live(fn,
a=(0, 1, 0.05), …)` with `continuous_update=False` sliders. The published page replaces the frozen output with a
"run it in Colab" note, so never make a live cell the only way to see an idea.

## 6. Budgets
Whole notebook < 5 min on a laptop (use `FAST`), page < 15 MB, each animation < 6 MB (the helper warns), no base64
data blobs other than figures/animations.

## 7. Lessons (the knowledge-keeper appends; one line each, dated by chapter)
- (ch01) A `live` widget is always paired with a static figure or a `slider_figure` carrying the same idea (C51); the
  published page freezes the widget.
- (ch01) When a quantity has two sign conventions, put both values and both inequalities in the slider trace name or
  legend ("dT/dz = −6.5 > −9.8 K/km ⇔ Γ_met = 6.5 < 9.8 K/km"), so they move together (C54).
- (ch01) Look at every saved figure image before the review: arrows crossing long labels (book map), a label under a
  curve (C36 "c → ∞"), unexplained point letters and a data range that contradicts its comment were all round-2 findings.
- (ch01) `player="frames"` for processes the reader should stop at (zoom by decades, leg-by-leg routes); `"video"` for
  smooth transients (diffusion, parcels); keep `dpi` ≤ 80 to stay inside the page budget.
