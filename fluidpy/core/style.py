"""One look for every figure in the project, and the notebook set-up every chapter notebook calls first.

    from fluidpy.core.style import setup_notebook, COLORS, savefig
    FAST = setup_notebook()          # style + plotly renderer + print options; returns the FAST flag

The palette matches the interactive explainers (``assets/viz_base.css``) and the site theme
(``assets/clean-educational.css``), so a matplotlib figure, a plotly figure and an explainer sitting next to each
other in a chapter page read as one family.
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

COLORS = {
    "accent": "#6c5ce7",   # purple — the main curve / the quantity being taught
    "teal": "#0fa39c",     # second quantity, "key idea", positive values
    "orange": "#f97316",   # highlight, the point to watch
    "rose": "#e11d6a",     # warnings, negative values, backflow
    "blue": "#2563eb",     # particles, water
    "amber": "#d97706",
    "ink": "#1a1f36",
    "muted": "#697386",
    "grid": "#e3e6ee",
}
CYCLE = [COLORS[k] for k in ("accent", "teal", "orange", "blue", "rose", "amber")]


def use_style(dpi: int = 110) -> None:
    """Apply the house matplotlib style (labels with units, light grid, colour cycle shared with the explainers)."""
    import matplotlib as mpl
    from cycler import cycler

    mpl.rcParams.update({
        "figure.dpi": dpi,
        "savefig.dpi": 150,
        "figure.figsize": (7.0, 4.2),
        "figure.constrained_layout.use": True,
        "axes.prop_cycle": cycler(color=CYCLE),
        "axes.grid": True,
        "grid.color": COLORS["grid"],
        "grid.linewidth": 0.8,
        "axes.edgecolor": "#c9cedb",
        "axes.labelcolor": COLORS["ink"],
        "axes.titleweight": "semibold",
        "axes.titlesize": 12,
        "axes.labelsize": 11,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "xtick.color": COLORS["muted"],
        "ytick.color": COLORS["muted"],
        "legend.frameon": False,
        "lines.linewidth": 2.0,
        "font.size": 10.5,
        "mathtext.fontset": "cm",
        "image.cmap": "viridis",
        "animation.html": "jshtml",
    })


def plotly_renderer() -> str:
    """Pick the plotly renderer that works in the current host (and survives nbconvert's HTML export)."""
    import plotly.io as pio

    if "google.colab" in sys.modules:
        pio.renderers.default = "colab"
    else:
        # mimetype for JupyterLab/VS Code + a self-loading HTML copy (plotly.js from the CDN) for the exported page
        pio.renderers.default = "plotly_mimetype+notebook_connected"
    return pio.renderers.default


def setup_notebook(fast: bool | None = None, dpi: int = 110) -> bool:
    """Everything a chapter notebook needs before its first figure. Returns the FAST flag.

    ``FAST`` shrinks grids and frame counts so the published run stays quick; it defaults to the environment variable
    ``FLUIDPY_FAST`` ("1" → True) and to False otherwise.
    """
    import numpy as np

    use_style(dpi)
    try:
        plotly_renderer()
    except ImportError:
        pass
    np.set_printoptions(precision=4, suppress=True, linewidth=110)
    if fast is None:
        fast = os.environ.get("FLUIDPY_FAST", "0") == "1"
    return bool(fast)


def savefig(fig, chapter: str, name: str, root: str | Path | None = None) -> Path:
    """Save ``fig`` to ``outputs/<chapter>/<name>.png`` (git-ignored) and return the path."""
    from .project import repo_root

    out = (Path(root) if root else repo_root()) / "outputs" / chapter
    out.mkdir(parents=True, exist_ok=True)
    path = out / (name if name.endswith(".png") else f"{name}.png")
    fig.savefig(path, bbox_inches="tight")
    return path


__all__ = ["COLORS", "CYCLE", "use_style", "plotly_renderer", "setup_notebook", "savefig"]
