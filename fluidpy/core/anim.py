"""Matplotlib animations that play everywhere: Jupyter, VS Code, Google Colab and the published HTML page.

Two players, one call::

    import numpy as np, matplotlib.pyplot as plt
    from fluidpy.core.anim import animate, show_animation

    x = np.linspace(0, 2*np.pi, 200)
    fig, ax = plt.subplots(figsize=(6, 3))
    (line,) = ax.plot(x, np.sin(x))

    def update(i):                          # called once per frame with the frame number
        line.set_ydata(np.sin(x - 0.1*i))   # move the wave to the right
        return (line,)

    show_animation(animate(update, frames=60, fig=fig, interval=50))              # MP4 video if ffmpeg exists
    show_animation(animate(update, frames=40, fig=fig), player="frames")          # frame-by-frame JS player

=============  ==========================================================================================
player         when to use
=============  ==========================================================================================
``"video"``    default when ffmpeg is available (``imageio-ffmpeg`` locally, system ffmpeg on Colab): an
               H.264 ``<video>`` with a scrub bar, ~10× smaller than PNG frames — best for smooth motion.
``"frames"``   ``FuncAnimation.to_jshtml()``: every frame is a PNG plus a JS player with ◀ ▮▮ ▶ and
               single-step buttons — best when the reader should step through time ("what happens at step 3?").
``"auto"``     video if ffmpeg is found, otherwise frames.
=============  ==========================================================================================

Both are self-contained HTML, so they survive nbconvert's export and run inside Colab's sandboxed output frame.
Keep ``frames`` ≲ 120, figures ≲ 7×4 in and ``dpi`` ≲ 90 (frames player) — the helper warns above ``warn_mb``.
"""
from __future__ import annotations

import shutil
import warnings
from functools import lru_cache
from pathlib import Path
from typing import Callable, Iterable

import matplotlib

MAX_FRAMES = 150


@lru_cache(maxsize=1)
def ffmpeg_path() -> str | None:
    """Path to an ffmpeg executable (system PATH, else the one bundled by ``imageio-ffmpeg``), or None."""
    exe = shutil.which("ffmpeg")
    if exe is None:
        try:
            import imageio_ffmpeg

            exe = imageio_ffmpeg.get_ffmpeg_exe()
        except Exception:  # noqa: BLE001
            exe = None
    if exe:
        matplotlib.rcParams["animation.ffmpeg_path"] = exe
    return exe


def animate(update: Callable[[int], object], frames: int | Iterable = 60, fig=None, interval: int = 50,
            blit: bool = False, init: Callable[[], object] | None = None, repeat: bool = True):
    """Build a ``FuncAnimation`` from an ``update(frame)`` function.

    Parameters
    ----------
    update : function called with each frame value; it changes the artists already drawn on ``fig``.
    frames : number of frames (``range(frames)``) or an iterable of frame values.
    fig : the figure to animate (default: the current figure).
    interval : delay between frames in milliseconds (50 ms = 20 frames per second).
    blit : only redraw changed artists (``update`` must then return them); off by default because it is fragile.
    init : optional function drawing the first frame.
    """
    import matplotlib.pyplot as plt
    from matplotlib.animation import FuncAnimation

    fig = fig if fig is not None else plt.gcf()
    return FuncAnimation(fig, update, frames=frames, init_func=init, interval=interval, blit=blit, repeat=repeat,
                         cache_frame_data=False)


def animation_html(anim, *, player: str = "auto", max_frames: int = MAX_FRAMES, dpi: int | None = 80,
                   default_mode: str = "loop", warn_mb: float = 6.0) -> str:
    """Render a ``FuncAnimation`` to self-contained HTML (see the module docstring) and close its figure."""
    import matplotlib.pyplot as plt

    fig = anim._fig  # noqa: SLF001
    if dpi is not None:
        fig.set_dpi(dpi)
    if hasattr(anim, "_save_count") and anim._save_count and anim._save_count > max_frames:  # noqa: SLF001
        warnings.warn(f"animation has {anim._save_count} frames; embedding the first {max_frames}", stacklevel=2)  # noqa: SLF001
        anim._save_count = max_frames  # noqa: SLF001
    if player == "auto":
        player = "video" if ffmpeg_path() else "frames"
    if player == "video" and not ffmpeg_path():
        warnings.warn("ffmpeg not found — falling back to the frames player", stacklevel=2)
        player = "frames"
    old_limit = matplotlib.rcParams["animation.embed_limit"]
    matplotlib.rcParams["animation.embed_limit"] = max(old_limit, 60.0)  # MB
    try:
        if player == "video":
            out = anim.to_html5_video()
            out = out.replace("<video ", '<video style="max-width:100%;height:auto" playsinline muted ', 1)
        else:
            out = anim.to_jshtml(default_mode=default_mode)
    finally:
        matplotlib.rcParams["animation.embed_limit"] = old_limit
        plt.close(fig)  # otherwise Jupyter also shows a duplicate static last frame
    size_mb = len(out.encode("utf-8")) / 1e6
    if size_mb > warn_mb:
        warnings.warn(f"animation is {size_mb:.1f} MB in the page — reduce frames, figure size or dpi", stacklevel=2)
    return out


def show_animation(anim, **kwargs) -> None:
    """Display a ``FuncAnimation`` inline (``player="auto" | "video" | "frames"``; see :func:`animation_html`)."""
    from IPython.display import HTML, display

    display(HTML(animation_html(anim, **kwargs)))


def save_gif(anim, path: str | Path, fps: int = 20, dpi: int = 80) -> Path:
    """Save the animation as an animated GIF with Pillow (no ffmpeg). Use for README/figures, not for notebooks."""
    from matplotlib.animation import PillowWriter

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    anim.save(str(path), writer=PillowWriter(fps=fps), dpi=dpi)
    return path


__all__ = ["animate", "animation_html", "show_animation", "save_gif", "ffmpeg_path"]
