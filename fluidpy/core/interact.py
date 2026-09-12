"""Python interactive plots that keep working in the PUBLISHED page (not only while a kernel is running).

Two kinds of interactivity, and when to use which
-------------------------------------------------

=========================  ==============================  ===========================================================
tool                       works in                        use it for
=========================  ==============================  ===========================================================
:func:`slider_figure`      Jupyter, VS Code, Colab AND     "move one parameter and watch the curve(s) change": every
(plotly, precomputed)      the static HTML page            slider position is computed in Python up front, so the
                                                           page needs no kernel. Hover shows exact values.
:func:`animate_figure`     same                            a time evolution with play/pause and a time slider
:func:`live`               a running kernel only           exploring 2+ parameters freely (ipywidgets). In the
(ipywidgets)                                               published page it is replaced by a note that points to the
                                                           Colab button, so always pair it with a slider_figure.
=========================  ==============================  ===========================================================

Example — dispersion curves with a depth slider::

    import numpy as np
    from fluidpy.core.interact import slider_figure

    k = np.linspace(0.01, 2, 300)
    def curves(H):                                   # return {trace name: (x, y)}
        return {"c(k)": (k, np.sqrt(9.81 / k * np.tanh(k * H)))}
    fig = slider_figure(curves, "H", np.linspace(1, 50, 25), unit="m",
                        xlabel="wavenumber k [1/m]", ylabel="phase speed c [m/s]", title="Deeper water, faster waves")
    fig.show()

Keep figures light: ≤ 40 slider steps × ≤ 4 traces × ≤ 400 points is ~ 200 kB in the page.
"""
from __future__ import annotations

from typing import Callable, Iterable, Mapping, Sequence

import numpy as np

from .style import COLORS, CYCLE

LAYOUT = dict(
    template="plotly_white",
    font=dict(family="Inter, -apple-system, Segoe UI, Roboto, sans-serif", size=13, color=COLORS["ink"]),
    margin=dict(l=60, r=20, t=60, b=90),
    hovermode="x unified",
    legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="left", x=0),
    colorway=CYCLE,
)


def _fmt(v: float) -> str:
    return f"{v:.3g}"


def slider_figure(fn: Callable[[float], Mapping[str, tuple[Sequence[float], Sequence[float]]]], name: str,
                  values: Iterable[float], *, unit: str = "", xlabel: str = "", ylabel: str = "", title: str = "",
                  xrange: Sequence[float] | None = None, yrange: Sequence[float] | None = None, active: int | None = None,
                  height: int = 460, modes: Mapping[str, str] | None = None):
    """A plotly figure whose curves follow a slider over ``values`` of the parameter ``name``.

    Parameters
    ----------
    fn : function ``fn(value) -> {trace_name: (x, y)}`` — called once per slider value, in Python, now.
    name, unit : the parameter shown on the slider ("H", "m").
    values : the slider positions (≤ ~40).
    xrange, yrange : fix the axes so the reader sees the change, not the autoscale. Default: the global min/max.
    modes : optional ``{trace_name: "lines" | "markers" | "lines+markers"}``.
    """
    import plotly.graph_objects as go

    values = [float(v) for v in values]
    data = [fn(v) for v in values]
    names = list(data[0].keys())
    n = len(names)
    active = len(values) // 2 if active is None else active
    fig = go.Figure()
    for i, v in enumerate(values):
        for j, tname in enumerate(names):
            x, y = data[i][tname]
            fig.add_trace(go.Scatter(x=np.asarray(x), y=np.asarray(y), name=tname, visible=(i == active),
                                     mode=(modes or {}).get(tname, "lines"),
                                     line=dict(color=CYCLE[j % len(CYCLE)], width=2.5),
                                     hovertemplate=f"{tname}: %{{y:.4g}}<extra>{name} = {_fmt(v)} {unit}</extra>"))
    if xrange is None:
        xs = np.concatenate([np.asarray(d[t][0], float) for d in data for t in names])
        xrange = [float(np.nanmin(xs)), float(np.nanmax(xs))]
    if yrange is None:
        ys = np.concatenate([np.asarray(d[t][1], float) for d in data for t in names])
        lo, hi = float(np.nanmin(ys)), float(np.nanmax(ys))
        pad = 0.05 * (hi - lo or 1.0)
        yrange = [lo - pad, hi + pad]
    steps = []
    for i, v in enumerate(values):
        vis = [False] * (len(values) * n)
        vis[i * n:(i + 1) * n] = [True] * n
        steps.append(dict(method="update", args=[{"visible": vis}], label=_fmt(v)))
    fig.update_layout(**LAYOUT, height=height, title=dict(text=title, x=0),
                      xaxis=dict(title=xlabel, range=list(xrange)), yaxis=dict(title=ylabel, range=list(yrange)),
                      sliders=[dict(active=active, steps=steps, pad=dict(t=40),
                                    currentvalue=dict(prefix=f"{name} = ", suffix=f" {unit}", font=dict(size=14)))])
    return fig


def animate_figure(frame_fn: Callable[[float], Mapping[str, tuple[Sequence[float], Sequence[float]]]],
                   times: Iterable[float], *, time_label: str = "t", unit: str = "s", xlabel: str = "", ylabel: str = "",
                   title: str = "", xrange: Sequence[float] | None = None, yrange: Sequence[float] | None = None,
                   frame_ms: int = 60, height: int = 460):
    """A plotly animation (play / pause / time slider) of ``frame_fn(t) -> {trace_name: (x, y)}`` over ``times``.

    Axes are fixed from the global data range unless given — an autoscaling animation hides what moves.
    """
    import plotly.graph_objects as go

    times = [float(t) for t in times]
    data = [frame_fn(t) for t in times]
    names = list(data[0].keys())
    xs = np.concatenate([np.asarray(d[t][0], float) for d in data for t in names])
    ys = np.concatenate([np.asarray(d[t][1], float) for d in data for t in names])
    xrange = xrange or [float(np.nanmin(xs)), float(np.nanmax(xs))]
    if yrange is None:
        lo, hi = float(np.nanmin(ys)), float(np.nanmax(ys))
        pad = 0.05 * (hi - lo or 1.0)
        yrange = [lo - pad, hi + pad]

    def traces(d):
        return [go.Scatter(x=np.asarray(d[t][0]), y=np.asarray(d[t][1]), name=t, mode="lines",
                           line=dict(color=CYCLE[j % len(CYCLE)], width=2.5)) for j, t in enumerate(names)]

    frames = [go.Frame(data=traces(d), name=_fmt(t)) for d, t in zip(data, times)]
    fig = go.Figure(data=traces(data[0]), frames=frames)
    play = dict(label="▶ Play", method="animate",
                args=[None, dict(frame=dict(duration=frame_ms, redraw=False), transition=dict(duration=0), fromcurrent=True)])
    pause = dict(label="⏸ Pause", method="animate",
                 args=[[None], dict(frame=dict(duration=0, redraw=False), mode="immediate", transition=dict(duration=0))])
    fig.update_layout(
        **LAYOUT, height=height, title=dict(text=title, x=0),
        xaxis=dict(title=xlabel, range=list(xrange)), yaxis=dict(title=ylabel, range=list(yrange)),
        updatemenus=[dict(type="buttons", direction="left", x=0, y=-0.28, xanchor="left", yanchor="top",
                          showactive=False, buttons=[play, pause])],
        sliders=[dict(active=0, x=0.18, len=0.82, y=-0.22, pad=dict(t=10),
                      currentvalue=dict(prefix=f"{time_label} = ", suffix=f" {unit}"),
                      steps=[dict(method="animate", label=_fmt(t),
                                  args=[[_fmt(t)], dict(mode="immediate", frame=dict(duration=0, redraw=False),
                                                         transition=dict(duration=0))]) for t in times])])
    return fig


def live(fn: Callable[..., object], **widgets):
    """``ipywidgets.interact`` with sensible slider defaults — **only interactive while a kernel runs**.

    Tag the cell ``live-only`` (``tools/nbkit.py`` does it) so the published page replaces the frozen output with a
    note pointing to Colab. Pass widgets as ``name=(min, max, step)`` tuples or ipywidgets objects.
    """
    import ipywidgets as w

    controls = {}
    for k, spec in widgets.items():
        if isinstance(spec, tuple) and len(spec) == 3:
            lo, hi, step = spec
            controls[k] = w.FloatSlider(min=lo, max=hi, step=step, value=(lo + hi) / 2, continuous_update=False,
                                        description=k, style={"description_width": "initial"},
                                        layout=w.Layout(width="95%"))
        else:
            controls[k] = spec
    return w.interact(fn, **controls)


__all__ = ["slider_figure", "animate_figure", "live", "LAYOUT"]
