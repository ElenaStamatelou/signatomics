"""Charts for the statistics results, in the same style as every signatomics chart."""

from __future__ import annotations

import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from ..olink.plots import CRITICAL, INK_2, MUTED, SERIES, SURFACE, _style
from ._common import StatsResult


def plot_volcano(result: StatsResult, *, alpha: float = 0.05, labels: int = 10) -> go.Figure:
    """Volcano plot: the effect of every protein (x) against how sure it is (y, -log10 p).

    Works for ``compare_groups`` with two groups (x = difference in NPX, the log2 fold change),
    with an ordered trend (x = slope per group step) and for ``correlate`` (x = r). Proteins with
    an FDR below ``alpha`` are drawn in colour and the strongest ``labels`` of them are named;
    hover any dot for its numbers.
    """
    t = result.table
    if t is None or t.empty:
        raise ValueError("The result has no table to plot.")
    effect, axis = next(((c, a) for c, a in (("difference", "difference (NPX, = log2 fold change)"),
                                             ("slope", "slope (NPX per group step)"),
                                             ("r", "correlation (r)")) if c in t.columns), (None, None))
    if effect is None:
        raise ValueError("A volcano plot needs an effect per protein: compare two groups, test a "
                         "trend (ordered=True), or correlate.")
    t = t.dropna(subset=[effect, "p"]).copy()
    t["_y"] = -np.log10(t["p"].clip(lower=1e-300))
    hit = t["fdr"] < alpha
    up, down = hit & (t[effect] > 0), hit & (t[effect] < 0)

    fig = go.Figure()
    for mask, name, colour in ((~hit, f"FDR {alpha:g} or more", MUTED),
                               (up, f"higher (FDR below {alpha:g})", SERIES[1]),
                               (down, f"lower (FDR below {alpha:g})", SERIES[0])):
        part = t[mask]
        if part.empty:
            continue
        fig.add_scatter(
            x=part[effect], y=part["_y"], mode="markers", name=name,
            marker=dict(size=8, color=colour, opacity=0.85, line=dict(color=SURFACE, width=1)),
            customdata=np.column_stack([part["protein"].astype(str), part["ci_low"], part["ci_high"],
                                        part["p"], part["fdr"]]),
            hovertemplate=(f"<b>%{{customdata[0]}}</b><br>{effect} %{{x:+.2f}} "
                           "(%{customdata[1]:+.2f} to %{customdata[2]:+.2f})"
                           "<br>p %{customdata[3]:.2g}<br>FDR %{customdata[4]:.2g}<extra></extra>"))
    top = t[hit].sort_values("p").head(labels)
    for _, r in top.iterrows():
        fig.add_annotation(x=r[effect], y=r["_y"], text=str(r["protein"]), showarrow=False,
                           yshift=11, font=dict(size=11, color=INK_2))
    fig.add_vline(x=0, line=dict(color=INK_2, width=1, dash="dash"))
    if hit.any():                       # the p-value of the weakest protein that passed the FDR
        fig.add_hline(y=t.loc[hit, "_y"].min(), line=dict(color=CRITICAL, width=1, dash="dot"))
    fig.update_xaxes(title=axis)
    fig.update_yaxes(title="-log10 p (higher = more certain)")
    fig.update_layout(height=520, width=760, margin=dict(l=70, r=30, b=60))
    return _style(fig, f"Volcano plot of the {result.title.lower()}",
                  [f"{len(t)} proteins", f"{int(hit.sum())} with FDR below {alpha:g}",
                   "dotted = the FDR cutoff" if hit.any() else "",
                   "hover a dot for its numbers"])


def plot_volcano_pair(unadjusted: StatsResult, adjusted: StatsResult, *, alpha: float = 0.05,
                      labels: int = 6) -> go.Figure:
    """Two volcano plots side by side on the same axes: unadjusted (left), adjusted (right).

    Both results come from ``compare_groups`` on the same groups; the second is run with
    ``covariates=[...]``. Proteins with an FDR below ``alpha`` are drawn in colour and the
    strongest ``labels`` of each panel are named. A protein that keeps its place in the
    right-hand panel does not owe its difference to the covariates.
    """
    panels = []
    for res in (unadjusted, adjusted):
        t = res.table
        if t is None or t.empty:
            raise ValueError("A result has no table to plot.")
        effect = next((c for c in ("difference", "slope", "r") if c in t.columns), None)
        if effect is None:
            raise ValueError("A volcano plot needs an effect per protein: compare two groups or "
                             "test a trend (ordered=True).")
        panels.append((effect, t.dropna(subset=[effect, "p"]).assign(_y=lambda d: -np.log10(d["p"].clip(lower=1e-300)))))
    if panels[0][0] != panels[1][0]:
        raise ValueError("Both results must answer the same question (the same groups and test).")
    effect = panels[0][0]
    axis = {"difference": "difference (NPX, = log2 fold change)",
            "slope": "slope (NPX per group step)", "r": "correlation (r)"}[effect]

    fig = make_subplots(rows=1, cols=2, shared_yaxes=True, horizontal_spacing=0.04,
                        subplot_titles=("Unadjusted", "Adjusted"))
    seen = set()
    for col, (_, t) in enumerate(panels, start=1):
        hit = t["fdr"] < alpha
        up, down = hit & (t[effect] > 0), hit & (t[effect] < 0)
        for mask, name, colour in ((~hit, f"FDR {alpha:g} or more", MUTED),
                                   (up, f"higher (FDR below {alpha:g})", SERIES[1]),
                                   (down, f"lower (FDR below {alpha:g})", SERIES[0])):
            part = t[mask]
            if part.empty:
                continue
            fig.add_scatter(
                x=part[effect], y=part["_y"], mode="markers", name=name, legendgroup=name,
                showlegend=name not in seen,
                marker=dict(size=8, color=colour, opacity=0.85, line=dict(color=SURFACE, width=1)),
                customdata=np.column_stack([part["protein"].astype(str), part["ci_low"], part["ci_high"],
                                            part["p"], part["fdr"]]),
                hovertemplate=(f"<b>%{{customdata[0]}}</b><br>{effect} %{{x:+.2f}} "
                               "(%{customdata[1]:+.2f} to %{customdata[2]:+.2f})"
                               "<br>p %{customdata[3]:.2g}<br>FDR %{customdata[4]:.2g}<extra></extra>"),
                row=1, col=col)
            seen.add(name)
        for _, r in t[hit].sort_values("p").head(labels).iterrows():
            fig.add_annotation(x=r[effect], y=r["_y"], text=str(r["protein"]), showarrow=False,
                               yshift=11, font=dict(size=11, color=INK_2), row=1, col=col)
        fig.add_vline(x=0, line=dict(color=INK_2, width=1, dash="dash"), row=1, col=col)
        if hit.any():
            fig.add_hline(y=t.loc[hit, "_y"].min(), line=dict(color=CRITICAL, width=1, dash="dot"),
                          row=1, col=col)
        fig.update_xaxes(title=axis, row=1, col=col)
    fig.update_yaxes(title="-log10 p (higher = more certain)", row=1, col=1)
    xs = np.concatenate([p[effect].to_numpy() for _, p in panels])
    lim = float(np.nanmax(np.abs(xs))) * 1.08
    fig.update_xaxes(range=[-lim, lim])
    n_hit = [int((p["fdr"] < alpha).sum()) for _, p in panels]
    fig.update_layout(height=520, width=980, margin=dict(l=70, r=30, b=60))
    fig.for_each_annotation(lambda a: a.update(font=dict(size=13, color=INK_2)) if a.text in ("Unadjusted", "Adjusted") else None)
    return _style(fig, "Volcano plots before and after adjustment",
                  [f"FDR below {alpha:g}: {n_hit[0]} unadjusted, {n_hit[1]} adjusted",
                   "dotted = the FDR cutoff", "hover a dot for its numbers"])
