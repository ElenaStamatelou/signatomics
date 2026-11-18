"""Charts for the signature functions, in the same style as every signatomics chart."""

from __future__ import annotations

import numpy as np
import plotly.graph_objects as go
from sklearn.metrics import roc_curve

from ..olink.plots import INK_2, MUTED, SERIES, _style
from ._engine import model_label
from .checks import CheckResult
from .find import SignatureResult

def _rgba(hex_colour: str, alpha: float) -> str:
    h = hex_colour.lstrip("#")
    return f"rgba({int(h[0:2], 16)},{int(h[2:4], 16)},{int(h[4:6], 16)},{alpha})"


_MODEL_COLOR = {"regression": SERIES[0], "random_forest": SERIES[1]}


def plot_signature_performance(result: SignatureResult) -> go.Figure:
    """How well the signature predicts people it has not seen.

    Binary outcome: the ROC curve of every model, one thin line per repeat of the cross-validation.
    A number: observed against predicted for the first repeat. Predictions come from the folds in
    which the person was left out.
    """
    pr = result.predictions
    fig = go.Figure()
    if result.task == "binary":
        pos = result.classes[-1]
        grid = np.linspace(0, 1, 101)
        for model, part in pr.groupby("model", sort=False):
            colour = _MODEL_COLOR.get(model, MUTED)
            label = model_label(model, result.task)
            curves = []
            for _, g in part.groupby("repeat"):
                fpr, tpr, _ = roc_curve((g["observed"] == pos).astype(int), g["predicted"])
                curves.append(np.interp(grid, fpr, tpr))
            curves = np.array(curves)
            curves[:, 0] = 0.0
            fig.add_scatter(x=np.r_[grid, grid[::-1]], y=np.r_[curves.max(axis=0), curves.min(axis=0)[::-1]],
                            fill="toself", fillcolor=_rgba(colour, 0.15), line=dict(width=0),
                            legendgroup=model, showlegend=False, hoverinfo="skip")
            fig.add_scatter(x=grid, y=curves.mean(axis=0), mode="lines", line=dict(color=colour, width=2.5),
                            name=label, legendgroup=model,
                            hovertemplate=f"{label}<br>1 - specificity %{{x:.2f}}<br>sensitivity %{{y:.2f}}<extra></extra>")
        fig.add_scatter(x=[0, 1], y=[0, 1], mode="lines", line=dict(color=INK_2, width=1, dash="dash"),
                        name="no information", hoverinfo="skip")
        fig.update_xaxes(title="1 - specificity", range=[0, 1])
        fig.update_yaxes(title="sensitivity", range=[0, 1.02])
        auc = result.performance.query("metric == 'auc'")
        sub = [f"{model_label(r['model'], result.task)} AUC {r['mean']:.2f} ({r['ci_low']:.2f} to {r['ci_high']:.2f})"
               for _, r in auc.iterrows() if "covariates only" not in r["model"]]
        sub.append("line = mean over the repeats, band = their range")
        title = f"ROC curve for {result.classes[-1]}"
    elif result.task == "continuous":
        for model, part in pr[pr["repeat"] == 0].groupby("model", sort=False):
            fig.add_scatter(x=part["observed"], y=part["predicted"], mode="markers",
                            name=model_label(model, result.task),
                            marker=dict(size=7, color=_MODEL_COLOR.get(model, MUTED), opacity=0.7))
        lo, hi = float(pr["observed"].min()), float(pr["observed"].max())
        fig.add_scatter(x=[lo, hi], y=[lo, hi], mode="lines", line=dict(color=INK_2, width=1, dash="dash"),
                        name="perfect", hoverinfo="skip")
        fig.update_xaxes(title=f"observed {result.outcome}")
        fig.update_yaxes(title=f"predicted {result.outcome}")
        rho = result.performance.query("metric == 'spearman'")
        sub = [f"{model_label(r['model'], result.task)} Spearman r {r['mean']:.2f}" for _, r in rho.iterrows() if "covariates only" not in r["model"]]
        title = f"Predicted against observed {result.outcome}"
    else:
        raise ValueError("This chart is for a binary or a continuous outcome. For several classes read "
                         "result.performance and result.predictions.")
    fig.update_layout(height=470, width=620, margin=dict(l=70, r=30, b=60))
    return _style(fig, title, sub)


def plot_selection_frequency(result: SignatureResult, *, top: int = 20) -> go.Figure:
    """How often each protein was chosen: the share of the outer folds, strongest first."""
    p = result.proteins.head(top).iloc[::-1]
    if p.empty:
        raise ValueError("No protein was chosen in any fold.")
    cut = p["in_signature"]
    fig = go.Figure()
    for mask, name, colour in ((cut, "in the signature", SERIES[0]), (~cut, "chosen less often", MUTED)):
        part = p[mask]
        if part.empty:
            continue
        fig.add_scatter(x=part["frequency"], y=part["protein"], mode="markers", name=name,
                        marker=dict(size=10, color=colour),
                        hovertemplate="<b>%{y}</b><br>chosen in %{x:.0%} of the folds<extra></extra>")
    fig.update_xaxes(title="share of the folds that chose the protein", range=[0, 1.03], tickformat=".0%")
    fig.update_layout(height=max(320, 22 * len(p) + 180), width=620, margin=dict(l=90, r=30, b=60))
    return _style(fig, "How stable the signature is",
                  [f"{len(result.signature)} proteins in the signature", "hover a dot for its share"])


def plot_curve(check: CheckResult) -> go.Figure:
    """The learning curve or the selection curve: performance with its 95% interval."""
    t = check.table
    if "fraction" in t.columns:
        x, xt, title = t["n_train"], "people in the training part of each fold", "Learning curve"
    elif "proteins" in t.columns:
        x, xt, title = t["proteins"], "proteins in the signature", "Selection curve"
    else:
        raise ValueError("plot_curve draws the result of learning_curve or selection_curve.")
    name = "AUC" if check.metric == "auc" else "Spearman r"
    fig = go.Figure()
    fig.add_scatter(x=list(x) + list(x)[::-1], y=list(t["ci_high"]) + list(t["ci_low"])[::-1], fill="toself",
                    fillcolor="rgba(120,120,120,0.15)", line=dict(width=0), hoverinfo="skip", name="95% interval")
    fig.add_scatter(x=x, y=t["mean"], mode="lines+markers", name=name, line=dict(color=SERIES[0]),
                    marker=dict(size=9, color=SERIES[0]),
                    hovertemplate=f"%{{x}}<br>{name} %{{y:.2f}}<extra></extra>")
    fig.update_xaxes(title=xt)
    fig.update_yaxes(title=name)
    fig.update_layout(height=430, width=620, margin=dict(l=70, r=30, b=60))
    return _style(fig, title, ["every point repeats the whole nested cross-validation", "hover a point for its value"])


def plot_contributions(check: CheckResult) -> go.Figure:
    """What each protein adds: its value alone and the panel without it, strongest first."""
    t = check.table
    if "alone" not in t.columns:
        raise ValueError("plot_contributions draws the result of protein_contributions.")
    t = t.iloc[::-1]
    name = "AUC" if check.metric == "auc" else "Spearman r"
    fig = go.Figure()
    fig.add_scatter(x=t["alone"], y=t["protein"], mode="markers", name="protein alone",
                    marker=dict(size=10, color=SERIES[1]),
                    hovertemplate=f"<b>%{{y}}</b><br>alone: {name} %{{x:.2f}}<extra></extra>")
    fig.add_scatter(x=t["without"], y=t["protein"], mode="markers", name="panel without the protein",
                    marker=dict(size=10, color=SERIES[0]),
                    hovertemplate=f"<b>%{{y}}</b><br>without it: {name} %{{x:.2f}}<extra></extra>")
    if hasattr(check, "full"):
        fig.add_vline(x=check.full, line=dict(color=INK_2, width=1, dash="dash"))
    fig.update_xaxes(title=name)
    fig.update_layout(height=max(320, 26 * len(t) + 180), width=620, margin=dict(l=90, r=30, b=60))
    return _style(fig, "What each protein adds",
                  ["dashed = the full panel"])
