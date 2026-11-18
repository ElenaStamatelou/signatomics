"""Checks on a signature: what each protein adds, how performance grows with more people, and how
many proteins are enough. Each takes the result of ``find_signature`` and reruns the folds."""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pandas as pd

from ._engine import corrected_ci, model_key, model_label, run_cv
from .find import SignatureResult
from .._report import format_report


@dataclass
class CheckResult:
    """Result of a check: summary lines plus a table."""

    title: str = ""
    lines: list = field(default_factory=list)
    table: pd.DataFrame | None = None
    metric: str = ""

    def __str__(self) -> str:
        return format_report(self.title, self.lines)

    def __repr__(self) -> str:       # notebooks show the same text as print()
        return str(self)


def _inputs(result: SignatureResult) -> dict:
    if not result._inputs:
        raise ValueError("This result has no data attached; run find_signature first.")
    return result._inputs


def _metric(task: str) -> str:
    return "spearman" if task == "continuous" else "auc"


def _mean_metric(folds, model, metric):
    vals = [f.metrics[model][metric] for f in folds]
    lo, hi = corrected_ci(vals, np.mean([f.n_train for f in folds]), np.mean([f.n_test for f in folds]),
                              metric=metric)
    return float(np.nanmean(vals)), lo, hi


def _n_prot(k: int) -> str:
    return f"{k} protein{'s' if k != 1 else ''}"


def _classes(i):
    return i["classes"] if i["task"] == "multiclass" else None


def protein_contributions(result: SignatureResult, *, model: str = "regression", verbose: bool = True) -> CheckResult:
    """What each signature protein adds: the panel without it, and the protein on its own.

    The folds are rerun with the signature proteins fixed (the same proteins in every fold). The
    numbers compare the proteins with each other. Because the panel was found with all people, the
    full-panel value is higher than the honest estimate in ``result.performance``.

    Returns a ``CheckResult``; ``.table`` has one row per protein: ``alone`` (the protein by
    itself), ``without`` (the rest of the panel) and ``drop`` (full panel minus without).
    """
    i = _inputs(result)
    model = model_key(model)
    panel = result.signature
    if not panel:
        raise ValueError("The signature is empty, so there is nothing to compare.")
    metric = _metric(i["task"])
    kw = dict(models=(model,), n_splits=i["n_splits"], n_repeats=i["n_repeats"], seed=i["seed"],
              classes=_classes(i))

    def run(cols):
        return _mean_metric(run_cv(i["P"], i["C"], i["y"], i["task"], fixed_panel=cols, **kw), model, metric)

    full = run(panel)
    rows = []
    for p in panel:
        alone = run([p])
        without = run([q for q in panel if q != p]) if len(panel) > 1 else (np.nan, np.nan, np.nan)
        rows.append({"protein": p, "alone": alone[0], "alone_ci_low": alone[1], "alone_ci_high": alone[2],
                     "without": without[0], "drop": full[0] - without[0]})
    table = pd.DataFrame(rows).sort_values("drop", ascending=False).reset_index(drop=True)
    name = "AUC" if metric == "auc" else "Spearman r"
    top = table.iloc[0]
    alone = table.sort_values("alone", ascending=False).iloc[0]
    lines = [f"full panel: {name} {full[0]:.2f}",
             f"largest drop: {top['protein']} ({name} {top['drop']:+.3f} when left out)",
             f"strongest alone: {alone['protein']} ({name} {alone['alone']:.2f})"]
    out = CheckResult(f"Protein contributions: {_n_prot(len(panel))}, {model_label(model, i['task'])}",
                      lines, table, metric)
    out.full = full[0]
    if verbose:
        print(out)
    return out


def learning_curve(result: SignatureResult, *, fractions=(0.3, 0.5, 0.7, 0.85, 1.0), model: str = "regression",
                   verbose: bool = True) -> CheckResult:
    """Performance as the training part grows. The test people stay the same; only the training
    people in each fold are reduced to a share of their number, and the proteins are chosen again
    on every smaller training set. A curve that still rises at the right edge says more people would help."""
    i = _inputs(result)
    model = model_key(model)
    metric = _metric(i["task"])
    rows = []
    for fr in fractions:
        folds = run_cv(i["P"], i["C"], i["y"], i["task"], models=(model,), n_splits=i["n_splits"],
                       n_repeats=i["n_repeats"], seed=i["seed"], train_fraction=fr, classes=_classes(i))
        m, lo, hi = _mean_metric(folds, model, metric)
        rows.append({"fraction": fr, "n_train": float(np.mean([f.n_train for f in folds])), "mean": m,
                     "ci_low": lo, "ci_high": hi,
                     "proteins_chosen": float(np.mean([len(f.selected) for f in folds]))})
    table = pd.DataFrame(rows)
    name = "AUC" if metric == "auc" else "Spearman r"
    first, last = table.iloc[0], table.iloc[-1]
    lines = [f"{name}: {first['mean']:.2f} with {first['n_train']:.0f} training people; "
             f"{last['mean']:.2f} with {last['n_train']:.0f}",
             f"last step: {table['mean'].iloc[-1] - table['mean'].iloc[-2]:+.3f}"]
    out = CheckResult("Learning curve", lines, table, metric)
    if verbose:
        print(out)
    return out


def selection_curve(result: SignatureResult, *, sizes=(1, 2, 3, 5, 8, 12, 20), model: str = "regression",
                    verbose: bool = True) -> CheckResult:
    """Performance when the signature holds the k strongest proteins of each training part.

    In every fold the proteins are ranked on the training people (a ridge model) and the top k are
    kept, so the left-out people play no part in the ranking. Shows how many proteins are enough."""
    i = _inputs(result)
    model = model_key(model)
    metric = _metric(i["task"])
    sizes = [s for s in sizes if s <= i["P"].shape[1]]
    rows = []
    for k in sizes:
        folds = run_cv(i["P"], i["C"], i["y"], i["task"], models=(model,), n_splits=i["n_splits"],
                       n_repeats=i["n_repeats"], seed=i["seed"], top_k=k, classes=_classes(i))
        m, lo, hi = _mean_metric(folds, model, metric)
        rows.append({"proteins": k, "mean": m, "ci_low": lo, "ci_high": hi})
    table = pd.DataFrame(rows)
    name = "AUC" if metric == "auc" else "Spearman r"
    best = table.loc[table["mean"].idxmax()]
    first_k, last_k = int(table["proteins"].iloc[0]), int(table["proteins"].iloc[-1])
    lines = [f"{name}: {table['mean'].iloc[0]:.2f} with {_n_prot(first_k)}; "
             f"{table['mean'].iloc[-1]:.2f} with {last_k}",
             f"highest: {name} {best['mean']:.2f} with {_n_prot(int(best['proteins']))}"]
    out = CheckResult("Selection curve", lines, table, metric)
    if verbose:
        print(out)
    return out
