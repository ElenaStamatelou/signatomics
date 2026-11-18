"""find_signature: a small set of proteins that predicts an outcome, with an honest estimate.

The estimate comes from nested cross-validation. In every fold the proteins are chosen on the
training people only (the lasso, with its own inner cross-validation picking the penalty by the
one-standard-error rule: the sparsest model within one standard error of the best), a
model is fitted on them, and the people left out are predicted. Nothing about the left-out
people is used before that. The proteins that are chosen in many folds form the signature.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pandas as pd
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import StandardScaler

from ..stats._common import _n
from ..stats.adjust import _design_covariates
from ._engine import MODELS, TASKS, _predict, _predictor, _quiet, _select, corrected_ci, model_key, model_label, run_cv
from .._report import format_report, table_lines


@dataclass
class SignatureResult:
    """Result of ``find_signature``: summary lines plus the tables behind them.

    ``performance``  one row per model (and for the covariates alone, when you gave covariates):
                     the metric averaged over all outer folds with a 95% interval that allows for
                     the folds sharing people
    ``proteins``     every protein chosen at least once: the share of folds that chose it
                     (``frequency``), its mean coefficient size, and ``in_signature``
    ``signature``    the proteins chosen in at least ``min_frequency`` of the folds
    ``folds``        the metric of every outer fold
    ``predictions``  the prediction for every person from the folds where the person was left out
    ``predict(table)`` predictions for new people from the final model (fitted on everyone, on the
                     signature proteins); only the performance table is an honest estimate
    """

    title: str = "Biomarker signature"
    lines: list = field(default_factory=list)
    task: str = ""
    outcome: str = ""
    classes: list = field(default_factory=list)
    performance: pd.DataFrame | None = None
    proteins: pd.DataFrame | None = None
    signature: list = field(default_factory=list)
    folds: pd.DataFrame | None = None
    predictions: pd.DataFrame | None = None
    _final: object = None
    _inputs: dict = field(default_factory=dict, repr=False)   # what the checks need to rerun the folds

    def __str__(self) -> str:
        return format_report(self.title, self.lines)

    def __repr__(self) -> str:       # notebooks show the same text as print()
        return str(self)

    def predict(self, table: pd.DataFrame, model: str = "regression") -> pd.Series | pd.DataFrame:
        """Predictions for new people from the final model (proteins and covariates as columns)."""
        if self._final is None:
            raise ValueError("No signature was found, so there is no model to predict with.")
        return self._final.predict(table, model_key(model))


class _Final:
    """The final model: fitted on all people, on the signature proteins and the covariates."""

    def __init__(self, P, C, y, task, signature, models, seed, classes):
        self.task, self.signature, self.classes = task, signature, classes
        self.cov_columns = None if C is None else list(C.columns)
        Pv = P[signature].to_numpy(dtype=float)
        self.p_imp = SimpleImputer(strategy="median").fit(Pv)
        self.p_sc = StandardScaler().fit(self.p_imp.transform(Pv))
        blocks = [self.p_sc.transform(self.p_imp.transform(Pv))]
        if C is not None and C.shape[1]:
            self.c_imp = SimpleImputer(strategy="median").fit(C.to_numpy(dtype=float))
            self.c_sc = StandardScaler().fit(self.c_imp.transform(C.to_numpy(dtype=float)))
            blocks.append(self.c_sc.transform(self.c_imp.transform(C.to_numpy(dtype=float))))
        else:
            self.c_imp = None
        self.est = {m: _predictor(m, task, seed).fit(np.hstack(blocks), y) for m in models}
        self.covariates = None

    def predict(self, table, model):
        if model not in self.est:
            raise ValueError(f"Model '{model}' was not fitted. Fitted: {', '.join(self.est)}")
        missing = [p for p in self.signature if p not in table.columns]
        if missing:
            raise KeyError(f"Signature proteins not in the table: {', '.join(missing)}")
        blocks = [self.p_sc.transform(self.p_imp.transform(table[self.signature].to_numpy(dtype=float)))]
        if self.c_imp is not None:
            C = _design_covariates(table, self.covariates).reindex(columns=self.cov_columns, fill_value=0.0)
            blocks.append(self.c_sc.transform(self.c_imp.transform(C.to_numpy(dtype=float))))
        pr = _predict(self.est[model], np.hstack(blocks), self.task)
        if self.task == "multiclass":
            return pd.DataFrame(pr, index=table.index, columns=[f"p_{c}" for c in self.classes])
        return pd.Series(pr, index=table.index, name="predicted")


def prepare_xy(table: pd.DataFrame, outcome: str, *, positive, task: str, proteins, covariates):
    """Shared by find_signature and the checks: the protein table, covariates and outcome codes."""
    if outcome not in table.columns:
        raise KeyError(f"Column '{outcome}' not found in the table.")
    if task not in ("auto",) + TASKS:
        raise ValueError(f"task must be 'auto' or one of: {', '.join(TASKS)}")
    covariates = [covariates] if isinstance(covariates, str) else list(covariates or [])
    absent = [c for c in covariates if c not in table.columns]
    if absent:
        raise KeyError(f"Covariate columns not in the table: {', '.join(absent)}")
    if proteins is None:
        proteins = table.attrs.get("proteins")
        if proteins is None:
            raise ValueError("Which columns are proteins? Make the table with to_wide, or pass "
                             "proteins=[...] (a list of column names).")
    proteins = [p for p in proteins if p != outcome and p not in covariates]
    raw = table[outcome]
    if task == "auto":
        if pd.api.types.is_numeric_dtype(raw) and not pd.api.types.is_bool_dtype(raw) and raw.nunique() > 5:
            task = "continuous"
        else:
            task = "binary" if raw.dropna().nunique() == 2 else "multiclass"
    keep = raw.notna()
    if covariates:
        keep &= table[covariates].notna().all(axis=1)
    info = {"no_outcome": int(raw.isna().sum()),
            "no_covariate": int((raw.notna() & ~keep).sum())}
    data = table.loc[keep]
    y_raw = data[outcome]
    classes = []
    if task == "continuous":
        y = pd.to_numeric(y_raw).to_numpy(dtype=float)
    else:
        classes = sorted(y_raw.unique(), key=str) if task == "multiclass" else None
        if task == "binary":
            levels = sorted(y_raw.unique(), key=str)
            if len(levels) != 2:
                raise ValueError(f"'{outcome}' has {len(levels)} classes; a binary outcome has two.")
            if positive is None:
                opts = [f"positive={g!r}" for g in levels]
                raise ValueError(f"Which class is the one to predict? Pass {opts[0]} or {opts[1]}.")
            if positive not in levels:
                raise ValueError(f"positive={positive!r} not found. Classes: {', '.join(map(str, levels))}")
            classes = [g for g in levels if g != positive] + [positive]
            y = (y_raw == positive).to_numpy(dtype=int)
        else:
            y = pd.Categorical(y_raw, categories=classes).codes.astype(int)
    P = data[proteins].apply(pd.to_numeric, errors="coerce")
    C = _design_covariates(data, covariates) if covariates else None
    return data, P, C, y, task, classes, covariates, info


def _summary(folds, models, task, label_extra=""):
    metrics = sorted({m for f in folds for m in f.metrics[models[0]]})
    rows = []
    for model in models:
        for metric in metrics:
            vals = [f.metrics[model][metric] for f in folds]
            lo, hi = corrected_ci(vals, np.mean([f.n_train for f in folds]), np.mean([f.n_test for f in folds]),
                              metric=metric)
            rows.append({"model": model + label_extra, "metric": metric, "mean": float(np.nanmean(vals)),
                         "ci_low": lo, "ci_high": hi})
    return rows


def find_signature(table: pd.DataFrame, outcome: str, *, positive=None, task: str = "auto",
                   proteins=None, covariates=None, models=MODELS, n_splits: int = 5,
                   n_repeats: int = 5, min_frequency: float = 0.5, seed: int = 0,
                   verbose: bool = True) -> SignatureResult:
    """Find a protein signature for an outcome and estimate how well it predicts.

    Parameters
    ----------
    table : the analysis table from ``to_wide`` (one row per sample, one column per protein, your
        sample metadata included).
    outcome : the column to predict. Two classes (binary), three or more classes (multiclass) or
        a number (continuous); ``task="auto"`` decides from the column.
    positive : binary only, required: the class to predict, e.g. ``"T2DM"``. Predictions are the
        probability of this class.
    task : ``"auto"`` (default), ``"binary"``, ``"multiclass"`` or ``"continuous"``.
    proteins : the protein columns; default from ``to_wide``.
    covariates : columns that go into every model next to the proteins (e.g. ``["age", "sex"]``)
        and are never chosen or dropped. The report then also gives the model on the covariates
        alone: the proteins add what the gap between the two shows.
    models : ``"regression"`` (logistic regression; ridge regression for a number) and
        ``"random_forest"`` (random forest). Both are reported; no model is picked as the best.
        The names before version 0.2.0, ``"linear"`` and ``"forest"``, still work.
    n_splits, n_repeats : outer folds and repeats (default 5 x 5 = 25 folds).
    min_frequency : a protein belongs to the signature when at least this share of the folds chose
        it (default 0.5).
    seed : makes the folds reproducible.
    verbose : print the report.

    Returns a ``SignatureResult``.
    """
    models = tuple(model_key(m) for m in ([models] if isinstance(models, str) else models))
    data, P, C, y, task, classes, covariates, info = prepare_xy(
        table, outcome, positive=positive, task=task, proteins=proteins, covariates=covariates)
    n = len(y)
    if task != "continuous":
        counts = np.bincount(y, minlength=len(classes))
        if counts.min() < n_splits:
            raise ValueError(f"The smallest class has {counts.min()} people; with n_splits={n_splits} "
                             f"every fold needs at least one. Lower n_splits or merge classes.")
    elif n < 2 * n_splits:
        raise ValueError(f"Only {n} people for {n_splits} folds.")

    folds = run_cv(P, C, y, task, models=models, n_splits=n_splits, n_repeats=n_repeats, seed=seed,
                   classes=classes if task == "multiclass" else None)
    perf = _summary(folds, models, task)
    if covariates:
        base = run_cv(P, C, y, task, models=models, n_splits=n_splits, n_repeats=n_repeats, seed=seed,
                      fixed_panel=[], classes=classes if task == "multiclass" else None)
        perf += _summary(base, models, task, label_extra=" (covariates only)")
    performance = pd.DataFrame(perf)

    # how often each protein was chosen
    n_folds = len(folds)
    freq = pd.Series([p for f in folds for p in f.selected]).value_counts() / n_folds
    size = pd.DataFrame([f.coef for f in folds]).abs().mean() if any(f.coef for f in folds) else pd.Series(dtype=float)
    prot = pd.DataFrame({"protein": freq.index, "frequency": freq.to_numpy()})
    prot["mean_coefficient_size"] = prot["protein"].map(size).to_numpy()
    prot["in_signature"] = prot["frequency"] >= min_frequency
    prot = prot.sort_values(["frequency", "mean_coefficient_size"], ascending=False).reset_index(drop=True)
    signature = list(prot.loc[prot["in_signature"], "protein"])

    fold_rows = [{"repeat": f.repeat, "fold": f.fold, "n_train": f.n_train, "n_test": f.n_test,
                  "n_selected": len(f.selected), "model": m, **f.metrics[m]} for f in folds for m in models]
    pred_rows = []
    for f in folds:
        for m in models:
            pr = f.predictions[m]
            part = pd.DataFrame({"sample": data.index[f.test_index], "repeat": f.repeat, "model": m,
                                 "observed": data[outcome].to_numpy()[f.test_index]})
            if task == "multiclass":
                for j, c in enumerate(classes):
                    part[f"p_{c}"] = pr[:, j]
            else:
                part["predicted"] = pr
            pred_rows.append(part)
    predictions = pd.concat(pred_rows, ignore_index=True)

    res = SignatureResult(task=task, outcome=outcome, classes=classes, performance=performance,
                          proteins=prot, signature=signature, folds=pd.DataFrame(fold_rows),
                          predictions=predictions)
    res._inputs = dict(P=P, C=C, y=y, task=task, classes=classes, models=models, n_splits=n_splits,
                       n_repeats=n_repeats, seed=seed)
    if signature:
        res._final = _Final(P, C, y, task, signature, models, seed, classes)
        res._final.covariates = covariates
    lines = _report_lines(res, n, P.shape[1], n_folds, n_splits, n_repeats, min_frequency, info,
                          covariates, y if task != "continuous" else None, classes)
    res.lines = lines
    if verbose:
        print(res)
    return res


def _report_lines(res, n, n_prot, n_folds, n_splits, n_repeats, min_frequency, info, covariates, y, classes):
    task = res.task
    lines = []
    if task == "binary":
        what = f"{res.outcome} = {classes[-1]}"
    elif task == "multiclass":
        what = f"{res.outcome} ({len(classes)} classes)"
    else:
        what = f"{res.outcome} (a number)"
    res.title = f"Biomarker signature: {what}"
    if task == "binary":
        k = int(y.sum())
        lines.append(f"people: {n} ({k} {classes[-1]}, {n - k} {classes[0]})")
    elif task == "multiclass":
        lines.append(f"people: {n} (" + ", ".join(f"{c} {int((y == i).sum())}" for i, c in enumerate(classes))
                     + ")")
    else:
        lines.append(f"people: {n}")
    lines.append(f"proteins: {n_prot}")
    if info["no_outcome"]:
        lines.append(f"left out: {_n(info['no_outcome'], 'sample')} without a value in '{res.outcome}'")
    if info["no_covariate"]:
        lines.append(f"left out: {_n(info['no_covariate'], 'sample')} with a missing covariate")
    lines.append(f"method: nested cross-validation ({n_splits} folds x {n_repeats} repeats)")
    lines.append("selection: lasso, inside each training fold")
    if covariates:
        lines.append(f"covariates: {', '.join(covariates)} (in every model, never dropped)")
    perf = res.performance
    if task == "binary":
        head, extra = "AUC (95% interval)", None
    elif task == "multiclass":
        head, extra = "AUC (95% interval), balanced accuracy", "balanced_accuracy"
    else:
        head, extra = "Spearman r (95% interval), RMSE", "rmse"
    primary = "spearman" if task == "continuous" else "auc"
    rows = []
    for _, r in perf[perf["metric"] == primary].iterrows():
        row = [model_label(r["model"], task), f"{r['mean']:.2f}", f"{r['ci_low']:.2f} to {r['ci_high']:.2f}"]
        if extra:
            e = perf[(perf["metric"] == extra) & (perf["model"] == r["model"])]["mean"]
            row.append(f"{e.iloc[0]:.2f}" if len(e) else "")
        rows.append(row)
    lines.append(f"performance: {head}")
    lines += table_lines(rows)
    empty = int((res.folds[res.folds["model"] == res.folds["model"].iloc[0]]["n_selected"] == 0).sum())
    if empty:
        lines.append(f"{_n(empty, 'fold')} chose no protein; those folds predict what the training people "
                     f"show, and count in the averages")
    sig = res.signature
    if sig:
        top = res.proteins[res.proteins["in_signature"]].head(8)
        lines.append(f"signature: {_n(len(sig), 'protein')}, chosen in at least {min_frequency:.0%} of the folds")
        lines += table_lines([[str(r["protein"]), f"{r['frequency']:.0%}"] for _, r in top.iterrows()])
        if len(sig) > 8:
            lines.append(f"  ... {len(sig)} in total, in .signature")
    else:
        lines.append(f"signature: none (no protein chosen in at least {min_frequency:.0%} of the folds)")
    return lines
