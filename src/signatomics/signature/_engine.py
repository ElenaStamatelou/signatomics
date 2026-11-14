"""The cross-validation engine behind ``find_signature`` and its checks.

Everything that learns from the outcome happens inside the training part of each outer fold:
imputing, scaling, choosing the proteins and fitting the model. The test part of the fold is only
predicted. The chosen proteins of every fold are kept, so their stability can be reported.
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass, field

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor
from sklearn.impute import SimpleImputer
from sklearn.linear_model import Lasso, LogisticRegression, Ridge
from sklearn.metrics import balanced_accuracy_score, roc_auc_score
from sklearn.model_selection import KFold, RepeatedKFold, RepeatedStratifiedKFold, StratifiedKFold, cross_val_score
from sklearn.preprocessing import StandardScaler

MODELS = ("linear", "forest")
TASKS = ("binary", "multiclass", "continuous")


@dataclass
class Fold:
    repeat: int
    fold: int
    n_train: int
    n_test: int
    selected: list
    coef: dict
    test_index: np.ndarray
    predictions: dict            # model -> array (binary: P(positive); multiclass: n x classes; continuous: value)
    metrics: dict = field(default_factory=dict)   # model -> {metric: value}


def _quiet():
    warnings.filterwarnings("ignore", category=FutureWarning, module="sklearn")
    warnings.filterwarnings("ignore", category=UserWarning, module="sklearn")
    warnings.filterwarnings("ignore", category=Warning, module="sklearn.linear_model")


def _l1_logreg(C: float, multiclass: bool):
    solver = "saga" if multiclass else "liblinear"
    try:
        return LogisticRegression(C=C, l1_ratio=1.0, solver=solver, max_iter=3000)
    except TypeError:                       # an older scikit-learn
        return LogisticRegression(C=C, penalty="l1", solver=solver, max_iter=3000)


def _one_se(scores: np.ndarray) -> int:
    """Index of the sparsest setting (the first, settings run from sparse to dense) whose mean score
    is within one standard error of the best: the usual 'one-standard-error rule'."""
    mean = scores.mean(axis=1)
    se = scores.std(axis=1, ddof=1) / np.sqrt(scores.shape[1])
    best = int(np.nanargmax(mean))
    return int(np.flatnonzero(mean >= mean[best] - se[best])[0])


def _select(Z: np.ndarray, y: np.ndarray, task: str, top_k, inner: int, seed: int):
    """Choose the proteins on training data: a sparse (L1) model, or the top_k of a ridge model.

    Returns the chosen column numbers and a score per column (the size of its coefficient)."""
    p = Z.shape[1]
    if top_k is not None:
        if task == "continuous":
            coef = Ridge(alpha=10.0).fit(Z, y).coef_
        else:
            coef = LogisticRegression(C=0.05, max_iter=5000).fit(Z, y).coef_
        coef = np.atleast_2d(coef)
        score = np.abs(coef).max(axis=0)
        order = np.argsort(-score)[:min(top_k, p)]
        return sorted(order.tolist()), score
    folds = list((KFold if task == "continuous" else StratifiedKFold)(inner, shuffle=True, random_state=seed)
                 .split(Z, y))
    if task == "continuous":
        grid = np.logspace(0.2, -2.5, 12) * np.abs(Z.T @ (y - y.mean())).max() / len(y)    # strong to weak penalty
        scores = np.array([[-np.mean((y[te] - Lasso(alpha=a, max_iter=20000).fit(Z[tr], y[tr]).predict(Z[te])) ** 2)
                            for tr, te in folds] for a in grid])
        coef = np.atleast_2d(Lasso(alpha=grid[_one_se(scores)], max_iter=20000).fit(Z, y).coef_)
    else:
        grid = np.logspace(-2.2, 0.3, 10)                                                  # strong to weak penalty
        multi = task == "multiclass"
        scoring = "roc_auc_ovr" if multi else "roc_auc"
        scores = np.array([cross_val_score(_l1_logreg(c, multi), Z, y, cv=folds, scoring=scoring) for c in grid])
        coef = np.atleast_2d(_l1_logreg(grid[_one_se(scores)], multi).fit(Z, y).coef_)
    score = np.abs(coef).max(axis=0)
    return np.flatnonzero(score > 0).tolist(), (coef[0] if coef.shape[0] == 1 else score)


def _predictor(model: str, task: str, seed: int):
    if model == "linear":
        return Ridge(alpha=1.0) if task == "continuous" else LogisticRegression(C=1.0, max_iter=5000)
    if task == "continuous":
        return RandomForestRegressor(300, min_samples_leaf=3, random_state=seed, n_jobs=1)
    return RandomForestClassifier(300, min_samples_leaf=3, class_weight="balanced", random_state=seed, n_jobs=1)


def _predict(est, M, task):
    if task == "continuous":
        return est.predict(M)
    proba = est.predict_proba(M)
    return proba[:, 1] if task == "binary" else proba


def _score(task, y_true, pred, classes):
    if task == "binary":
        return {"auc": roc_auc_score(y_true, pred)}
    if task == "multiclass":
        return {"auc": roc_auc_score(y_true, pred, multi_class="ovr", average="macro",
                                     labels=list(range(len(classes)))),
                "balanced_accuracy": balanced_accuracy_score(y_true, pred.argmax(axis=1))}
    rho = stats.spearmanr(y_true, pred)[0] if np.ptp(pred) > 0 else np.nan
    return {"rmse": float(np.sqrt(np.mean((y_true - pred) ** 2))), "spearman": float(rho)}


def run_cv(P: pd.DataFrame, C: pd.DataFrame | None, y: np.ndarray, task: str, *, models=MODELS,
           n_splits=5, n_repeats=5, inner=5, seed=0, top_k=None, train_fraction=1.0,
           fixed_panel: list | None = None, classes=None) -> list[Fold]:
    """Run the repeated outer cross-validation and return one ``Fold`` per fold.

    P: proteins (rows = samples), C: covariates already numeric (or None), y: 0/1 (binary), class
    codes (multiclass) or values (continuous). ``fixed_panel`` skips the selection and uses these
    protein columns in every fold (only for comparing proteins: its AUC is optimistic).
    """
    _quiet()
    rng = np.random.default_rng(seed)
    Pv = P.to_numpy(dtype=float)
    Cv = None if C is None or C.shape[1] == 0 else C.to_numpy(dtype=float)
    cv = (RepeatedKFold if task == "continuous" else RepeatedStratifiedKFold)(
        n_splits=n_splits, n_repeats=n_repeats, random_state=seed)
    names = list(P.columns)
    folds = []
    for i, (tr, te) in enumerate(cv.split(Pv, y)):
        if train_fraction < 1.0:                 # a smaller training set, for the learning curve
            if task == "continuous":
                tr = rng.choice(tr, size=max(int(round(len(tr) * train_fraction)), 10), replace=False)
            else:                                # every class stays in the smaller training set
                tr = np.concatenate([
                    rng.choice(tr[y[tr] == c], size=max(int(round((y[tr] == c).sum() * train_fraction)), 3),
                               replace=False) for c in np.unique(y[tr])])
        imp = SimpleImputer(strategy="median").fit(Pv[tr])
        sc = StandardScaler().fit(imp.transform(Pv[tr]))
        Ztr, Zte = sc.transform(imp.transform(Pv[tr])), sc.transform(imp.transform(Pv[te]))
        if fixed_panel is not None:
            sel = [names.index(c) for c in fixed_panel]
            score = np.zeros(len(names))
        else:
            sel, score = _select(Ztr, y[tr], task, top_k, inner, seed + i)
        blocks_tr, blocks_te = [], []
        if sel:
            blocks_tr.append(Ztr[:, sel])
            blocks_te.append(Zte[:, sel])
        if Cv is not None:
            cimp = SimpleImputer(strategy="median").fit(Cv[tr])
            csc = StandardScaler().fit(cimp.transform(Cv[tr]))
            blocks_tr.append(csc.transform(cimp.transform(Cv[tr])))
            blocks_te.append(csc.transform(cimp.transform(Cv[te])))
        preds, mets = {}, {}
        for model in models:
            if not blocks_tr:                              # nothing chosen: predict what the training data say
                if task == "multiclass":
                    pr = np.tile(np.bincount(y[tr], minlength=len(classes)) / len(tr), (len(te), 1))
                else:
                    pr = np.full(len(te), y[tr].mean())
            else:
                est = _predictor(model, task, seed + i).fit(np.hstack(blocks_tr), y[tr])
                pr = _predict(est, np.hstack(blocks_te), task)
                if task == "multiclass" and pr.shape[1] != len(classes):     # a class missing in training
                    full = np.zeros((len(te), len(classes)))
                    full[:, est.classes_] = pr
                    pr = full
            preds[model] = pr
            mets[model] = _score(task, y[te], pr, classes)
        coef = {names[j]: float(score[j]) for j in sel}
        folds.append(Fold(repeat=i // n_splits, fold=i % n_splits, n_train=len(tr), n_test=len(te),
                          selected=[names[j] for j in sel], coef=coef, test_index=te,
                          predictions=preds, metrics=mets))
    return folds


# the range each metric can take: its interval is clipped to it (an AUC cannot pass 1)
METRIC_RANGE = {"auc": (0.0, 1.0), "balanced_accuracy": (0.0, 1.0), "spearman": (-1.0, 1.0),
                "rmse": (0.0, None)}


def corrected_ci(values, n_train, n_test, level=0.95, metric=None):
    """Interval for a metric averaged over repeated k-fold CV (Nadeau and Bengio 2003): the folds
    share people, so the spread of the fold values is inflated by (1/J + n_test/n_train).
    With ``metric`` the interval is clipped to the range that metric can take."""
    v = np.asarray(values, dtype=float)
    v = v[~np.isnan(v)]
    j = len(v)
    if j < 2:
        return np.nan, np.nan
    se = np.sqrt((1 / j + n_test / n_train) * np.var(v, ddof=1))
    q = stats.t.ppf(0.5 + level / 2, j - 1)
    lo, hi = float(v.mean() - q * se), float(v.mean() + q * se)
    low, high = METRIC_RANGE.get(metric, (None, None))
    if low is not None:
        lo = max(lo, low)
    if high is not None:
        hi = min(hi, high)
    return lo, hi
