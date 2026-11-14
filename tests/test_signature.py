"""Tests for signatomics.signature (find_signature and its checks) on made-up data with known signal."""

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import pytest

from signatomics.signature import (find_signature, learning_curve, plot_contributions, plot_curve,
                                   plot_selection_frequency, plot_signature_performance,
                                   protein_contributions, selection_curve)

FAST = dict(n_splits=4, n_repeats=1, models=("linear",), verbose=False)


def _table(n=120, p=40, signal=(0, 1, 2), effect=1.0, seed=0, kind="binary"):
    rng = np.random.default_rng(seed)
    X = rng.normal(size=(n, p))
    if kind == "binary":
        y = np.array(["control", "case"])[(rng.random(n) < 0.5).astype(int)]
        for j in signal:
            X[y == "case", j] += effect
    elif kind == "multiclass":
        y = rng.choice(["a", "b", "c"], size=n)
        for j, c in zip(signal, ["a", "b", "c"]):
            X[y == c, j] += 1.5
    else:
        y = X[:, list(signal)].sum(axis=1) + rng.normal(0, 0.5, n)
    t = pd.DataFrame(X, columns=[f"P{j}" for j in range(p)])
    t["outcome"] = y
    t["age"] = rng.normal(60, 8, n)
    t.attrs["proteins"] = [f"P{j}" for j in range(p)]
    return t


def test_binary_finds_the_built_in_signal():
    r = find_signature(_table(), "outcome", positive="case", **FAST)
    auc = r.performance.query("metric == 'auc'")["mean"].iloc[0]
    assert auc > 0.8
    assert {"P0", "P1", "P2"} <= set(r.signature)
    assert len(r.signature) <= 12
    assert r.proteins["frequency"].between(0, 1).all()


def test_noise_gives_no_skill_because_selection_stays_inside_the_folds():
    # 150 proteins of pure noise and 50 people: choosing proteins on everyone first would give a high AUC
    t = _table(n=50, p=150, signal=(), seed=4)
    r = find_signature(t, "outcome", positive="case", n_splits=5, n_repeats=2, models=("linear",), verbose=False)
    row = r.performance.query("metric == 'auc'").iloc[0]
    assert row["ci_low"] < 0.5 < row["ci_high"] and row["mean"] < 0.65


def test_the_positive_class_is_required_and_checked():
    t = _table()
    with pytest.raises(ValueError, match="positive='case' or positive='control'"):
        find_signature(t, "outcome", **FAST)
    with pytest.raises(ValueError, match="not found"):
        find_signature(t, "outcome", positive="sick", **FAST)


def test_multiclass_and_continuous():
    m = find_signature(_table(kind="multiclass", n=150), "outcome", **FAST)
    assert m.task == "multiclass" and m.performance.query("metric == 'auc'")["mean"].iloc[0] > 0.7
    c = find_signature(_table(kind="continuous"), "outcome", **FAST)
    assert c.task == "continuous" and c.performance.query("metric == 'spearman'")["mean"].iloc[0] > 0.7
    assert {"P0", "P1", "P2"} <= set(c.signature)


def test_covariates_get_a_baseline_and_are_never_chosen():
    t = _table()
    t["age"] = t["age"] + (t["outcome"] == "case") * 10                       # age predicts the outcome too
    r = find_signature(t, "outcome", positive="case", covariates=["age"], **FAST)
    models = set(r.performance["model"])
    assert models == {"linear", "linear (covariates only)"}
    assert "age" not in r.proteins["protein"].tolist()
    base = r.performance.query("model == 'linear (covariates only)' and metric == 'auc'")["mean"].iloc[0]
    full = r.performance.query("model == 'linear' and metric == 'auc'")["mean"].iloc[0]
    assert full > base > 0.6


def test_missing_outcome_and_covariate_rows_are_left_out_and_reported():
    t = _table()
    t.loc[t.index[:3], "outcome"] = np.nan
    t.loc[t.index[3:5], "age"] = np.nan
    r = find_signature(t, "outcome", positive="case", covariates=["age"], **FAST)
    text = str(r)
    assert "3 samples without a value" in text and "2 samples with a missing covariate" in text


def test_predict_gives_probabilities_for_new_people():
    t = _table()
    r = find_signature(t, "outcome", positive="case", **FAST)
    p = r.predict(t.head(10))
    assert len(p) == 10 and p.between(0, 1).all()
    p_forest = find_signature(t, "outcome", positive="case", n_splits=4, n_repeats=1, verbose=False).predict(t.head(5), "forest")
    assert len(p_forest) == 5
    with pytest.raises(KeyError):
        r.predict(t.drop(columns=r.signature[0]))


def test_same_seed_same_result():
    t = _table()
    a = find_signature(t, "outcome", positive="case", **FAST).performance
    b = find_signature(t, "outcome", positive="case", **FAST).performance
    pd.testing.assert_frame_equal(a, b)


def test_checks_and_charts():
    t = _table()
    r = find_signature(t, "outcome", positive="case", **FAST)
    pc = protein_contributions(r, verbose=False)
    assert set(pc.table["protein"]) == set(r.signature) and pc.table["alone"].between(0, 1).all()
    lc = learning_curve(r, fractions=(0.4, 1.0), verbose=False)
    assert lc.table["n_train"].iloc[0] < lc.table["n_train"].iloc[1]
    sc = selection_curve(r, sizes=(1, 3, 10), verbose=False)
    assert sc.table["mean"].iloc[1] > sc.table["mean"].iloc[0] - 0.05
    for fig in (plot_signature_performance(r), plot_selection_frequency(r), plot_curve(lc), plot_curve(sc),
                plot_contributions(pc)):
        assert isinstance(fig, go.Figure)
    with pytest.raises(ValueError):
        plot_curve(pc)


def test_too_small_a_class_is_refused():
    t = _table(n=40)
    t.loc[t.index[:38], "outcome"] = "control"
    with pytest.raises(ValueError, match="smallest class"):
        find_signature(t, "outcome", positive="case", n_splits=5, n_repeats=1, verbose=False)


def test_interval_stays_inside_the_range_of_the_metric():
    from signatomics.signature._engine import corrected_ci
    lo, hi = corrected_ci([0.99, 1.0, 0.95, 1.0, 0.97], 40, 10, metric="auc")
    assert hi == 1.0 and lo < 1.0
    lo, hi = corrected_ci([-0.99, -1.0, -0.9], 40, 10, metric="spearman")
    assert lo == -1.0
