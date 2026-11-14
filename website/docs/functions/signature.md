---
title: "Biomarker signatures from Olink data in Python: find_signature"
description: "signatomics.signature: find_signature (nested cross-validation for a binary, multiclass or continuous outcome), protein_contributions, learning_curve, selection_curve and their charts."
---

# Signature

Find a small set of proteins that predicts an outcome, with an honest estimate of how well it
predicts people the model has not seen. Then check it.
{ .lede }

```python
from signatomics.signature import (find_signature, protein_contributions, learning_curve,
                                   selection_curve, plot_signature_performance,
                                   plot_selection_frequency, plot_contributions, plot_curve)
```

Every function takes the table from [`to_wide`](analysis-table.md#to_wide) (or any table with
`proteins=[...]` given). The checks take the result of `find_signature` and rerun its folds.

## `find_signature`

Find a protein signature for an outcome and estimate how well it predicts.

```python
find_signature(table, outcome, *, positive=None, task="auto", proteins=None, covariates=None,
               models=("linear", "forest"), n_splits=5, n_repeats=5, min_frequency=0.5,
               seed=0, verbose=True) -> SignatureResult
```

| Argument | Type | Default | Meaning |
|---|---|---|---|
| `table` | DataFrame | — | the analysis table |
| `outcome` | str | — | the column to predict: two classes (binary), three or more (multiclass) or a number (continuous). People without a value are left out and counted |
| `positive` | | required for binary | the class to predict, e.g. `"case"`; predictions are the probability of this class |
| `task` | str | `"auto"` | `"binary"`, `"multiclass"` or `"continuous"`; `"auto"` decides from the column (a number with more than 5 values is continuous) |
| `proteins` | list | `None` | the assay columns; default from `to_wide` |
| `covariates` | list | `None` | columns that go into every model next to the proteins, e.g. `["age", "sex"]`, and are never chosen or dropped. The report then also gives the model on the covariates alone |
| `models` | tuple | `("linear", "forest")` | `"linear"`: logistic regression (ridge for a number); `"forest"`: random forest. Both are reported; none is picked as the best |
| `n_splits`, `n_repeats` | int | `5`, `5` | outer folds and repeats: 5 × 5 = 25 folds |
| `min_frequency` | float | `0.5` | a protein belongs to the signature when at least this share of the folds chose it |
| `seed` | int | `0` | makes the folds the same every time |
| `verbose` | bool | `True` | print the report |

**How it works.** The people are split into 5 parts; each part in turn is left out and the rest
is the training part. Inside each training part an L1 model chooses the proteins (its penalty set
by an inner cross-validation, with the one-standard-error rule), then the models are fitted on
those proteins and predict the people left out. This is repeated 5 times with new splits. The
left-out people play no part in choosing the proteins, so the performance is an honest estimate.
Its 95% interval allows for the folds sharing people (Nadeau and Bengio's correction).

**Output** (`SignatureResult`): it prints as a report, and holds

| Attribute | What it holds |
|---|---|
| `.performance` | one row per model and metric: `mean` over the folds with `ci_low`, `ci_high`. Binary: AUC; multiclass: AUC and balanced accuracy; continuous: Spearman r and RMSE. With covariates, also the covariates alone |
| `.proteins` | every protein chosen at least once: `frequency` (share of folds), `mean_coefficient_size`, `in_signature` |
| `.signature` | the proteins chosen in at least `min_frequency` of the folds |
| `.folds` | the metric of every fold |
| `.predictions` | every person's prediction from the folds where the person was left out |
| `.predict(table, model="linear")` | predictions for new people from the final model, fitted on everyone with the signature proteins |

Only `.performance` is an honest estimate: the final model behind `.predict` has seen everyone.

**Errors.** `KeyError` when the outcome or a covariate is not in the table; `ValueError` for a
binary outcome without `positive`, a class smaller than `n_splits`, or too few people.

## `protein_contributions`

What each signature protein adds: the panel without it, and the protein on its own.

```python
protein_contributions(result, *, model="linear", verbose=True) -> CheckResult
```

`.table`: one row per protein with `alone` (with its interval), `without` (the rest of the panel)
and `drop` (full panel minus without). The folds are rerun with the signature fixed, so the
numbers compare the proteins with each other; they sit above the honest estimate in
`result.performance`.

## `learning_curve`

Performance as the training part grows.

```python
learning_curve(result, *, fractions=(0.3, 0.5, 0.7, 0.85, 1.0), model="linear",
               verbose=True) -> CheckResult
```

The left-out people stay the same; the training part shrinks to each fraction and the proteins
are chosen again. `.table`: `fraction`, `n_train`, `mean`, `ci_low`, `ci_high`,
`proteins_chosen`. A curve still rising at the right says more people would help.

## `selection_curve`

Performance with the 1, 2, 3, ... strongest proteins.

```python
selection_curve(result, *, sizes=(1, 2, 3, 5, 8, 12, 20), model="linear",
                verbose=True) -> CheckResult
```

In every fold the proteins are ranked on the training part and the top k kept. `.table`:
`proteins`, `mean`, `ci_low`, `ci_high`. Shows how many proteins are enough; the intervals
overlap, so read the shape of the curve.

## Charts

| Chart | Draws |
|---|---|
| `plot_signature_performance(result)` | binary: the ROC curve of each model, one line per repeat; continuous: predicted against observed |
| `plot_selection_frequency(result, *, top=20)` | how often each protein was chosen, signature proteins coloured |
| `plot_contributions(check)` | the result of `protein_contributions`: each protein alone and the panel without it |
| `plot_curve(check)` | the result of `learning_curve` or `selection_curve`, with its 95% interval |

<figure class="sg-chart" markdown>
<div class="sg-plot" data-src="../../assets/figures/example_signature_performance.json"><img src="../../assets/figures/example_signature_performance.svg" alt="ROC curves of the signature on the synthetic data" loading="lazy"></div>
<figcaption markdown>`plot_signature_performance(result)` on the synthetic data. See [Find a signature](../learn/find-signature.md) for how to read it.</figcaption>
</figure>

**Example**

```python
result = find_signature(table, "group", positive="case")
check = protein_contributions(result)
plot_contributions(check)
```
