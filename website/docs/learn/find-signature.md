---
title: "Find a biomarker signature in Olink data with Python"
description: "find_signature chooses a small set of proteins that predicts an outcome, with nested cross-validation for an honest estimate of how well it predicts new people."
hide:
  - toc
---

<p class="eyebrow">Machine learning · step 1</p>

# Find a signature

<ol class="sg-route" aria-label="Machine learning route"><li class="done here"><span class="dot"></span>Find a signature</li><li class="done"><span class="dot"></span><a href="../check-signature/">Check the signature</a></li></ol>

**Use this when** you want a small set of proteins that predicts an outcome (cases or controls,
a disease stage, a lab value) and an honest answer to the question: how well would it predict
people it has never seen?
{ .lede }

```python
from signatomics.signature import find_signature, plot_signature_performance, plot_selection_frequency

result = find_signature(table, "group", positive="case")
```

`table` is the [analysis table](analysis-table.md). `positive` names the class to predict: the
predictions are the chance that a person is a case.

<div class="sg-example" markdown>
<p class="tag">On the synthetic example files, P07 and P33 removed in quality control</p>

```text
Biomarker signature
  predicting 'group' = case (binary)
  people: 58 (29 case, 29 control); 92 proteins
  method: nested cross-validation, 5 folds x 5 repeats = 25 folds; in every training part an L1 model (penalty chosen by inner cross-validation, the one-standard-error rule) chooses the proteins
  linear: AUC 0.98 (95% interval 0.96 to 1.00)
  forest: AUC 0.97 (95% interval 0.94 to 1.00)
  signature: 4 proteins chosen in at least 50% of the folds: IL-17C (100%), OSM (100%), IL6 (96%), TWEAK (52%)
  the performance above comes from people left out of every step that chose proteins; the final model (.predict) is fitted on everyone
```

</div>

### How it stays honest

The people are split into 5 parts. Each part in turn is left out; on the other 4 the proteins are
chosen and the model is fitted, and then the left-out people are predicted. The whole split is
repeated 5 times, so 25 folds in all. The people that measure the performance never help to
choose the proteins.

Choosing the proteins on everyone first and then cross-validating is the usual mistake: the
proteins have then already seen the test people, and the performance comes out too good. On
data with no signal at all, it can look excellent.

| Report line | What it means |
|---|---|
| `linear`, `forest` | two kinds of model, both reported: a logistic regression and a random forest. When they agree, the result does not hang on the choice of model |
| `AUC` | the chance that the model ranks a random case above a random control: 0.5 is a coin toss, 1 is perfect |
| `95% interval` | the range the AUC is compatible with, widened because the folds share people |
| `signature` | the proteins chosen in at least half of the folds, with how often |

<figure class="sg-chart" markdown>
<div class="sg-plot" data-src="../../assets/figures/example_signature_performance.json"><img src="../../assets/figures/example_signature_performance.svg" alt="ROC curves of both models, close to the top-left corner" loading="lazy"></div>
<figcaption markdown>`plot_signature_performance(result)`: the ROC curve, one thin line per repeat. The closer to the top-left corner, the better; the dashed diagonal is a coin toss. Every point comes from people left out of the fold.</figcaption>
</figure>

<figure class="sg-chart" markdown>
<div class="sg-plot" data-src="../../assets/figures/example_selection_frequency.json"><img src="../../assets/figures/example_selection_frequency.svg" alt="How often each protein was chosen across the 25 folds" loading="lazy"></div>
<figcaption markdown>`plot_selection_frequency(result)`: how often each protein was chosen. IL-17C and OSM in every fold, IL6 in 24 of 25: a stable signature. TWEAK just passes the cut at half of the folds.</figcaption>
</figure>

The synthetic data were built with IL-17C, OSM and IL6 clearly higher in cases, and the
signature finds exactly these. TWEAK is chosen in about half of the folds: a borderline protein,
which [Check the signature](check-signature.md) looks at more closely.

### What is in the result

| Attribute | What it holds |
|---|---|
| `result.performance` | the AUC per model with its interval: the honest estimate |
| `result.proteins` | every protein chosen at least once, with how often and how strongly |
| `result.signature` | the list of signature proteins |
| `result.predictions` | each person's prediction from the folds where they were left out |
| `result.predict(new_table)` | predictions for new people, from a final model fitted on everyone |

??? "Can I keep age and sex in every model?"

    ```python
    result = find_signature(table, "group", positive="case", covariates=["age", "sex"])
    ```

    The covariates go into every model and are never dropped. The report then adds a model on
    the covariates alone: on the example it reaches an AUC of 0.44, and with the proteins 0.99,
    so the proteins carry the prediction. When the two are close, the proteins add little; the
    [PAD000040 example](../examples/pad000040.md#9-find-a-signature) shows such a case.

??? "My outcome has three groups, or is a number"

    The same call works: `task="auto"` sees three or more classes (multiclass, scored by AUC and
    balanced accuracy) or a number with many values (continuous, scored by Spearman r and the
    RMSE). `positive` is needed only for two classes.

??? "Is the AUC from `result.predict` on my own data also honest?"

    Only on new people. The final model has seen everyone in `table`, so its predictions for
    them are too good. The honest estimate is `result.performance`; a second cohort tells you
    more.

??? question "Check yourself: why are the proteins chosen again in every fold?"

    So that the people who score the model never help to choose its proteins. Choosing once on
    everyone leaks the test people into the choice, and the AUC comes out too high.

Every option: [`find_signature`](../functions/signature.md#find_signature) in Functions.

<div class="sg-lessonnav"><a href="../correlate/">&larr; Correlate</a><a href="../check-signature/">Check the signature &rarr;</a></div>
