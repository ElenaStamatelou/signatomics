---
title: "Find a biomarker signature in Olink data with Python"
description: "find_signature chooses a small set of proteins that predicts an outcome, with nested cross-validation for an honest estimate of how well it predicts new people."
hide:
  - toc
---

<p class="eyebrow">Machine learning · step 1</p>

# Find a signature

<ol class="sg-route" aria-label="Machine learning route"><li class="done here"><span class="dot"></span>Find a signature</li><li class="done"><span class="dot"></span><a href="../check-signature/">Check the signature</a></li></ol>

`find_signature` finds the few proteins that together predict an outcome, such as case or
control, and measures how well they predict people the model has not seen.

## 1. Run it

```python
from signatomics.signature import find_signature, plot_signature_performance, plot_selection_frequency

result = find_signature(table, "group", positive="case")
plot_signature_performance(result)
plot_selection_frequency(result)
```

`table` is the [analysis table](analysis-table.md). `positive` is the group to predict.

## 2. What you get

<div class="sg-example" markdown>
<p class="tag">Output</p>

```text
--8<-- "outputs/find_signature.txt"
```

</div>

<figure class="sg-chart" markdown>
<div class="sg-plot" data-src="../../assets/figures/example_signature_performance.json"><img src="../../assets/figures/example_signature_performance.svg" alt="ROC curves of logistic regression and random forest, both close to the top-left corner" loading="lazy"></div>
<figcaption markdown>How well the signature predicts: the ROC curve of each model. Line: the mean over the 5 repeats; band: their range.</figcaption>
</figure>

<figure class="sg-chart" markdown>
<div class="sg-plot" data-src="../../assets/figures/example_selection_frequency.json"><img src="../../assets/figures/example_selection_frequency.svg" alt="How often each protein was chosen across the 25 folds; IL-17C, OSM, IL6 and TWEAK in the signature" loading="lazy"></div>
<figcaption markdown>Which proteins were chosen, and how often over the 25 folds. Blue dots: the signature.</figcaption>
</figure>

## 3. How to read it

- **AUC 0.98**: the logistic regression ranks a random case above a random control
  98 times in 100. An AUC of 0.5 is a coin toss; 1 is perfect.
- **Two models, one answer**: logistic regression and random forest reach almost the same AUC.
- **The signature**: IL-17C and OSM are chosen in every fold, IL6 in almost every fold, TWEAK in
  half of them. The synthetic data were built with IL-17C, OSM and IL6 higher in cases.
- **The ROC curve**: the closer to the top-left corner, the better the prediction.

??? "How it works"

    --8<-- "find-signature.html"

Adjust for age and sex, predict three or more groups or a number, or predict new people: see all
options of [`find_signature`](../functions/signature.md#find_signature).

<div class="sg-lessonnav"><a href="../correlate/">&larr; Correlate</a><a href="../check-signature/">Check the signature: what each protein adds &rarr;</a></div>
