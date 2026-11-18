---
title: "Correlate Olink proteins with a continuous variable in Python"
description: "Correlate every protein with a continuous variable (age, BMI, a lab value) with correlate: Spearman by default, with 95% intervals and an FDR over all proteins."
hide:
  - toc
---

<p class="eyebrow">Statistical analysis · step 2</p>

# Correlate

<ol class="sg-route" aria-label="Statistical analysis route"><li class="done"><span class="dot"></span><a href="../compare-groups/">Compare groups</a></li><li class="done here"><span class="dot"></span>Correlate</li><li><span class="dot"></span>Before and after (coming)</li></ol>

`correlate` tells you which proteins rise or fall with a number, such as age, BMI or a lab value.

## 1. Run it

```python
from signatomics.stats import correlate, plot_volcano

result = correlate(table, "age")
plot_volcano(result)
```

## 2. What you get

<div class="sg-example" markdown>
<p class="tag">Output</p>

```text
--8<-- "outputs/correlate.txt"
```

</div>

<figure class="sg-chart" markdown>
<div class="sg-plot" data-src="../../assets/figures/example_volcano_age.json"><img src="../../assets/figures/example_volcano_age.svg" alt="Volcano plot of the correlation with age: no protein passes the FDR" loading="lazy"></div>
<figcaption markdown>Each dot is one protein: its correlation with age against how certain it is. All dots are grey: none passes the FDR.</figcaption>
</figure>

## 3. How to read it

- **r**: the Spearman correlation, from -1 (falls with age) to +1 (rises with age).
- **FDR**: the p-value corrected for testing 92 proteins at once.
- **The result**: no protein is linked to age. In the synthetic data age plays no role.
- `result.table` holds every protein, with its r, 95% interval, p-value and FDR.

??? "How it works"

    --8<-- "correlate.html"

Pearson instead of Spearman, samples without a value, and the columns of `result.table`: see all
options of [`correlate`](../functions/stats.md#correlate).

<div class="sg-lessonnav"><a href="../compare-groups/">&larr; Compare groups</a><a href="../find-signature/">Find a signature: predict the outcome &rarr;</a></div>
