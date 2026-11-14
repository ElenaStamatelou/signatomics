---
title: "Correlate Olink assays with a continuous variable in Python"
description: "Correlate every assay with a continuous variable (age, BMI, a lab value) with correlate: Spearman by default, with 95% intervals and an FDR over all assays."
hide:
  - toc
---

<p class="eyebrow">Statistical analysis · step 2</p>

# Correlate

<ol class="sg-route" aria-label="Statistical analysis route"><li class="done"><span class="dot"></span><a href="../compare-groups/">Compare groups</a></li><li class="done here"><span class="dot"></span>Correlate</li><li><span class="dot"></span>Before and after (coming)</li></ol>

**Use this when** you want to know which assays rise or fall with a continuous variable, such
as age, BMI or a lab value: every assay is tested, and the results are corrected for testing
many assays at once.
{ .lede }

```python
from signatomics.stats import correlate, plot_volcano

result = correlate(table, "age")
plot_volcano(result)
```

<div class="sg-example" markdown>
<p class="tag">On the synthetic example files, P07 and P33 removed in quality control</p>

```text
Correlation
  test: Spearman's rank correlation with age, 95% interval
  samples: 58 with a value
  tested: 92 proteins
  FDR below 0.05: none (lowest FDR 0.5)
```

</div>

### Reading the result

```text
protein  n      r  ci_low  ci_high      p   fdr
  CXCL9 58 -0.334   -0.55  -0.0753 0.0105   0.5
   CST5 58  0.318   0.059    0.538 0.0148   0.5
TNFRSF9 58  0.314  0.0545    0.534 0.0163   0.5
```

| Column | What it means |
|---|---|
| `r` | the correlation, from -1 to +1: positive means the assay rises with age |
| `ci_low`, `ci_high` | the 95% interval for `r` |
| `p`, `fdr` | as in [Compare groups](compare-groups.md): read the FDR, corrected for 92 assays |

In the synthetic data age plays no role, and the result says exactly that: CXCL9 has a p-value
of 0.011, yet its FDR is 0.5. Among 92 assays, a few p-values this small are expected by chance
alone. This is what the FDR protects you from.

<figure class="sg-chart" markdown>
<div class="sg-plot" data-src="../../assets/figures/example_volcano_age.json"><img src="../../assets/figures/example_volcano_age.svg" alt="Volcano plot of the correlation with age: no assay passes the FDR" loading="lazy"></div>
<figcaption markdown>`plot_volcano(result)`: every assay's correlation with age (left to right) and how certain it is (bottom to top). All dots are grey: no assay passes the FDR.</figcaption>
</figure>

??? "Spearman or Pearson?"

    Spearman (the default) follows any relation that keeps rising or falling, and a few extreme
    values hardly move it. Pearson measures a straight-line relation:

    ```python
    result = correlate(table, "age", method="pearson")
    ```

    Choose before you look at the results, and use the same method for every assay.

??? "Some samples have no value for the variable. What happens?"

    They are left out of that analysis and counted in the report ("left out: ... without a
    value").

??? question "Check yourself: why can a p-value of 0.011 still mean nothing?"

    Testing 92 assays is 92 chances. Even without any real relation, about one in a hundred tests
    gives p below 0.01 by chance. The FDR takes the number of tests into account.

Every option: [`correlate` and `plot_volcano`](../functions/stats.md) in Functions.

<div class="sg-lessonnav"><a href="../compare-groups/">&larr; Compare groups</a><a href="../find-signature/">Find a signature &rarr;</a></div>
