---
title: "Correlate Olink assays with a continuous variable in Python"
description: "Correlate every assay with a continuous variable (age, BMI, a lab value) with correlate: Spearman by default, with 95% intervals and an FDR over all assays."
hide:
  - toc
---

<p class="eyebrow">Statistical analysis · step 3</p>

# Correlate

<ol class="sg-route" aria-label="Statistical analysis route"><li class="done"><span class="dot"></span><a href="../analysis-table/">Analysis table</a></li><li class="done"><span class="dot"></span><a href="../compare-groups/">Compare groups</a></li><li class="done here"><span class="dot"></span>Correlate</li></ol>

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
<p class="tag">On the synthetic example files, without P07 and P33</p>

```text
Correlation
  test: Spearman's rank correlation with age, 95% interval
  samples: 57 with a value
  tested: 92 proteins
  FDR below 0.05: none (lowest FDR 0.47)
```

</div>

### Reading the result

```text
protein  n      r  ci_low  ci_high      p   fdr
  CXCL9 57 -0.334  -0.552  -0.0735  0.011 0.467
TNFRSF9 57  0.319  0.0573     0.54 0.0155 0.467
   CST5 57  0.315  0.0528    0.537  0.017 0.467
```

| Column | What it means |
|---|---|
| `r` | the correlation, from -1 to +1: positive means the assay rises with age |
| `ci_low`, `ci_high` | the 95% interval for `r` |
| `p`, `fdr` | as in [Compare groups](compare-groups.md): read the FDR, corrected for 92 assays |

In the synthetic data age plays no role, and the result says exactly that: CXCL9 has a p-value
of 0.011, yet its FDR is 0.47. Among 92 assays, a few p-values this small are expected by chance
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

<div class="sg-lessonnav"><a href="../compare-groups/">&larr; Compare groups</a><span>Next: adjust for confounders (coming)</span></div>
