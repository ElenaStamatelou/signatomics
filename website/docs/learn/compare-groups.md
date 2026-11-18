---
title: "Compare groups in Olink data with Python"
description: "Test every protein between groups with compare_groups: Welch's t-test by default, effect sizes with 95% intervals, an FDR over all proteins, and adjustment for age, sex or other covariates; plot_volcano shows them all at once."
hide:
  - toc
---

<p class="eyebrow">Statistical analysis · step 1</p>

# Compare groups

<ol class="sg-route" aria-label="Statistical analysis route"><li class="done here"><span class="dot"></span>Compare groups</li><li class="done"><span class="dot"></span><a href="../correlate/">Correlate</a></li><li><span class="dot"></span>Before and after (coming)</li></ol>

`compare_groups` tests every protein between groups, such as cases and controls, and tells you
which ones differ. `plot_volcano` shows all proteins in one chart.

## 1. Run it

```python
from signatomics.stats import compare_groups, plot_volcano

result = compare_groups(table, "group", reference="control")
plot_volcano(result)
```

`table` is the [analysis table](analysis-table.md). `reference` is the group the others are
compared with.

## 2. What you get

<div class="sg-example" markdown>
<p class="tag">Output</p>

```text
--8<-- "outputs/compare_groups.txt"
```

</div>

<figure class="sg-chart" markdown>
<div class="sg-plot" data-src="../../assets/figures/example_volcano.json"><img src="../../assets/figures/example_volcano.svg" alt="Volcano plot: the difference of every protein against how certain it is; three proteins above the dotted FDR line on the right" loading="lazy"></div>
<figcaption markdown>Each dot is one protein: its difference (right = higher in cases) against how certain it is (higher = more certain). Above the dotted line the FDR is below 0.05; hover a dot for its numbers.</figcaption>
</figure>

## 3. How to read it

- **difference**: mean NPX of cases minus mean NPX of controls. NPX is on a log2 scale, so +1
  means twice as much.
- **FDR**: the p-value corrected for testing 92 proteins at once. Below 0.05 counts as a finding.
- **The result**: IL-17C, OSM and IL6 are higher in cases.
- `result.table` holds every protein, with its difference, 95% interval, p-value and FDR.

## 4. Adjust for age and sex

```python
from signatomics.stats import plot_volcano_pair

adjusted = compare_groups(table, "group", reference="control", covariates=["age", "sex"])
plot_volcano_pair(result, adjusted)
```

<figure class="sg-chart" markdown>
<div class="sg-plot" data-src="../../assets/figures/example_volcano_pair.json"><img src="../../assets/figures/example_volcano_pair.svg" alt="Two volcano plots side by side: unadjusted and adjusted for age and sex" loading="lazy"></div>
<figcaption markdown>Unadjusted (left) and adjusted for age and sex (right). Dotted line: FDR 0.05. The same three proteins stand out.</figcaption>
</figure>

Every protein gets the same linear model, `NPX ~ group + age + sex`. Choose the covariates before
looking at the results. On real data adjusting can change the conclusion: see the
[PAD000040 example](../examples/pad000040.md#8-adjust-for-covariates).

??? "How it works"

    --8<-- "compare-groups.html"

Three or more groups, a rank-based test, a trend, and the share below the LOD per group: see all
options of [`compare_groups`](../functions/stats.md#compare_groups).

<div class="sg-lessonnav"><a href="../analysis-table/">&larr; Analysis table</a><a href="../correlate/">Correlate: proteins and a number &rarr;</a></div>
