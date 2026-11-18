---
title: "The analysis table: Olink data with one row per sample in Python"
description: "to_wide turns the long NPX table into the analysis table: one row per sample, one column per protein, your sample metadata in front. Every statistics and machine learning function uses it."
hide:
  - toc
---

<p class="eyebrow">Preprocessing · step 5</p>

# Analysis table

<ol class="sg-route" aria-label="Preprocessing route"><li class="done"><span class="dot"></span><a href="../data-import/">Data import</a></li><li class="done"><span class="dot"></span><a href="../quality-control/">Quality control</a></li><li class="done opt"><span class="dot"></span><a href="../bridging/">Bridging</a></li><li class="done"><span class="dot"></span><a href="../sample-metadata/">Sample metadata</a></li><li class="done here"><span class="dot"></span>Analysis table</li></ol>

`to_wide` turns your data into the analysis table: one row per sample, one column per protein,
your sample metadata in front. Every statistics and machine learning function uses this table.

## 1. Run it

```python
table = so.to_wide(data)
table.head()
```

## 2. What you get

<div class="sg-example" markdown>
<p class="tag">Output</p>

```text
--8<-- "outputs/to_wide.txt"
```

</div>

```text
--8<-- "outputs/analysis_table_head.txt"
```

## 3. How to read it

- **result**: 58 samples, each one row, so every sample counts once in the statistics.
- **averaged**: the 11 bridge samples were measured on both plates; their two bridged values are
  averaged.
- Olink's control samples are left out.

??? "How it works"

    --8<-- "analysis-table.html"

The table of values below the LOD, and every option: see [`to_wide`](../functions/analysis-table.md#to_wide).

<div class="sg-lessonnav"><a href="../sample-metadata/">&larr; Sample metadata</a><a href="../compare-groups/">Compare groups: which proteins differ &rarr;</a></div>
