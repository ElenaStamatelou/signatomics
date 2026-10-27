---
title: "The analysis table: Olink data ready for statistics in Python"
description: "Turn preprocessed Olink NPX data into the analysis table with to_matrix: one row per sample, one column per assay, sample metadata first."
hide:
  - toc
---

<p class="eyebrow">Statistical analysis · step 1</p>

# Analysis table

<ol class="sg-route" aria-label="Statistical analysis route"><li class="done here"><span class="dot"></span>Analysis table</li><li class="done"><span class="dot"></span><a href="../compare-groups/">Compare groups</a></li><li class="done"><span class="dot"></span><a href="../correlate/">Correlate</a></li></ol>

**Use this when** your data has gone through preprocessing and you want to analyse it: the
analysis table has one row per sample and one column per assay, with your sample metadata in
front. Every statistics function works on it.
{ .lede }

```python
import signatomics.olink as so

data = so.prepare([so.example_file(1), so.example_file(2)], so.example_samples(), on="sample_id")
table = so.to_matrix(data)
```

<div class="sg-example" markdown>
<p class="tag">On the synthetic example files</p>

```text
Analysis table
  averaged: 12 samples measured more than once (bridge samples or repeats) get the average of their NPX values
  result: 59 samples x 92 assays, plus 4 columns of sample metadata (sample_id, group, age, sex)
```

</div>

The long table becomes wide: each sample is now one row, so every sample counts once in the
statistics. The 12 bridge samples were measured on both plates; after bridging, their two values
are averaged. Controls are left out.

```text
       sample_id    group  age     sex   IL8  VEGFA  CD8A  ...
sample
P01          P01  control   73  female  4.28   7.47  7.16
P02          P02     case   60  female  6.75   8.00  8.20
P03          P03  control   63  female  6.13   7.42  7.73
```

### Remove what quality control flagged

Quality control found two broken samples, P07 and P33: both read about 5 NPX too low on every
assay. `prepare` leaves that decision to you. Before the statistics, remove them:

```python
table = table.drop(["P07", "P33"])
```

[Compare groups](compare-groups.md) shows what happens when they stay: both are controls, so
almost every assay would look higher in cases.

??? "Can I keep P33's good measurement?"

    Yes. P33 is a bridge sample: only its measurement on Plate 2 is broken. Remove just that one
    before bridging, and its Plate 1 value stays:

    ```python
    clean = so.apply_qc(data, drop_samples={"P33": "Plate 2"})
    ```

??? "Can I see which values are below the LOD, per sample?"

    ```python
    below = so.to_matrix(data, values="below_lod")    # True where the value is below the LOD
    ```

    Give it to `compare_groups(..., below_lod=below)` to see the share below the LOD per group.

??? question "Check yourself: why does the analysis table come last, after bridging?"

    Each bridge sample has two measurements, one per plate. Only after bridging are they on the
    same scale, so their average is meaningful.

Every option: [`to_matrix`](../functions/analysis-table.md#to_matrix) in Functions.

<div class="sg-lessonnav"><a href="../sample-metadata/">&larr; Sample metadata</a><a href="../compare-groups/">Compare groups &rarr;</a></div>
