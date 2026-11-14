---
title: "The analysis table: Olink data ready for statistics in Python"
description: "Turn preprocessed Olink NPX data into the analysis table with to_wide: one row per sample, one column per assay, sample metadata first."
hide:
  - toc
---

<p class="eyebrow">Preprocessing · last step</p>

# Analysis table

<ol class="sg-route" aria-label="Preprocessing route"><li class="done"><span class="dot"></span><a href="../data-import/">Data import</a></li><li class="done"><span class="dot"></span><a href="../quality-control/">Quality control</a></li><li class="done opt"><span class="dot"></span><a href="../bridging/">Bridging</a></li><li class="done"><span class="dot"></span><a href="../sample-metadata/">Sample metadata</a></li><li class="done here"><span class="dot"></span>Analysis table</li></ol>

**Use this when** your data has gone through preprocessing and you want to analyse it: the
analysis table has one row per sample and one column per assay, with your sample metadata in
front. Every statistics and machine learning function works on it.
{ .lede }

The whole route in a few lines, with the decisions from [quality control](quality-control.md#3-decide-remove-what-quality-control-flagged)
and [bridging](bridging.md): P07 removed, and P33 only on Plate 2.

```python
import signatomics.olink as so

data = so.prepare([so.example_file(1), so.example_file(2)], so.example_samples(), on="sample_id",
                  drop_samples={"P07": "Plate 1", "P33": "Plate 2"})
table = so.to_wide(data)
```

<div class="sg-example" markdown>
<p class="tag">On the synthetic example files</p>

```text
Prepared data
  read 2 files: 60 samples + 8 controls, 92 assays, 2 plates
  your choice: removed 2 measurements: P07, P33
  quality control: removed 16 control samples: 6 NEGATIVE_CONTROL, 6 PLATE_CONTROL, 4 SAMPLE_CONTROL
  quality control: removed 1 measurement (QC status FAIL): P20
  bridging: Plate 2 put on the scale of Plate 1 (92 assays adjusted, 11 bridge samples)
  sample metadata: matched: 58 of 58 sample IDs in the NPX data
  result: 58 samples, 92 assays

Analysis table
  averaged: 11 samples measured more than once (bridge samples or repeats) get the average of their NPX values
  result: 58 samples x 92 assays, plus 4 columns of sample metadata (sample_id, group, age, sex)
```

</div>

The long table becomes wide: each sample is now one row, so every sample counts once in the
statistics. The 11 bridge samples were measured on both plates; after bridging, their two values
are averaged. P33 keeps its one good value, from Plate 1. Controls are left out.

```text
       sample_id    group  age     sex   IL8  VEGFA  CD8A  ...
sample
P01          P01  control   73  female  4.28   7.47  7.16
P02          P02     case   60  female  6.75   8.00  8.20
P03          P03  control   63  female  6.13   7.42  7.73
```

??? "Can I see which values are below the LOD, per sample?"

    ```python
    below = so.to_wide(data, values="below_lod")    # True where the value is below the LOD
    ```

    Give it to `compare_groups(..., below_lod=below)` to see the share below the LOD per group.

??? question "Check yourself: why does the analysis table come last, after bridging?"

    Each bridge sample has two measurements, one per plate. Only after bridging are they on the
    same scale, so their average is meaningful.

Every option: [`to_wide`](../functions/analysis-table.md#to_wide) in Functions.

<div class="sg-lessonnav"><a href="../sample-metadata/">&larr; Sample metadata</a><a href="../compare-groups/">Compare groups &rarr;</a></div>
