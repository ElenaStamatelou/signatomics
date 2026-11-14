---
title: "Example Olink data to download and try: synthetic NPX files"
description: "Synthetic Olink Target 96 data to practise with: two NPX files, sample metadata and a ready analysis table. Built into signatomics and free to download."
hide:
  - toc
---

<p class="eyebrow">Learn</p>

# Example data

A small, synthetic Olink study to practise on: built into signatomics, and here as files you can
download. Every lesson uses it, so you get the same results.
{ .lede }

The values are synthetic: no person, no study. They look like an Olink Target 96 study of
60 people (30 cases, 30 controls), with 92 assays named after the Target 96 Inflammation panel,
measured on two plates. People P29 to P40 are on both plates: the bridge samples.

## Download

| File | What it holds |
|---|---|
| [signatomics_example_plate1.csv](../assets/data/signatomics_example_plate1.csv) | NPX file of Plate 1: P01 to P40 and 8 controls, one row per sample and assay |
| [signatomics_example_plate2.csv](../assets/data/signatomics_example_plate2.csv) | NPX file of Plate 2: P29 to P60 and 8 controls |
| [signatomics_example_samples.xlsx](../assets/data/signatomics_example_samples.xlsx) | sample metadata: `sample_id`, `group` (case / control), `age`, `sex` |
| [signatomics_example_analysis_table.xlsx](../assets/data/signatomics_example_analysis_table.xlsx) | the analysis table, ready for statistics: one row per sample, one column per assay, after quality control and bridging |

The NPX files have the columns of an Olink NPX file, so you can also try them in other software.
The analysis table is the quickest start if you want to practise statistics or machine learning
in any tool.

## In Python

The same data comes with the package:

```python
import signatomics.olink as so

so.example_file(1)          # path to the NPX file of Plate 1 (and example_file(2))
so.example_samples()        # the sample metadata, as a table
so.example_data()           # both plates in one table
```

## What is built in

The data has problems on purpose, so the checks have something to find:

- Plate 1: P05 with QC status WARN, P20 with FAIL, and P07 far too low on every assay (an outlier);
- Plate 2: P41 with QC status WARN, and P33 far too low on every assay (a broken bridge sample);
- Plate 2 reads higher than Plate 1, by a different amount per assay;
- ten assays mostly below the LOD, and a few around it;
- IL-17C detectable mostly in cases;
- higher in cases: IL6, OSM and IL-17C clearly, HGF and CXCL10 a little.

The lessons find each of these. The analysis table to download already has P20 (FAIL), P07 and
P33's Plate 2 measurement removed, the decisions made in [Quality control](quality-control.md).

The example data is free to use for any purpose, in teaching, courses or other software
(CC0 public domain dedication): no permission or citation needed.

<div class="sg-lessonnav"><a href="../">&larr; Learn</a><a href="../data-import/">Data import &rarr;</a></div>
