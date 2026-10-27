---
title: "to_wide: the analysis table for Olink data"
description: "to_wide turns preprocessed Olink NPX data into the analysis table: one row per sample, one column per assay, sample metadata first."
---

# Analysis table

The hand-off between preprocessing and statistics: one row per sample, one column per assay.
{ .lede }

## `to_wide`

Turn the long NPX table into the analysis table.

```python
to_wide(npx, *, values="NPX", verbose=True) -> pandas.DataFrame
```

| Argument | Type | Default | Meaning |
|---|---|---|---|
| `npx` | DataFrame | — | the preprocessed long table (from `prepare`, or after `add_sample_info`) |
| `values` | str | `"NPX"` | `"NPX"`: the NPX values; `"below_lod"`: True where the value is below the LOD, for `compare_groups(below_lod=...)` |
| `verbose` | bool | `True` | print the report |

**What it does**

- control samples and control assays are left out;
- a sample measured more than once (a bridge sample, or a repeat) gets the **average** of its
  NPX values, so every sample counts once;
- each assay becomes one column, named after the assay; an assay name used by more than one
  OlinkID (the same assay on several Explore panels) gets the OlinkID added, e.g. `IL6 (OID20101)`;
- the columns of your sample metadata come first.

**Output.** A table with one row per sample (index `sample`, the tidy sample ID). The assay
columns are listed in `table.attrs["proteins"]`, so the statistics functions know which columns
hold NPX values; the report is in `table.attrs["wide_report"]`.

**Errors.** `ValueError` without an `NPX` column, for `values="below_lod"` without a `LOD`
column, or when a sample-metadata column has the name of an assay.

**Example**

```python
table = so.to_wide(so.prepare([so.example_file(1), so.example_file(2)],
                                so.example_samples(), on="sample_id"))
```

```text
Analysis table
  averaged: 12 samples measured more than once (bridge samples or repeats) get the average of their NPX values
  result: 59 samples x 92 assays, plus 4 columns of sample metadata (sample_id, group, age, sex)
```

The lesson: [Analysis table](../learn/analysis-table.md).
