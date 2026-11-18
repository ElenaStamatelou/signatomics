---
title: "to_wide: the analysis table for Olink data"
description: "to_wide turns preprocessed Olink NPX data into the analysis table: one row per sample, one column per protein, sample metadata first."
---

# Analysis table

The last step of preprocessing: one row per sample, one column per protein, ready for statistics and machine learning.
{ .lede }

## `to_wide`

Turn the long NPX table into the analysis table.

```python
to_wide(npx, *, values="NPX", verbose=True) -> pandas.DataFrame
```

| Argument | Type | Default | Meaning |
|---|---|---|---|
| `npx` | DataFrame | — | the preprocessed long table (after `add_sample_info`) |
| `values` | str | `"NPX"` | `"NPX"`: the NPX values; `"below_lod"`: True where the value is below the LOD, for `compare_groups(below_lod=...)` |
| `verbose` | bool | `True` | print the report |

**What it does**

- control samples and control assays are left out;
- a sample measured more than once (a bridge sample, or a repeat) gets the **average** of its
  NPX values, so every sample counts once;
- each protein becomes one column, named after the protein; a protein name used by more than one
  OlinkID (the same protein on several Explore panels) gets the OlinkID added, e.g. `IL6 (OID20101)`;
- the columns of your sample metadata come first.

**Output.** A table with one row per sample (index `sample`, the tidy sample ID). The protein
columns are listed in `table.attrs["proteins"]`, so the statistics functions know which columns
hold NPX values; the report is in `table.attrs["wide_report"]`.

**Errors.** `ValueError` without an `NPX` column, for `values="below_lod"` without a `LOD`
column, or when a sample-metadata column has the name of a protein.

**Example**

```python
clean = so.bridge_plates(so.apply_qc(so.example_data()))
data = so.add_sample_info(clean, so.example_samples(), on="sample_id")
table = so.to_wide(data)
```

```text
--8<-- "outputs/to_wide_no_drops.txt"
```

The lesson: [Analysis table](../learn/analysis-table.md).
