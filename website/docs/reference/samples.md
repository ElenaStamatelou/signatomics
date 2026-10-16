---
title: "Add sample metadata to Olink data (add_sample_info)"
description: "Sample metadata (demographics, clinical data, lab values) to Olink NPX data in Python and get a report of every sample ID without a match."
---

# Sample metadata

Sample metadata (demographics, clinical data, lab values, …) to the NPX data, and see exactly
which IDs have no match.
{ .lede }

## `add_sample_info`

Add the columns of your sample metadata to the NPX data, and report every ID without a match.

```python
add_sample_info(npx, samples, *, on, npx_on="SampleID",
                verbose=True) -> pandas.DataFrame
```

**Inputs**

| Argument | Type | Default | Meaning |
|---|---|---|---|
| `npx` | DataFrame | — | table from [`read_npx`](reading.md#read_npx) |
| `samples` | DataFrame, or path to a `.csv` / `.xlsx` file | — | your sample metadata, one row per sample |
| `on` | str | — | the column in your sample metadata holding the same IDs as the NPX data (e.g. `"FIF ID"`) |
| `npx_on` | str | `"SampleID"` | the ID column in the protein data |
| `verbose` | bool | `True` | print the matching report |

**Output.** The protein table with:

- every column of your sample metadata added. NPX rows without a match are **kept**, with empty
  sample metadata;
- a column `sample` with the tidy ID used for matching: spaces at the ends removed, capitals,
  and spaces, underscores and hyphens treated the same, so `"308 BIS"`, `"308_bis"` and
  `"308-BIS"` all become `"308_BIS"`. Leading zeros are kept (`"058"` stays `"058"`).
  `apply_qc` and the bridging functions match IDs the same way;
- the ending `_samplelist` on sample-list columns whose name already exists in the protein data
  (e.g. `PlateID_samplelist`).

The matching report is stored in `df.attrs["sample_report"]`. It lists:

- how many samples matched, and how many only after tidying the IDs;
- which IDs are only in the protein data (marked WARNING);
- which IDs are only in your sample metadata.

**Errors**

- `KeyError`: `on` is missing from your sample metadata, or `npx_on` from the NPX data.
- `ValueError`: your sample metadata have the same ID twice (after tidying). Each sample may appear once.

**Example**

```python
data = so.add_sample_info(plates, "2025.05.08_sample_manifest.xlsx", on="FIF ID")
```

```text
Matching report
  matched on: SampleID (NPX data) = FIF ID (your sample metadata)
  matched: 142 of 143 sample IDs in the NPX data
  WARNING missing from your sample metadata (kept, without sample metadata): 249
  columns in both tables (yours renamed with the ending _samplelist): PlateID
```
