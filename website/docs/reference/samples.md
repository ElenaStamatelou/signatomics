# Sample information

Attach your clinical data (diabetes status, age, sex, …) to the protein data, and see exactly
what did not match.
{ .lede }

## `add_sample_info`

Add the columns of a sample list to the protein data, and report what did not match.

```python
add_sample_info(npx, samples, *, on, npx_on="SampleID", verbose=True) -> pandas.DataFrame
```

**Inputs**

| Argument | Type | Default | Meaning |
|---|---|---|---|
| `npx` | DataFrame | — | table from [`read_npx`](reading.md#read_npx) |
| `samples` | DataFrame, or path to a csv / Excel file | — | the sample list, one row per sample |
| `on` | str | — | the column in the sample list holding the same IDs as the protein file (e.g. `"FIF ID"`) |
| `npx_on` | str | `"SampleID"` | the ID column in the protein data |
| `verbose` | bool | `True` | print the matching report |

**Output.** The protein table with:

- every column of the sample list added. Protein rows without a match are **kept**, with empty
  sample information;
- a column `sample` with the tidy ID used for matching (see
  [`tidy_sample_ids`](reading.md#tidy_sample_ids));
- the ending `_samplelist` on sample-list columns whose name already exists in the protein data
  (e.g. `PlateID_samplelist`).

The matching report is stored in `df.attrs["sample_report"]`. It lists:
- how many samples matched, and how many only after tidying the IDs;
- which IDs are only in the protein data (marked WARNING);
- which IDs are only in the sample list.

**Errors**

- `KeyError`: `on` is not a column of the sample list, or `npx_on` is not in the protein data.
- `ValueError`: the sample list has the same ID twice (after tidying). Each sample may appear once.

**Example**

```python
data = so.add_sample_info(plates, "2025.05.08_sample_manifest.xlsx", on="FIF ID")
```

```text
Matching report
  protein data: 143 sample IDs; sample list: 142 IDs (matched on 'SampleID' = 'FIF ID')
  matched: 142
  WARNING in protein data but NOT in the sample list (kept, without sample information): 249
  columns that exist in both files got the ending _samplelist: PlateID
```
