---
title: "Read Olink NPX files in Python (read_npx)"
description: "read_npx reads any Olink export (long or wide NPX Manager layout; csv, Excel, parquet) into one standard table, with a reading report."
---

# Reading Olink files

Turn any Olink export into the [standard table](index.md#the-standard-table).
{ .lede }

## `read_npx`

Read an Olink export into the standard table.

```python
read_npx(path, *, sep=None, column_map=None, verbose=True) -> pandas.DataFrame
```

`read_npx` reads csv, txt and tsv files, Excel (`.xlsx`, `.xls`) and parquet. It handles the
current **long** layout and older **wide** exports in the NPX Manager style. In those, proteins
sit side by side, with a header block (Panel, Assay, Uniprot ID, OlinkID) and a bottom block
(LOD, Missing Data freq., Plate LOD, …). Wide files are converted to the long layout
automatically. The reader uses fixed rules, with no AI and no internet, so the same file always
gives the same table.

**Inputs**

| Argument | Type | Default | Meaning |
|---|---|---|---|
| `path` | str or Path | — | the file to read |
| `sep` | str | `None` | column separator of a text file; detected automatically when `None` |
| `column_map` | dict | `None` | your own renames for unusual files, e.g. `{"Patient": "SampleID"}` |
| `verbose` | bool | `True` | print the reading report |

**Output.** The [standard table](index.md#the-standard-table). What the reader did is stored in
`df.attrs["read_report"]` (a list of lines).

**What it handles for you**

- column names in their different spellings → standard names (`PlateLOD`/`LODNPX` → `LOD`,
  `SampleQC` → `QC_Warning`, `Sample_ID` → `SampleID`, …)
- sample IDs always read as text (`058` stays `058`, `308_BIS` is fine)
- numbers written as text: decimal commas (`0,5`), percentages (`100%` → 1.0), and "No Data",
  "NA", "-", … → empty
- true/false in other languages (`VERO`/`FALSO`, `WAHR`/`FALSCH`, `WAAR`/`ONWAAR`, …) → True/False
- text files in UTF-8 or Latin-1; separators comma, semicolon, tab or `|`
- control samples and control assays are counted in the report (removed later by
  [`apply_qc`](quality.md#apply_qc))

**Errors**

- `FileNotFoundError`: the file does not exist.
- `ValueError`: the file has no sample ID, value or protein column. The message lists the columns
  that *were* found and suggests `column_map`.

**Example**

```python
npx = so.read_npx("inflam_PRA_dt2.raw.csv")
```

```text
Reading report
  file: inflam_PRA_dt2.raw.csv
  file type: text, separated by comma
  layout: long (one row per sample and protein)
  standard column names: PlateLOD -> LOD
  NPX: converted 184 values written as No Data -> empty
  MissingFreq: converted percentages
  result: 8096 rows, 88 sample IDs, 92 proteins, 88 sample measurements, 1 plates
```

To read several plates, read each file and stack them:

```python
npx = pd.concat([so.read_npx("plate1.csv"), so.read_npx("plate2.csv")], ignore_index=True)
```

---

## `tidy_sample_ids`

Make sample IDs comparable between files.

```python
tidy_sample_ids(ids: pandas.Series) -> pandas.Series
```

Removes spaces at the ends, uses capitals, and treats spaces, underscores and hyphens the same.
Leading zeros are kept.

| In | Out |
|---|---|
| `"308 BIS"`, `"308_bis"`, `"308-BIS"` | `"308_BIS"` |
| `" 058 "` | `"058"` |

**Input.** A pandas Series of IDs of any type; values are turned into text.
**Output.** A Series of tidy IDs, with the same length and order.

[`add_sample_info`](samples.md#add_sample_info) and the bridging functions use it internally. It
is also useful on its own, for comparing ID lists.

```python
so.tidy_sample_ids(pd.Series(["308 BIS", "a-1"]))     # -> "308_BIS", "A_1"
```
