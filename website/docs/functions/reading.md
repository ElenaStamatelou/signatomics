---
title: "Read Olink NPX files in Python (read_npx)"
description: "read_npx reads any NPX file (long or wide NPX Manager layout; .csv, .txt, .xlsx, .parquet) into one standard table, with a reading report."
---

# Data import

Turn any NPX file into the [standard table](index.md#the-standard-table), in long format.
{ .lede }

--8<-- "standard-table.html"

## `read_npx`

Read an NPX file into the standard table.

```python
read_npx(path, *, sep=None, column_map=None, verbose=True) -> pandas.DataFrame
```

`read_npx` reads `.csv`, `.txt` and `.tsv` files, `.xlsx` and `.xls` files, and `.parquet` files. It handles the
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
- `ValueError`: a needed column is missing (sample ID, protein, value). The message lists your
  file's columns and suggests a `column_map` to copy, based on the column names.

The report also warns when a file looks unusual, as seen in public NPX files:

- numbers damaged by a spreadsheet's thousands separator (e.g. `1.205.603`), and NPX values far
  outside the usual range;
- a `SampleID` column that seems to hold something else (each ID has only a few of the assays),
  with the columns that may hold the sample IDs; `column_map={"SampleName": "SampleID"}` then
  keeps the file's own column as `SampleID_in_file`;
- a `Sample_Type` column holding groups (such as `AD`, `HC`) instead of Olink's sample types: it
  is kept as sample information;
- `LOD` or `MissingFreq` columns that are empty in the file: they are left out.

Missing optional columns (`LOD`, `PlateID`, `QC_Warning`) get a line in the reading report, with
what each is needed for and, when a column name hints at one, a `column_map` suggestion.

**Example**

```python
npx = so.read_npx("inflam_PRA_dt2.raw.csv")
```

```text
Reading report
  file: inflam_PRA_dt2.raw.csv
  file type: .csv
  layout: long format (one row per sample and assay)
  columns renamed to the names signatomics uses: PlateLOD -> LOD
  NPX: 184 cells said No Data, read as missing values
  MissingFreq (missing frequency): percentages ("12%") read as fractions (0.12)
  result: 8096 rows = 88 samples x 92 assays
```

To read several plates, read each file and stack them:

```python
npx = pd.concat([so.read_npx("plate1.csv"), so.read_npx("plate2.csv")],
                ignore_index=True)
```
