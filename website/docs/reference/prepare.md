---
title: "prepare and example data: Olink data in one call"
description: "prepare runs reading, quality control, bridging and sample metadata with the defaults in one call; example_file, example_data and example_samples give made-up data to try."
---

# In one call

`prepare` runs the usual route with each step's default choices and prints one short report.
`example_file`, `example_data` and `example_samples` give generated example data (made-up values) to try everything on.
{ .lede }

## `prepare`

Read, quality control, bridging and sample metadata in one call.

```python
prepare(data, samples=None, *, on=None, group=None, show_ids=True,
        verbose=True) -> pandas.DataFrame
```

| Argument | Type | Default | Meaning |
|---|---|---|---|
| `data` | path, list of paths, or table(s) | — | an NPX file, several NPX files (e.g. one per plate), or tables from `read_npx` / `example_data` |
| `samples` | path or DataFrame | `None` | your sample metadata (demographics, clinical data, ...); needs `on` |
| `on` | str | `None` | the column in your sample metadata that holds the same IDs as `SampleID` in the NPX data |
| `group` | str | `None` | a column of your sample metadata: also checks whether a group is missing from a plate, and whether the plate shift is the same in every group |
| `show_ids` | bool | `True` | list sample IDs in the report; `False` gives counts only |
| `verbose` | bool | `True` | print the report |

**What it does**, in this order, each with its defaults:

1. reads the files ([`read_npx`](reading.md#read_npx));
2. quality control: looks ([`qc_report`](quality.md#qc_report)), then removes control samples,
   measurements with QC status FAIL (WARN is kept), empty measurements and extra copies of a
   sample measured twice on one plate; all assays are kept ([`apply_qc`](quality.md#apply_qc));
3. with several plates and bridge samples: puts every plate on the scale of the reference plate
   ([`bridge_check`](bridging.md#bridge_check), [`bridge_plates`](bridging.md#bridge_plates));
4. adds your sample metadata, if given ([`add_sample_info`](samples.md#add_sample_info)).

For other choices, run the steps yourself.

**Output.** The preprocessed table. The short report is in `df.attrs["prepare_report"]`, every
step's full report in `df.attrs["prepare_details"]`. Things that need your decision, such as an
outlier or a group missing from a plate, are listed under "to look at": `prepare` leaves them
to you.

**Errors.** `ValueError` when `samples` is given without `on`; `TypeError` for data other than a
path or a table; the errors of the steps it runs.

**Example**

```python
clean = so.prepare([so.example_file(1), so.example_file(2)], so.example_samples(),
                   on="sample_id", group="group")
```

```text
Prepared data
  read 2 files: 60 samples, 24 assays, 2 plates
  quality control: removed 1 measurement (QC status FAIL): P20
  bridging: Plate 2 put on the scale of Plate 1 (24 assays adjusted, 12 bridge samples)
  sample metadata: matched: 59 of 59 sample IDs in the NPX data
  result: 59 samples, 24 assays
  to look at (prepare leaves these to you):
    outliers (median or spread over all assays more than 3 SD from the others; see plot_sample_qc): P07, P33
    bridge samples that shift differently from the others (a mix-up or a degraded tube?): P33
    assays whose shift differs by more than 0.5 NPX between groups: IL4, IL8, IFN-gamma, TGF-alpha - ...
```

---

## `example_file`

The path to a small, generated NPX file that comes with the package (made-up values). Read it like your own file.

```python
example_file(plate=1) -> str
```

| Argument | Type | Default | Meaning |
|---|---|---|---|
| `plate` | int | `1` | `1` or `2`: the example has two plates, one file each |

```python
data = so.read_npx(so.example_file())        # plate 1: 40 samples, 24 assays
plate2 = so.read_npx(so.example_file(2))     # plate 2: 32 samples, 12 of them also on plate 1
```

---

## `example_data`

Both example plates in one table, in the [standard table](index.md#the-standard-table), as
`read_npx` returns them. The values are invented and the same every time.

```python
example_data() -> pandas.DataFrame
```

60 people (30 cases, 30 controls), 24 assays, two plates: P01–P40 on Plate 1, P29–P60 on
Plate 2, so P29–P40 are bridge samples. Plate 2 reads higher, by a different amount per assay.
Built in on purpose, so the checks have something to find:

- Plate 1: QC status WARN for P05 and FAIL for P20; P07 far too low on every assay (an outlier);
- Plate 2: QC status WARN for P41; P33 far too low on every assay (a broken bridge sample);
- three assays mostly below the LOD (IL13, IL4, IL5);
- four assays higher in cases (IL6, CXCL10, HGF, OSM).

---

## `example_samples`

The sample metadata that go with `example_data` and `example_file`: one row per person, with `sample_id`, `group`
(`"case"` / `"control"`), `age` and `sex`. Made-up values.

```python
example_samples() -> pandas.DataFrame
```

```python
data = so.add_sample_info(so.example_data(), so.example_samples(), on="sample_id")
```
