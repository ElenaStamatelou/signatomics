---
title: "Olink quality control in Python: QC report, LOD, outliers"
description: "Quality control for Olink NPX data in Python: QC warnings, detection limit (LOD), repeats, outliers, with every removal reported."
---

# Quality

The quality step has two halves: first **look** ([`qc_report`](#qc_report) and the
[charts](charts.md)), then **decide** ([`apply_qc`](#apply_qc)).
{ .lede }

## `qc_report`

Look at data quality without changing anything.

```python
qc_report(npx, *, group=None, max_below_lod=0.5, show_ids=True, verbose=True) -> QCReport
```

It checks:
- control samples and control assays;
- rows without a sample ID;
- quality warnings;
- completely empty measurements;
- repeats and bridging samples;
- proteins that are mostly below LOD.

With `group`, it also checks whether a group is mixed up with the plates, and whether
detectability differs between groups.

**Inputs**

| Argument | Type | Default | Meaning |
|---|---|---|---|
| `npx` | DataFrame | — | table from `read_npx` (optionally after `add_sample_info`) |
| `group` | str | `None` | a sample-information column, e.g. `"diabetes"` |
| `max_below_lod` | float 0–1 | `0.5` | proteins below LOD in more than this share are *reported* (not removed) |
| `show_ids` | bool | `True` | list sample IDs; `False` gives counts only |
| `verbose` | bool | `True` | print the report |

**Output.** A `QCReport` object. It prints as a readable report and displays neatly in a
notebook. The tables behind the report are its attributes:

| Attribute | What it holds |
|---|---|
| `.lines` | the report lines (list of text) |
| `.measurements` | one row per measurement: the ID columns, `values`, `missing`, `qc_warning`, and True/False flags `empty`, `warning`, `repeat`, `bridging` |
| `.plates` | one row per plate: `measurements`, `warnings`, `empty` |
| `.proteins` | one row per protein: `below_lod` (share 0–1), `missing` (share empty) |
| `.group_by_plate` | with `group`: number of samples of each group per plate (column `both/several plates` for bridging samples) |
| `.below_lod_by_group` | with `group`: the table from [`below_lod_by_group`](#below_lod_by_group) |

**Warnings in the report**

- **Group mixed up with plate:** a group that was never measured on a plate holding at least
  5 samples. Plate differences could then look like group differences.
- **Detectability differs:** proteins whose share below LOD differs by at least 30 percentage
  points between groups. Keep these: being undetectable in one group can itself be a signal.

Samples without sample information (no match in the sample list) are not counted as a group.

**Errors.** `KeyError` when `group` is not a column (did you run `add_sample_info`?).

**Example**

```python
rep = so.qc_report(data, group="diabetes", show_ids=False)
rep.group_by_plate          # the table behind the group-vs-plate warning
```

```text
Quality report
  178 sample measurements of 143 sample IDs, 92 proteins, 2 plates
  quality warnings: 17 measurements
  WARNING completely empty measurements (no values at all): 2 (IDs hidden)
  on more than one plate (bridging samples): 35
  proteins below the detection limit in more than 50% of measurements: 14 of 92
  WARNING group mixed up with plate: no 'Prediabetes' measured only on Inflammation F.I.F 21/03/2024. ...
```

---

## `apply_qc`

Remove low-quality data. Every choice is a visible argument with a sensible default, and every
removal is reported.

```python
apply_qc(npx, *, drop_controls=True, drop_warnings=True, drop_empty=True, max_below_lod=None,
         repeats="prefer_pass_then_last", drop_samples=None, show_ids=True, verbose=True)
    -> pandas.DataFrame
```

**Inputs**

| Argument | Type | Default | Meaning |
|---|---|---|---|
| `npx` | DataFrame | — | table from `read_npx` / `add_sample_info` |
| `drop_controls` | bool | `True` | remove control samples and control assays |
| `drop_samples` | list or dict | `None` | samples **you** decided to remove (e.g. an outlier). A list removes all measurements of those samples; a dict `{sample ID: PlateID}` removes only the one on that plate. IDs are compared after tidying |
| `drop_empty` | bool | `True` | remove measurements without any value |
| `drop_warnings` | bool | `True` | remove measurements with a quality warning |
| `repeats` | str | `"prefer_pass_then_last"` | for a sample measured more than once on the same plate: keep a passed measurement, and of those the last one (usually the rerun). Other options: `"first"`, `"last"`, `"keep_all"` |
| `max_below_lod` | float 0–1 | `None` | remove proteins below LOD in more than this share. Default: **keep all proteins** |
| `show_ids` | bool | `True` | list sample IDs in the report |
| `verbose` | bool | `True` | print the report |

The steps run in this order: controls → your samples → empty → warnings → repeats → proteins.
**Bridging samples are never removed as repeats**: the same sample on *different* plates is what
bridging needs.

**Output.** The cleaned table, with the same columns. The list of removals is stored in
`df.attrs["qc_applied"]`.

!!! note "Why all proteins are kept by default"
    A protein below LOD in one group and above it in another is informative. Proteins that carry
    no information drop out of the biomarker search anyway. Look at `qc_report(..., group=...)` or
    `plot_below_lod(..., group=...)` before choosing a threshold.

**Errors.** `ValueError` for an unknown `repeats` rule.

**Example**

```python
clean = so.apply_qc(data, show_ids=False)
```

```text
Quality steps applied
  removed 2 measurements (no values): 2 (IDs hidden)
  removed 15 measurements (quality warning): 14 (IDs hidden)
  kept: 161 measurements, 132 sample IDs, 92 proteins
```

A sensitivity check is one argument away: `so.apply_qc(data, drop_warnings=False)`.

---

## `below_lod_by_group`

Share of measurements below the detection limit, per protein and per group.

```python
below_lod_by_group(npx, group) -> pandas.DataFrame
```

**Inputs.** `npx`: a table with a `LOD` column. `group`: a sample-information column.

**Output.** One row per protein (index: `Assay`) and one column per group, holding the share
below LOD (0–1). A `difference` column holds the largest minus the smallest share. Rows are
sorted by `difference`, largest first.

**Errors.** `ValueError` without a `LOD` column; `KeyError` when `group` is not a column.

```python
so.below_lod_by_group(data, "diabetes").head()
```

---

## `sample_qc_table`

The median and spread (IQR) of each measurement's values over all proteins, with outlier flags.
This is the table behind [`plot_sample_qc`](charts.md#plot_sample_qc).

```python
sample_qc_table(npx, *, sd=3.0) -> pandas.DataFrame
```

**Inputs.** `npx`: the standard table. `sd`: how many standard deviations from the mean count as
an outlier (per panel).

**Output.** One row per measurement: the ID columns, `median`, `iqr`, the limits `median_low`,
`median_high`, `iqr_low` and `iqr_high`, and `outlier` (True/False).

An outlier points to a technical problem (dilution, haemolysis, a pipetting error), not to
biology. If you decide to remove one, do it with `apply_qc(drop_samples=...)` so that the
removal is reported.

```python
t = so.sample_qc_table(clean)
t[t["outlier"]]
```
