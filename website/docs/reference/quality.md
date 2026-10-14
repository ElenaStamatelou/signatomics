---
title: "Olink quality control in Python: QC report, LOD, outliers"
description: "Quality control for Olink NPX data in Python: QC warnings, detection limit (LOD), repeats, outlier samples, with every removal reported, plus Plotly charts."
---

# Quality

The quality step has two halves: first **look** ([`qc_report`](#qc_report) and the
[charts](#charts-for-this-step)), then **decide** ([`apply_qc`](#apply_qc)).
{ .lede }

## `qc_report`

Look at data quality without changing anything.

```python
qc_report(npx, *, group=None, max_below_lod=0.5, outlier_sd=3.0, show_ids=True, verbose=True)
    -> QCReport
```

It checks:
- control samples and control assays;
- rows without a sample ID;
- quality warnings;
- completely empty measurements;
- repeats and bridging samples;
- proteins that are mostly below LOD;
- sample measurements that look odd: their median or spread over all proteins is far from the
  others (Olink's standard sample check).

With `group`, it also checks whether a group is mixed up with the plates, and whether
detectability differs between groups.

**Inputs**

| Argument | Type | Default | Meaning |
|---|---|---|---|
| `npx` | DataFrame | — | table from `read_npx` (optionally after `add_sample_info`) |
| `group` | str | `None` | a sample-information column, e.g. `"diabetes"` |
| `max_below_lod` | float 0–1 | `0.5` | proteins below LOD in more than this share are *reported* (not removed) |
| `outlier_sd` | float | `3.0` | a measurement whose median or IQR over all proteins is more than this many standard deviations from the mean (per panel) is reported as odd |
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
| `.samples` | one row per measurement: `median` and `iqr` of its values over all proteins, the limits `median_low`, `median_high`, `iqr_low`, `iqr_high`, and `outlier` (True/False). The table behind [`plot_sample_qc`](#plot_sample_qc) |
| `.group_by_plate` | with `group`: number of samples of each group per plate (column `both/several plates` for bridging samples) |
| `.below_lod_by_group` | with `group`: one row per protein (index `Assay`), one column per group with the share below LOD (0–1), and `difference` (largest minus smallest share), sorted by `difference` |

**Warnings in the report**

- **Group mixed up with plate:** a group that was never measured on a plate holding at least
  5 samples. Plate differences could then look like group differences.
- **Detectability differs:** proteins whose share below LOD differs by at least 30 percentage
  points between groups. Keep these: being undetectable in one group can itself be a signal.
- **Odd sample measurements** point to a technical problem (dilution, haemolysis, a pipetting
  error), not to biology. If you decide to remove one, do it with `apply_qc(drop_samples=...)`
  so that the removal is reported.

Samples without sample information (no match in the sample list) are not counted as a group.

**Errors.** `KeyError` when `group` is not a column (did you run `add_sample_info`?).

**Example**

```python
rep = so.qc_report(data, group="diabetes", show_ids=False)
rep.group_by_plate          # the table behind the group-vs-plate warning
rep.samples[rep.samples["outlier"]]
```

```text
Quality report
  178 sample measurements of 143 sample IDs, 92 proteins, 2 plates
  quality warnings: 17 measurements
  WARNING completely empty measurements (no values at all): 2 (IDs hidden)
  on more than one plate (bridging samples): 35
  proteins below the detection limit in more than 50% of measurements: 14 of 92
  sample measurements that look odd (median or spread over all proteins more than 3 SD from
    the others; see plot_sample_qc): 1 (IDs hidden)
  WARNING group mixed up with plate: no 'Prediabetes' measured only on Inflammation F.I.F 21/03/2024. ...
```

The odd measurement is sample 249 on plate 2: about 5 NPX too low on almost every protein,
although it passed Olink's own QC. The dataset's authors removed it as well.

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
| `drop_samples` | list or dict | `None` | samples **you** decided to remove (e.g. an odd sample from `qc_report`). A list removes all measurements of those samples; a dict `{sample ID: PlateID}` removes only the one on that plate. IDs are compared after tidying, as in [`add_sample_info`](samples.md#add_sample_info) |
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

## Charts for this step

Every chart returns a Plotly figure: interactive in a notebook, `fig.write_image("chart.png")` for
a static file, `fig.write_html("chart.html")` for a file you can share. See
[the rules for charts](index.md#rules-that-hold-for-every-function).

### `plot_below_lod`

Horizontal bars showing the % of measurements below the detection limit, per protein.

```python
plot_below_lod(npx, *, group=None, max_below_lod=None, only_with_below=True) -> Figure
```

| Argument | Default | Meaning |
|---|---|---|
| `npx` | — | table with a `LOD` column |
| `group` | `None` | one bar per group (at most 3), proteins sorted by the difference between groups |
| `max_below_lod` | `None` | draw a threshold line; proteins above it in red (what `apply_qc(max_below_lod=...)` would remove) |
| `only_with_below` | `True` | leave out proteins that are never below LOD |

**Errors.** `ValueError` without a `LOD` column, or with more than 3 groups.

### `plot_npx_strip`

Every measurement as a dot, per protein, with groups side by side and the LOD as a dark dash.
Values below LOD are drawn **hollow**.

```python
plot_npx_strip(npx, *, group=None, proteins=None, per_row=30, order="name", show_ids=False, seed=0)
    -> Figure
```

| Argument | Default | Meaning |
|---|---|---|
| `group` | `None` | colour by this column (at most 3 groups) |
| `proteins` | `None` | list of protein names to show; default all |
| `per_row` | `30` | proteins per row; long lists are split over rows |
| `order` | `"name"` | `"name"` (alphabetical) or `"below_lod"` (most below LOD first) |
| `show_ids` | `False` | sample IDs in the hover text |
| `seed` | `0` | the random side-to-side spread of the dots (same seed = same picture) |

**Errors.** `ValueError` when no protein is left to plot, or with more than 3 groups.

### `plot_group_by_plate`

Stacked bars showing how many people of each group were measured on each plate. If a group is
missing from a plate, plate differences can look like group differences.

```python
plot_group_by_plate(npx, group) -> Figure
```

**Input.** `group`: a sample-information column. Uses
[`qc_report(...).group_by_plate`](#qc_report).

### `plot_sample_qc`

Each measurement as a dot: its median (x) against its spread (IQR, y) over all proteins. This is
Olink's standard sample check. Dashed lines mark `sd` standard deviations, and outliers are drawn
as large red rings. Colour and marker shape show the plate; there is one panel per Olink panel
(up to 4).

```python
plot_sample_qc(npx, *, sd=3.0, show_ids=False) -> Figure
```

The numbers behind it: [`qc_report(...).samples`](#qc_report).
