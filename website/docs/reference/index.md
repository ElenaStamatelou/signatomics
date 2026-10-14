---
title: "Function reference: Python tools for Olink data"
description: "Every function of signatomics.olink: reading Olink exports, sample matching, quality control, plate bridging and charts."
---

# Function reference

<p class="lede">Every public function in <code>signatomics.olink</code>: what it does, what goes
in, what comes out, the errors it can raise, and an example on real data.</p>

<p class="small">Version 0.0.1 · updated 13 Oct 2026</p>

```python
import signatomics.olink as so      # every function below is available as so.<function>
```

## Functions by step

| Step | Functions |
|---|---|
| [1. Reading](reading.md) | [`read_npx`](reading.md#read_npx), [`tidy_sample_ids`](reading.md#tidy_sample_ids) |
| [2. Sample information](samples.md) | [`add_sample_info`](samples.md#add_sample_info) |
| [3. Quality](quality.md) | [`qc_report`](quality.md#qc_report), [`apply_qc`](quality.md#apply_qc), [`below_lod_by_group`](quality.md#below_lod_by_group), [`sample_qc_table`](quality.md#sample_qc_table) |
| [4. Bridging plates](bridging.md) | [`bridge_check`](bridging.md#bridge_check), [`bridge_factors`](bridging.md#bridge_factors), [`bridge_plates`](bridging.md#bridge_plates), [`bridge_differences`](bridging.md#bridge_differences) |
| [5. Charts](charts.md) | [`plot_below_lod`](charts.md#plot_below_lod), [`plot_npx_strip`](charts.md#plot_npx_strip), [`plot_group_by_plate`](charts.md#plot_group_by_plate), [`plot_sample_qc`](charts.md#plot_sample_qc), [`plot_bridge_shift`](charts.md#plot_bridge_shift), [`plot_bridge_pca`](charts.md#plot_bridge_pca), [`plot_bridge_scatter`](charts.md#plot_bridge_scatter) |

## Rules that hold for every function

- **Nothing is changed silently.** Functions that change data print a report and also store it in
  the result's `df.attrs` (for example `df.attrs["qc_applied"]`). `verbose=False` stops the
  printing, but the report is still stored.
- **"Look" functions never change data.** These are `qc_report`, `bridge_check` and all charts.
  **"Do" functions** (`apply_qc`, `bridge_plates`) return a *new* table and leave the input as it was.
- **Sample IDs can be hidden.** Every function that can list IDs has a `show_ids` argument. With
  `show_ids=False` it gives counts only, which are safe to share.
- **Olink's control samples and control assays** are recognised from `SampleType` / `AssayType`.
  They are left out of all counts, checks and charts, and `apply_qc` removes them.

## The standard table

Every function reads and returns the same kind of table, a pandas `DataFrame` with **one row per
sample measurement and protein** ("long" format). [`read_npx`](reading.md#read_npx) makes it from
any Olink export.

| Column | What it holds | Always there? |
|---|---|---|
| `SampleID` | the sample ID as text, exactly as in the file (spaces at the ends removed) | yes |
| `NPX` | the value (log2 scale). Older exports may hold `Ct` or `Quantified_value` instead | yes (one of the three) |
| `Assay` and/or `OlinkID` | the protein: name (`IL8`) and/or Olink's ID (`OID00471`) | at least one |
| `UniProt` | UniProt ID of the protein | if in the file |
| `Panel`, `Panel_Version` | the Olink panel | if in the file |
| `PlateID` | the plate the sample was measured on | if in the file |
| `Index`, `WellID` | position of the measurement; tells repeats of the same sample apart | if in the file |
| `QC_Warning` | sample quality flag (`Pass` / `Warning`) | if in the file |
| `Assay_Warning` | protein quality flag | if in the file |
| `LOD`, `MaxLOD` | limit of detection (per protein, often per plate), in NPX | if in the file |
| `MissingFreq` | share of samples below LOD, as a number 0–1 | if in the file |
| `SampleType`, `AssayType` | marks control samples / control assays (newer exports) | if in the file |
| `Normalization`, `Unit`, … | any other column in the file is kept unchanged | — |

Columns added by later steps:

| Column | Added by | What it holds |
|---|---|---|
| `sample` | `add_sample_info` | the tidy sample ID used for matching (see `tidy_sample_ids`) |
| *(your sample-list columns)* | `add_sample_info` | e.g. `diabetes`, `age`, `sex` |
| `NPX_unbridged`, `LOD_unbridged` | `bridge_plates` | the values before bridging |
| `Bridge_adjustment` | `bridge_plates` | NPX added to this row by bridging (0 on the reference plate) |

## Words used in this reference

measurement
:   one sample in one well of one plate (for several panels: per panel). A sample measured twice
    has two measurements.

repeat
:   the same sample measured more than once *on the same plate* (for example a rerun).

bridging sample
:   the same sample measured on *different* plates, on purpose, to line up the plates.

below LOD
:   `NPX < LOD`: the signal cannot be told apart from background noise.
