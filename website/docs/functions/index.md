---
title: "Functions: every signatomics function for Olink data"
description: "Every function of signatomics.olink: reading NPX files, adding sample metadata, quality control, plate bridging, charts, and the one-call route."
hide:
  - toc
---

<p class="eyebrow">Functions</p>

# Functions

Every function with all its arguments, its output and its errors: look things up while you
work. To learn what each step does and what its output means, see [Learn](../learn/index.md).
{ .lede }

```python
import signatomics.olink as so      # every function below: so.<function>
```

| Step | Look (your data stays as it is) | Act (returns a new table) | Charts |
|---|---|---|---|
| [In one call](prepare.md) | | [`prepare`](prepare.md#prepare) | |
| [Data import](reading.md) | | [`read_npx`](reading.md#read_npx) | |
| [Quality control](quality.md) | [`qc_report`](quality.md#qc_report) | [`apply_qc`](quality.md#apply_qc) | [`plot_npx_strip`](quality.md#plot_npx_strip), [`plot_sample_qc`](quality.md#plot_sample_qc), [`plot_below_lod`](quality.md#plot_below_lod), [`plot_group_by_plate`](quality.md#plot_group_by_plate) |
| [Bridging](bridging.md) | [`bridge_check`](bridging.md#bridge_check) | [`bridge_plates`](bridging.md#bridge_plates) | [`plot_bridge_shift`](bridging.md#plot_bridge_shift), [`plot_bridge_bland_altman`](bridging.md#plot_bridge_bland_altman), [`plot_bridge_scatter`](bridging.md#plot_bridge_scatter), [`plot_bridge_pca`](bridging.md#plot_bridge_pca) |
| [Sample metadata](samples.md) | | [`add_sample_info`](samples.md#add_sample_info) | |
| [Analysis table](analysis-table.md) | | [`to_matrix`](analysis-table.md#to_matrix) | |
| [Statistics](stats.md) | | [`compare_groups`](stats.md#compare_groups), [`correlate`](stats.md#correlate) | [`plot_volcano`](stats.md#plot_volcano) |
| [Example data](prepare.md#example_file) | | [`example_file`](prepare.md#example_file), [`example_data`](prepare.md#example_data), [`example_samples`](prepare.md#example_samples) | |

--8<-- "look-act.html"

??? "Rules every function follows"

    - **Nothing is changed silently.** Functions that change data print a report and also store it in
      the result's `df.attrs` (for example `df.attrs["qc_applied"]`). `verbose=False` stops the
      printing, but the report is still stored.
    - **"Look" functions never change data.** These are `qc_report`, `bridge_check` and all charts.
      **"Do" functions** (`apply_qc`, `bridge_plates`) return a *new* table and leave the input as it was.
    - **Sample IDs can be hidden.** Every function that can list IDs has a `show_ids` argument. With
      `show_ids=False` it gives counts only, which are safe to share.
    - **Olink's control samples and control assays** are recognised from `SampleType` / `AssayType`.
      They are left out of all counts, checks and charts, and `apply_qc` removes them.
    - **Charts** return a Plotly figure: interactive in a notebook, `fig.write_image("chart.png")` for
      a static file, `fig.write_html("chart.html")` for a file you can share. Colours come from a
      colour-blind-checked palette in a fixed order, and identity is never shown by colour alone
      (legends, labels, hollow versus filled dots). Sample IDs appear in the hover text only
      with `show_ids=True`.
    - **Every chart looks the same way:** a title saying what is shown, a subtitle saying how to
      read it, the legend above the plot, dots for samples, axis titles with their units, and
      Olink's words (assay, bridge samples, LOD). On this website every chart is interactive:
      hover for values, drag to zoom, double-click to reset.

## The standard table

Every function reads and returns the same kind of table, a pandas `DataFrame` with **one row per
sample and assay** (long format). [`read_npx`](reading.md#read_npx) makes it from
any NPX file.

| Column | What it holds | Always there? |
|---|---|---|
| `SampleID` | the sample ID as text, exactly as in the file (spaces at the ends removed) | yes |
| `NPX` | the value (log2 scale). Older exports may hold `Ct` or `Quantified_value` instead | yes (one of the three) |
| `Assay` and/or `OlinkID` | the protein: name (`IL8`) and/or Olink's ID (`OID00471`) | at least one |
| `UniProt` | UniProt ID of the protein | if in the file |
| `Panel`, `Panel_Version` | the Olink panel | if in the file |
| `PlateID` | the plate the sample was measured on | if in the file |
| `Index`, `WellID` | position of the measurement; tells repeats of the same sample apart | if in the file |
| `QC_Warning` | QC status per sample (`PASS` / `WARN` / `FAIL`; older exports `Pass` / `Warning`) | if in the file |
| `Assay_Warning` | protein quality flag | if in the file |
| `LOD`, `MaxLOD` | limit of detection (per protein, often per plate), in NPX | if in the file |
| `MissingFreq` | share of samples below LOD, as a number 0–1 | if in the file |
| `SampleType`, `AssayType` | marks control samples / control assays (newer exports) | if in the file |
| `Normalization`, `Unit`, … | any other column in the file is kept unchanged | — |

Columns added by later steps:

| Column | Added by | What it holds |
|---|---|---|
| `sample` | `add_sample_info` | the tidy sample ID used for matching ("308 BIS" and "308_bis" both become "308_BIS"; see [`add_sample_info`](samples.md#add_sample_info)) |
| *(your sample-list columns)* | `add_sample_info` | e.g. `diabetes`, `age`, `sex` |
| `NPX_unbridged`, `LOD_unbridged` | `bridge_plates` | the values before bridging |
| `Bridge_adjustment` | `bridge_plates` | NPX added to this row by bridging (0 on the reference plate) |
