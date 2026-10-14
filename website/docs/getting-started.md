---
title: "Getting started: Olink NPX data in Python"
description: "Install signatomics and go from two raw Olink plate files to quality-checked, bridged NPX data in about ten lines of Python."
---

# Getting started

<p class="lede">Install signatomics, download a free public Olink dataset and go from two raw
plate files to bridged, quality-checked data in about ten lines. Every output on this page
is real; a few long lines are shortened.</p>

## 1. Install

signatomics needs **Python 3.10 or newer**. It is not on the Python package index yet, so you
install it straight from GitHub:

```bash
pip install "signatomics[olink] @ git+https://github.com/ElenaStamatelou/signatomics"
```

The `[olink]` part also installs what the Olink extension needs: `pyarrow` (parquet files),
`openpyxl` (Excel files) and `plotly` (charts).

??? tip "Working in a separate environment (recommended)"

    A separate environment keeps signatomics and its dependencies away from your other projects:

    ```bash
    python -m venv .venv
    .venv\Scripts\activate          # Windows
    source .venv/bin/activate       # macOS / Linux
    pip install "signatomics[olink] @ git+https://github.com/ElenaStamatelou/signatomics"
    ```

Check that it works:

```python
import signatomics.olink as so
print(so.__all__[:3])     # ['read_npx', 'add_sample_info', 'qc_report']
```

## 2. Get example data

This guide uses **PAD000040**, a public Olink dataset from the PRIDE archive. It is free to use
(CC0). It holds the **Target 96 Inflammation** panel (92 proteins) for **142 people with severe
obesity**, measured on **two plates**, with their diabetes status.

<div class="example-box" markdown>
<p class="tag">Dataset</p>

Milito MM, …, Giraudi PJ. *Inflammatory proteomic heterogeneity beyond glycemia status in severe
obesity.* Int J Mol Sci 2026;27(9):4152. [doi:10.3390/ijms27094152](https://doi.org/10.3390/ijms27094152)

Download from [PRIDE: PAD000040](https://www.ebi.ac.uk/pride/archive/projects/PAD000040):

| File | What it is |
|---|---|
| `inflam_NASH_dt1.raw.csv` | plate 1 (NASH study), raw |
| `inflam_PRA_dt2.raw.csv` | plate 2 (PRAESIIDIUM study), raw |
| `2025.05.08_sample_manifest.xlsx` | the sample list: 142 people with diabetes status |

</div>

Put the three files in one folder and start Python there.

## 3. Read the two plates

```python
import pandas as pd
import signatomics.olink as so

plates = pd.concat([so.read_npx("inflam_NASH_dt1.raw.csv"),
                    so.read_npx("inflam_PRA_dt2.raw.csv")], ignore_index=True)
```

`read_npx` prints a short **reading report** for each file:

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

The report shows what the reader fixed: a column renamed to its standard name, "No Data" cells
turned into empty values, and percentages turned into numbers. The result is the
[standard table](reference/index.md#the-standard-table) with one row per sample and protein.

## 4. Attach the sample list

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

`on="FIF ID"` names the column in the sample list that holds the same IDs as the protein file.
Sample 249 has protein data but no row in the sample list. It is **kept** and reported, not
silently dropped.

## 5. Look at quality, then decide

```python
so.qc_report(data, group="diabetes", show_ids=False)
```

```text
Quality report
  178 sample measurements of 143 sample IDs, 92 proteins, 2 plates
  quality warnings: 17 measurements
  WARNING completely empty measurements (no values at all): 2 (IDs hidden)
  on more than one plate (bridging samples): 35
  proteins below the detection limit in more than 50% of measurements: 14 of 92
  WARNING group mixed up with plate: no 'Prediabetes' measured only on Inflammation F.I.F 21/03/2024.
    Plate differences can look like group differences — bridge/normalise the plates and check
    before comparing groups.
```

The last warning matters for this dataset. People with prediabetes were measured only on plate 2
or on both plates, so a plain plate difference would look like a diabetes difference. That's why
we bridge before we compare. Now remove what should go:

```python
clean = so.apply_qc(data, show_ids=False)
```

```text
Quality steps applied
  removed 2 measurements (no values): 2 (IDs hidden)
  removed 15 measurements (quality warning): 14 (IDs hidden)
  kept: 161 measurements, 132 sample IDs, 92 proteins
```

!!! note "All proteins are kept by default"
    A protein that is undetectable in one group but measurable in another can itself be a
    biomarker. `apply_qc` removes proteins only if you give a threshold, for example
    `max_below_lod=0.75`.

## 6. Line up the plates

First look at how different the plates are:

```python
so.bridge_check(clean, group="diabetes", show_ids=False)
```

```text
Bridging check
  reference: Inflammation F.I.F 21/03/2024
  Inflammation PRAESIIDIUM Giraudi 04_03_2025: 29 bridging samples; reads +0.52 NPX vs the
    reference (median over proteins; per protein -0.96 to +2.30; 1.44 times)
    proteins shifted by more than 0.1 NPX: 83 of 92 - bridging needed
    noise: the bridge differences of one protein typically spread 0.44 NPX (IQR)
    bridging samples that shift differently from the others: 2 (IDs hidden)
  'diabetes' among the bridging samples: No diabetes: 9 of 62; T2DM: 9 of 46; Prediabetes: 10 of 23
    overall shift per group: No diabetes +0.56; Prediabetes +0.57; T2DM +0.53
```

Plate 2 reads higher, by **−0.96 to +2.30 NPX depending on the protein**. A single correction for
the whole plate would not be enough. The shift is the same in every diabetes group, so bridging
will not distort the comparison. Now bridge:

```python
bridged = so.bridge_plates(clean)
so.plot_bridge_pca(clean, bridged)
```

The PCA chart shows the bridging samples before and after. Before bridging, a sample's two
measurements lie **6.5** apart (median distance in the PCA). After bridging they lie **1.8**
apart, and the two plates overlap.

??? check "Check yourself: why is the PCA the *independent* check?"
    `bridge_plates` fits one shift **per protein** on the bridging samples. After bridging, those
    per-protein shifts are about zero by construction, so `bridge_check` will always look
    good. The PCA looks at **all proteins together** and asks a different question: do a sample's
    two measurements now land in the same place?

## What next?

<div class="grid cards" markdown>

-   __[Function reference](reference/index.md)__

    Every function, with its inputs, outputs, errors and an example.

-   __[Tutorials](tutorials/index.md)__

    The full PAD000040 analysis, step by step, with interactive figures.

-   __[Concepts](concepts/index.md)__

    NPX, detection limits and bridging, explained without code.

</div>
