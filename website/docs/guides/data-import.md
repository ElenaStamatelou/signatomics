---
title: "Data import: read an Olink NPX file into Python"
description: "Read an NPX file (.csv, .txt, .xlsx or .parquet; long or wide) into long format with read_npx."
hide:
  - toc
---

<p class="eyebrow">Preprocessing · step 1</p>

# Data import

<ol class="sg-route" aria-label="Preprocessing route"><li class="done here"><span class="dot"></span>Data import</li><li class="done"><span class="dot"></span><a href="../quality-control/">Quality control</a></li><li class="done opt"><span class="dot"></span><a href="../bridging/">Bridging</a></li><li class="done"><span class="dot"></span><a href="../sample-metadata/">Sample metadata</a></li></ol>

**Use this when** you have an NPX file, the data file that Olink's software exports
(`.csv`, `.txt`, `.xlsx` or `.parquet`), and want it as a table in Python.
{ .lede }

```python
import signatomics.olink as so

data = so.read_npx(so.example_file())      # your own file: so.read_npx("my_npx_file.csv")
```

`so.example_file()` gives the path to a small, **synthetic** NPX file that comes with
signatomics: one Target 96 plate with 40 samples, 8 controls and 92 assays. The values are
synthetic, so you can follow every guide without data of your own. Real data appears in
[Example: a full analysis](example-analysis.md).

<div class="sg-example" markdown>
<p class="tag">On the synthetic example file</p>

```text
Reading report
  file: example_plate1.csv
  file type: .csv
  layout: long format (one row per sample and assay)
  control samples: 8 (marked in SampleType); apply_qc removes them
  result: 4416 rows = (40 samples + 8 controls) x 92 assays
```

</div>

The plate holds 40 samples and 8 controls, each measured on 92 assays: one row per sample and
assay, so 48 × 92 = 4416 rows. This is the **long format**, the layout every step of signatomics
uses. The controls are Olink's own (negative, plate and sample controls); quality control
removes them.

--8<-- "standard-table.html"

??? "My file looks different. What does the report tell me?"

    When a file needs fixing, the report lists each change. Three examples from real NPX files:

    | Report line | Meaning |
    |---|---|
    | `columns renamed to the names signatomics uses: PlateLOD -> LOD` | the detection-limit column had another name; it is renamed `LOD`, so every step finds it |
    | `NPX: 184 cells said No Data, read as missing values` | 184 measurements held the text "No Data"; they become empty values |
    | `MissingFreq (missing frequency): percentages ("12%") read as fractions (0.12)` | the MissingFreq column, the share of samples below the LOD for each assay, was written as a percentage; it becomes a fraction |

??? "My file calls its columns something else, for example `Patient` instead of `SampleID`. What do I do?"

    signatomics looks for Olink's column names (`SampleID`, `Assay`, `NPX`, ...) and recognises
    the usual variants. When a needed column is missing, the message lists your columns and
    suggests a fix to copy:

    ```text
    my_npx_file.csv: the columns signatomics needs are missing.
      Needed: SampleID (the sample ID), Assay or OlinkID (the protein), NPX (the value)
      Your file has: Patient, Protein, Value, Plate
      Did you mean: read_npx("my_npx_file.csv", column_map={'Patient': 'SampleID', 'Protein': 'Assay', 'Value': 'NPX'})
    ```

    Optional columns (LOD, PlateID, QC_Warning) that are missing get a line in the reading
    report, with what each is needed for.

??? "Which NPX files have been tested?"

    Public NPX files from the [PRIDE affinity proteomics archive](https://www.ebi.ac.uk/pride/archive/affinity-proteomics)
    of Target 48, Target 96, Explore 384, Explore 3072, Explore HT and Reveal; long `.csv`
    files (commas or semicolons, decimal points or commas), `.xlsx` and `.parquet`; and older
    wide files from NPX Manager. Where a file deviates (a sample-ID column holding something
    else, numbers damaged by a spreadsheet), the report says so.

??? question "Check yourself: why does each sample get 92 rows?"

    Each measurement carries its own details: its plate, its QC status, and the LOD and
    MissingFreq of that assay on that plate. One row per sample and assay keeps them together.
    The table turns into one row per sample at the end of preprocessing, in the analysis table.

<div class="sg-lessonnav"><span></span><a href="../quality-control/">Quality control &rarr;</a></div>

