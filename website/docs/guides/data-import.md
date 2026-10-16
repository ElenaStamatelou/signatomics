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

`so.example_file()` is a small, **generated** NPX file of one plate that comes with signatomics:
the values are made up, so you can follow every guide without data of your own. Real data
appears in [Example: a full analysis](example-analysis.md).

<div class="sg-example" markdown>
<p class="tag">On the generated example file</p>

```text
Reading report
  file: example_plate1.csv
  file type: .csv
  layout: long format (one row per sample and assay)
  result: 960 rows = 40 samples x 24 assays
```

</div>

The file holds 40 samples, each measured on 24 assays: one row per sample and assay, so
40 × 24 = 960 rows. This is the **long format**, the layout every step of signatomics uses.

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

    Long `.csv` files (commas or semicolons, decimal points or commas), `.xlsx` and `.parquet`
    files of Target 96, Target 96 Extended and Reveal, and older wide files from NPX Manager.
    Explore HT and Explore 3072 files should work as well; they are yet to be tested on real data.

??? question "Check yourself: why does each sample get 24 rows?"

    Each measurement carries its own details: its plate, its QC status, and the LOD and
    MissingFreq of that assay on that plate. One row per sample and assay keeps them together.
    The table turns into one row per sample at the end of preprocessing, in the analysis table.

<div class="sg-lessonnav"><span></span><a href="../quality-control/">Quality control &rarr;</a></div>

