---
title: "Read an Olink NPX file in Python"
description: "read_npx reads the NPX file that Olink's software exports (.csv, .txt, .xlsx, .parquet), in long or wide format, and always returns one standard long table."
hide:
  - toc
---

<p class="eyebrow">Preprocessing · step 1</p>

# Data import

<ol class="sg-route" aria-label="Preprocessing route"><li class="done here"><span class="dot"></span>Data import</li><li class="done"><span class="dot"></span><a href="../quality-control/">Quality control</a></li><li class="done opt"><span class="dot"></span><a href="../bridging/">Bridging</a></li><li class="done"><span class="dot"></span><a href="../sample-metadata/">Sample metadata</a></li><li class="done"><span class="dot"></span><a href="../analysis-table/">Analysis table</a></li></ol>

`read_npx` reads the NPX file that Olink's software exports and gives back one standard table:
one row per sample and protein.

## 1. Run it

```python
import signatomics.olink as so

data = so.read_npx(so.example_file())      # your own file: so.read_npx("my_npx_file.csv")
```

`so.example_file()` is a synthetic NPX file that comes with signatomics: one Target 96 plate with
40 samples and 8 control samples.

## 2. What you get

<div class="sg-example" markdown>
<p class="tag">Output</p>

```text
--8<-- "outputs/read_plate1.txt"
```

</div>

## 3. How to read it

- **layout**: the file was already in long format, so it was read as it is. A wide file (one row
  per sample) is turned into long format.
- **control samples**: Olink's own control samples, recognised from `SampleType`.
- **result**: 48 samples and controls × 92 proteins = 4416 rows.
- When the file needs fixing (other column names, "No Data" cells, percentages), the report lists
  each change.

??? "How it works"

    --8<-- "standard-table.html"

Other column names, several plates, and the NPX files `read_npx` was tested on: see all options of
[`read_npx`](../functions/reading.md#read_npx).

<div class="sg-lessonnav"><span></span><a href="../quality-control/">Quality control: what to remove &rarr;</a></div>
