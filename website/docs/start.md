---
title: "Start in 5 minutes: Olink data in Python"
description: "Install signatomics, read an NPX file, look at its quality and preprocess it in one call: four steps."
hide:
  - navigation
  - toc
---

<p class="eyebrow">Start</p>

# Start in 5 minutes

Four steps on a generated example NPX file that comes with signatomics (made-up values), so you
get exactly the output shown here.
{ .lede }

### 1. Install

```bash
pip install "signatomics[olink] @ git+https://github.com/ElenaStamatelou/signatomics"
```

<p class="small">Needs Python 3.10 or newer. signatomics installs from GitHub until it is on the
Python package index; <code>[olink]</code> adds what the Olink part needs (reading
<code>.xlsx</code> and <code>.parquet</code> files, charts).</p>

### 2. Read an NPX file

```python
import signatomics.olink as so

data = so.read_npx(so.example_file())      # your own file: so.read_npx("my_npx_file.csv")
```

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

### 3. Look at its quality

```python
so.qc_report(data)
```

<div class="sg-example" markdown>
<p class="tag">On the generated example file</p>

```text
Quality report
  data: 40 samples, 24 assays, 1 plate
  QC status: 38 PASS, 1 WARN, 1 FAIL
  assays below the LOD in more than 50% of the measurements: 3 of 24
  outliers (median or spread over all assays more than 3 SD from the others; see plot_sample_qc): P07
```

</div>

This only looks: your data stays as it is.

### 4. Preprocess it in one call

```python
clean = so.prepare(so.example_file())
```

<div class="sg-example" markdown>
<p class="tag">On the generated example file</p>

```text
Prepared data
  read 1 file: 40 samples, 24 assays
  quality control: removed 1 measurement (QC status FAIL): P20
  result: 39 samples, 24 assays
  to look at (prepare leaves these to you):
    outliers (median or spread over all assays more than 3 SD from the others; see plot_sample_qc): P07
```

</div>

`clean` is ready for analysis. `prepare` also bridges plates and adds your sample metadata when you
give them: see [In one call](reference/prepare.md).

### Next

Each step can also be done on its own, with your own choices: follow [the route](guides/index.md).
