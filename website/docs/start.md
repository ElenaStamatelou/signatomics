---
title: "Start in 5 minutes: Olink data in Python"
description: "Install signatomics and run the four preprocessing steps on synthetic NPX data: data import, quality control, bridging and sample metadata."
hide:
  - navigation
  - toc
---

<p class="eyebrow">Start</p>

# Start in 5 minutes

The four preprocessing steps on synthetic NPX data that comes with signatomics, so you get
exactly the output shown here.
{ .lede }

### 1. Install

```bash
pip install "signatomics[olink] @ git+https://github.com/ElenaStamatelou/signatomics"
```

<p class="small">Needs Python 3.10 or newer. signatomics installs from GitHub until it is on the
Python package index; <code>[olink]</code> adds what the Olink part needs (reading
<code>.xlsx</code> and <code>.parquet</code> files, charts).</p>

### 2. Data import

```python
import signatomics.olink as so

data = so.read_npx(so.example_file())      # your own file: so.read_npx("my_npx_file.csv")
```

`so.example_file()` gives the path to a **synthetic** NPX file that comes with signatomics: one
Target 96 plate with 40 samples, 8 controls and 92 assays. Synthetic means the values are
generated to look like real Olink data, without any person or study behind them.

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

--8<-- "standard-table.html"

### 3. Quality control

```python
so.qc_report(data)                # look: your data stays as it is
so.plot_npx_strip(data)           # every value per assay, with the LOD
clean = so.apply_qc(data)         # act: removes controls and FAIL samples
```

<div class="sg-example" markdown>
<p class="tag">On the synthetic example file</p>

```text
Quality report
  controls (left out of the checks below; apply_qc removes them): 8 control samples (3 NEGATIVE_CONTROL, 3 PLATE_CONTROL, 2 SAMPLE_CONTROL)
  data: 40 samples, 92 assays, 1 plate
  QC status: 38 PASS, 1 WARN, 1 FAIL
  assays below the LOD in more than 75% of the samples: 11 of 92 (a sign of data quality; they stay in)
  outliers: 1 of 40 samples (P07); see plot_sample_qc
```

</div>

<figure class="sg-chart" markdown>
<div class="sg-plot" data-src="../assets/figures/example_npx_strip.json"><img src="../assets/figures/example_npx_strip.svg" alt="Every NPX value per assay as a dot, with the detection limit as a red dot" loading="lazy"></div>
<figcaption markdown>`so.plot_npx_strip(data)`. Each dot is one sample's NPX value for one assay; the red dot is the assay's detection limit (LOD), and hollow dots lie below it. Hover a dot for its value.</figcaption>
</figure>

### 4. Bridging <span class="sg-opt">optional: for samples on several plates</span>

The example has a second plate; 12 samples are on both.

```python
import pandas as pd

both = pd.concat([so.read_npx(so.example_file(1)), so.read_npx(so.example_file(2))],
                 ignore_index=True)
clean = so.apply_qc(both)
combined = so.bridge_plates(clean)     # plate 2 on the scale of plate 1
```

<div class="sg-example" markdown>
<p class="tag">On the synthetic example files</p>

```text
Bridging report
  reference (kept as it is): Plate 1
  Plate 2: 92 assays adjusted by their adjustment factor (median -0.51 NPX, from -1.67 to +1.05), each from 12 bridge samples
  ...
```

</div>

<figure class="sg-chart" markdown>
<div class="sg-plot" data-src="../assets/figures/example_bridge_pca.json"><img src="../assets/figures/example_bridge_pca.svg" alt="PCA of the bridge samples before and after bridging, coloured by plate" loading="lazy"></div>
<figcaption markdown>`so.plot_bridge_pca(clean, combined)`. Each line joins one bridge sample's two measurements: after bridging they lie together and the plates overlap. The point far from the rest is P33, the broken tube.</figcaption>
</figure>

### 5. Sample metadata

```python
data = so.add_sample_info(combined, so.example_samples(), on="sample_id")
```

<div class="sg-example" markdown>
<p class="tag">On the synthetic example files</p>

```text
Matching report
  matched on: SampleID (NPX data) = sample_id (your sample metadata)
  matched: 59 of 59 sample IDs in the NPX data
  missing from the NPX data: P20
```

</div>

Your sample metadata (demographics, clinical data, ...) now sits next to the NPX values.
P20 is missing because quality control removed it (QC status FAIL).

### All of it in one call

```python
data = so.prepare([so.example_file(1), so.example_file(2)], so.example_samples(), on="sample_id")
```

`prepare` runs the four steps with their default choices and prints one short report, with a
list of what to look at yourself (such as the outlier P07). See [In one call](functions/prepare.md).

### Next

- **[Learn](learn/index.md):** each step on its own, and what its output means.
- **[Functions](functions/index.md):** every function with all its options.
- **[Examples](examples/index.md):** the same steps on real Olink studies.
