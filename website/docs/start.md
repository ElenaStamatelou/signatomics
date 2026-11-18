---
title: "Start: your first Olink analysis in Python"
description: "Install signatomics, read an NPX file, take a first look, clean it, add your clinical data and compare groups with a volcano plot, on synthetic NPX data."
hide:
  - navigation
  - toc
---

<p class="eyebrow">Start</p>

# Your first analysis

From an NPX file and your clinical data to a first volcano plot, in six steps. The example uses
synthetic NPX data that comes with signatomics, so you get exactly the output shown here.
{ .lede }

## 1. Install

```bash
pip install signatomics
```

<p class="small">Needs Python 3.10 or newer.</p>

## 2. Read your NPX file

```python
import signatomics.olink as so

data = so.read_npx(so.example_file())      # your own file: so.read_npx("my_npx_file.csv")
```

<div class="sg-example" markdown>
<p class="tag">Output</p>

```text
--8<-- "outputs/read_plate1.txt"
```

</div>

One plate: 40 samples, 8 of Olink's control samples, 92 proteins.
[Data import](learn/data-import.md) explains the file formats.

## 3. Take a first look

```python
so.plot_npx_strip(data)
```

<figure class="sg-chart" markdown>
<div class="sg-plot" data-src="../assets/figures/example_npx_strip.json"><img src="../assets/figures/example_npx_strip.svg" alt="Every NPX value per protein as a dot, with the detection limit as a red dot" loading="lazy"></div>
<figcaption markdown>Each dot is one sample's NPX value for one protein. Red dot: the protein's detection limit (LOD); hollow dots lie below it.</figcaption>
</figure>

Proteins such as IL4, IL5 and IL13 lie mostly below the detection limit: there, the signal cannot
be told apart from background noise.

## 4. Clean

```python
so.qc_report(data)
clean = so.apply_qc(data, drop_samples=["P07"])
```

<div class="sg-example" markdown>
<p class="tag">Output of qc_report</p>

```text
--8<-- "outputs/qc_report_plate1.txt"
```

</div>

`apply_qc` removes Olink's control samples and the FAIL sample P20. `drop_samples` removes what
you decide: here the outlier P07. [Quality control](learn/quality-control.md) explains each line.

## 5. Add your clinical data

Your sample metadata is a table of your own: one row per sample, with its group, age, sex or
other clinical data.

```python
samples = so.example_samples(1)            # your own file: "my_sample_metadata.xlsx"
data = so.add_sample_info(clean, samples, on="sample_id")
```

<div class="sg-example" markdown>
<p class="tag">Output</p>

```text
--8<-- "outputs/start_add_sample_info.txt"
```

</div>

`on` names the column with the sample IDs. Every sample found its row; P07 and P20 were removed
in step 4.

## 6. Compare groups

```python
from signatomics.stats import compare_groups, plot_volcano

table = so.to_wide(data)                   # one row per sample, one column per protein
result = compare_groups(table, "group", reference="control")
plot_volcano(result)
```

<div class="sg-example" markdown>
<p class="tag">Output</p>

```text
--8<-- "outputs/start_compare_groups.txt"
```

</div>

<figure class="sg-chart" markdown>
<div class="sg-plot" data-src="../assets/figures/start_volcano.json"><img src="../assets/figures/start_volcano.svg" alt="Volcano plot of cases against controls on one plate: two proteins above the dotted FDR line on the right" loading="lazy"></div>
<figcaption markdown>Each dot is one protein: its difference (right = higher in cases) against how certain it is. Above the dotted line the FDR is below 0.05.</figcaption>
</figure>

IL-17C and OSM are higher in cases. [Compare groups](learn/compare-groups.md) explains the output
and shows how to adjust for age and sex.

## Next

- **[Learn](learn/index.md):** each step in detail, with a second plate (bridging), correlation
  and machine learning.
- **[Functions](functions/index.md):** every function with all its options.
- **[Examples](examples/index.md):** the same steps on a real Olink study.
