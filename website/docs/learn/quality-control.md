---
title: "Quality control of Olink data in Python: qc_report and apply_qc"
description: "qc_report shows what is wrong with your Olink data (QC status, proteins below the LOD, outlier samples); apply_qc removes the control samples, FAIL samples and the samples you choose."
hide:
  - toc
---

<p class="eyebrow">Preprocessing · step 2</p>

# Quality control

<ol class="sg-route" aria-label="Preprocessing route"><li class="done"><span class="dot"></span><a href="../data-import/">Data import</a></li><li class="done here"><span class="dot"></span>Quality control</li><li class="done opt"><span class="dot"></span><a href="../bridging/">Bridging</a></li><li class="done"><span class="dot"></span><a href="../sample-metadata/">Sample metadata</a></li><li class="done"><span class="dot"></span><a href="../analysis-table/">Analysis table</a></li></ol>

`qc_report` shows what is wrong with your data. `apply_qc` then removes Olink's control samples,
the samples that failed, and the samples you choose.

## 1. Run it

```python
data = so.read_npx(so.example_file())
so.qc_report(data)
so.plot_sample_qc(data)
```

## 2. What you get

<div class="sg-example" markdown>
<p class="tag">Output</p>

```text
--8<-- "outputs/qc_report_plate1.txt"
```

</div>

<figure class="sg-chart" markdown>
<div class="sg-plot" data-src="../../assets/figures/example_sample_qc.json"><img src="../../assets/figures/example_sample_qc.svg" alt="Each sample's median against its spread over all proteins; one outlier, P07, far to the left" loading="lazy"></div>
<figcaption markdown>Each dot is one sample: its median NPX against its spread. Dashed lines: 3 SD from the others. P07, far to the left, reads low on every protein.</figcaption>
</figure>

## 3. How to read it

- **checked**: the samples of the study. Olink's 8 control samples check the run itself and are
  left out.
- **QC status**: Olink's quality flag per sample. PASS: passed. WARN: one of Olink's internal
  controls lies outside its range. FAIL: failed Olink's quality control.
- **below the LOD**: proteins that mostly measure background noise.
- **outliers**: samples far from the others over all proteins, whatever their QC status. Here P07.

## 4. Clean

```python
clean = so.apply_qc(data, drop_samples=["P07"])
```

<div class="sg-example" markdown>
<p class="tag">Output</p>

```text
--8<-- "outputs/apply_qc_plate1_drop_p07.txt"
```

</div>

`apply_qc` removes the control samples and the FAIL samples, as Olink's own cleaning does.
`drop_samples` removes what you decided: here the outlier P07. WARN samples stay.

<figure class="sg-chart" markdown>
<div class="sg-plot" data-src="../../assets/figures/example_volcano_outliers.json"><img src="../../assets/figures/example_volcano_outliers.svg" alt="Volcano plot with the two broken samples still in: almost every protein leans to the right" loading="lazy"></div>
<figcaption markdown>The [group comparison](compare-groups.md) with P07 and P33 still in: almost every protein looks higher in cases.</figcaption>
</figure>

??? "How it works"

    --8<-- "quality-control.html"

Keep the FAIL samples, remove WARN samples or proteins below the LOD, reports without sample IDs,
and the other charts: see all options of [`qc_report` and `apply_qc`](../functions/quality.md).

<div class="sg-lessonnav"><a href="../data-import/">&larr; Data import</a><a href="../bridging/">Bridging: several plates &rarr;</a></div>
