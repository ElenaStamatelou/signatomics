---
title: "Olink proteomics in Python"
description: "signatomics: open-source Python tools for Olink (PEA) proteomics data. Read NPX files, check quality, bridge plates, with every step reported."
hide:
  - navigation
  - toc
---

<div class="sg-hero" markdown>

<p class="eyebrow">Python tools for Olink (PEA) proteomics · in development</p>

<h1>signatomics</h1>

<p class="lede">Read, check and bridge your Olink data in Python, with every step explained.
signatomics reads your exports, checks their quality, lines up the plates and reports
everything it does. Nothing is changed silently.</p>

[Get started](getting-started.md){ .md-button .md-button--primary }
[Function reference](reference/index.md){ .md-button }

</div>

<div class="sg-flow" role="img" aria-label="The signatomics workflow: read, match, check quality and bridge plates are ready; compare groups and find a signature are coming next.">
<svg viewBox="0 0 900 132" xmlns="http://www.w3.org/2000/svg">
  <line class="line done" x1="75" y1="40" x2="525" y2="40"/>
  <line class="line" x1="525" y1="40" x2="825" y2="40"/>
  <g class="done">
    <circle cx="75" cy="40" r="17"/><path class="tick" d="M67 40 l6 6 l10 -12" fill="none" stroke-width="3" stroke-linecap="round" stroke-linejoin="round"/>
    <text class="num" x="75" y="12" text-anchor="middle">1</text>
    <text class="lbl" x="75" y="82" text-anchor="middle">Read</text>
    <text class="fn" x="75" y="100" text-anchor="middle">read_npx</text>
  </g>
  <g class="done">
    <circle cx="225" cy="40" r="17"/><path class="tick" d="M217 40 l6 6 l10 -12" fill="none" stroke-width="3" stroke-linecap="round" stroke-linejoin="round"/>
    <text class="num" x="225" y="12" text-anchor="middle">2</text>
    <text class="lbl" x="225" y="82" text-anchor="middle">Match samples</text>
    <text class="fn" x="225" y="100" text-anchor="middle">add_sample_info</text>
  </g>
  <g class="done">
    <circle cx="375" cy="40" r="17"/><path class="tick" d="M367 40 l6 6 l10 -12" fill="none" stroke-width="3" stroke-linecap="round" stroke-linejoin="round"/>
    <text class="num" x="375" y="12" text-anchor="middle">3</text>
    <text class="lbl" x="375" y="82" text-anchor="middle">Check quality</text>
    <text class="fn" x="375" y="100" text-anchor="middle">qc_report · apply_qc</text>
  </g>
  <g class="done">
    <circle cx="525" cy="40" r="17"/><path class="tick" d="M517 40 l6 6 l10 -12" fill="none" stroke-width="3" stroke-linecap="round" stroke-linejoin="round"/>
    <text class="num" x="525" y="12" text-anchor="middle">4</text>
    <text class="lbl" x="525" y="82" text-anchor="middle">Bridge plates</text>
    <text class="fn" x="525" y="100" text-anchor="middle">bridge_plates</text>
  </g>
  <g class="todo">
    <circle cx="675" cy="40" r="17"/>
    <text class="num" x="675" y="12" text-anchor="middle">5</text>
    <text class="lbl" x="675" y="82" text-anchor="middle">Compare groups</text>
    <text class="fn" x="675" y="100" text-anchor="middle">coming next</text>
  </g>
  <g class="todo">
    <circle cx="825" cy="40" r="17"/>
    <text class="num" x="825" y="12" text-anchor="middle">6</text>
    <text class="lbl" x="825" y="82" text-anchor="middle">Biomarker signature</text>
    <text class="fn" x="825" y="100" text-anchor="middle">coming</text>
  </g>
</svg>
</div>
<p class="sg-flow-caption">Filled circles are steps you can use today; dashed circles are the next ones.</p>

## Why signatomics

<div class="grid cards" markdown>

-   :material-file-table-outline:{ .lg } __Reads what Olink gives you__

    ---

    Current long exports and older wide (NPX Manager) files. Reads csv, Excel and parquet, and
    copes with European decimal commas, "No Data" cells and Italian or Dutch TRUE/FALSE words.

-   :material-eye-check-outline:{ .lg } __Look first, then decide__

    ---

    "Look" functions never change your data. "Do" functions return a new table and print what
    they removed or changed. A sensitivity check takes one argument.

-   :material-swap-horizontal:{ .lg } __Plates that agree__

    ---

    Bridging per protein, with detection limits shifted together with NPX. An independent PCA
    check tells you whether it worked.

-   :material-shield-lock-outline:{ .lg } __Safe to share__

    ---

    Every report and chart can hide sample IDs (`show_ids=False`), so you can send results to a
    colleague without sending identities.

</div>

## The whole route in a few lines

```python
import pandas as pd
import signatomics.olink as so

npx   = pd.concat([so.read_npx("plate1.csv"), so.read_npx("plate2.csv")], ignore_index=True)
data  = so.add_sample_info(npx, "samples.xlsx", on="Sample ID")
so.qc_report(data, group="diabetes")            # look
clean = so.apply_qc(data)                       # decide
so.bridge_check(clean, group="diabetes")        # look at the plate difference
bridged = so.bridge_plates(clean)               # line up the plates
so.plot_bridge_pca(clean, bridged)              # check that it worked
```

[Walk through this on real public data →](getting-started.md)

## One family, Olink first

signatomics is a family of tools for finding **biomarker signatures in omics data**. The first
module, `signatomics.olink`, is for **Olink (PEA) proteomics**: Target 96/48, Explore, Explore HT
and Reveal exports. Later modules will build on it:

| Module | Data | Status |
|---|---|---|
| `signatomics.olink` | Olink (PEA) proteomics | in development |
| core | one standard table, combining data types, a small and honest biomarker signature (nested cross-validation) | planned |
| more modules | e.g. ELISA, microbiome | ideas |

See the [roadmap](about.md#roadmap).

<p class="small">signatomics is an independent project, not affiliated with or endorsed by Olink
or Thermo Fisher Scientific. Olink is a trademark of Thermo Fisher Scientific.</p>
