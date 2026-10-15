---
title: "signatomics: Olink proteomics in Python"
description: "signatomics: open-source Python tools for Olink (PEA) proteomics data. Read NPX files, check quality and bridge plates, with every step reported."
hide:
  - navigation
  - toc
---

# signatomics { .sg-title }

Python tools for Olink (PEA) proteomics data. Read your exports, check their quality and line up
the plates, with every step reported. Finding a biomarker signature comes next.
{ .lede }

**For** researchers who analyse Olink data in Python, especially in biomarker and
machine-learning studies, and in studies measured on several plates or projects.
{ .sg-for }

<div class="sg-map">
<svg viewBox="0 0 1000 250" role="group" aria-label="The signatomics functions as a graph. Available: read_npx, add_sample_info, qc_report, apply_qc, bridge_check, bridge_plates. Planned: one row per person, compare groups, find a signature, check the signature, adjust for confounders; later: other data types.">
  <text class="group" x="40" y="20">AVAILABLE · signatomics.olink</text>
  <text class="group" x="590" y="20">PLANNED</text>

  <path class="e" d="M60 130 H195"/>
  <path class="e" d="M195 130 C265 130 265 70 335 70"/>
  <path class="e" d="M195 130 H335"/>
  <path class="e look" d="M335 70 V130"/>
  <path class="e" d="M335 130 C405 130 405 70 475 70"/>
  <path class="e" d="M335 130 H475"/>
  <path class="e look" d="M475 70 V130"/>
  <path class="e plan" d="M475 130 H615"/>
  <path class="e plan" d="M615 205 V130"/>
  <path class="e plan" d="M615 130 C690 130 690 70 765 70"/>
  <path class="e plan" d="M615 130 H765"/>
  <path class="e plan" d="M765 130 C840 130 840 70 915 70"/>
  <path class="e plan" d="M765 130 H915"/>

  <a href="reference/reading/#read_npx"><g class="n done"><title>read_npx: read an Olink export</title><circle cx="60" cy="130" r="7"/><text x="60" y="156" text-anchor="middle">read_npx</text></g></a>
  <a href="reference/samples/#add_sample_info"><g class="n done"><title>add_sample_info: attach your sample list</title><circle cx="195" cy="130" r="7"/><text x="195" y="156" text-anchor="middle">add_sample_info</text></g></a>
  <a href="reference/quality/#qc_report"><g class="n done"><title>qc_report: look at quality (changes nothing)</title><circle cx="335" cy="70" r="7"/><text x="335" y="52" text-anchor="middle">qc_report</text></g></a>
  <a href="reference/quality/#apply_qc"><g class="n done"><title>apply_qc: remove low-quality data</title><circle cx="335" cy="130" r="7"/><text x="335" y="156" text-anchor="middle">apply_qc</text></g></a>
  <a href="reference/bridging/#bridge_check"><g class="n done"><title>bridge_check: look at the plate difference (changes nothing)</title><circle cx="475" cy="70" r="7"/><text x="475" y="52" text-anchor="middle">bridge_check</text></g></a>
  <a href="reference/bridging/#bridge_plates"><g class="n done"><title>bridge_plates: put every plate on one scale</title><circle cx="475" cy="130" r="7"/><text x="475" y="156" text-anchor="middle">bridge_plates</text></g></a>

  <g class="n plan"><title>Planned: one row per person, one column per protein</title><circle cx="615" cy="130" r="7"/><text x="615" y="112" text-anchor="middle">one row per person</text></g>
  <g class="n later"><title>Later: ELISA and microbiome data, ending in the same table</title><circle cx="615" cy="205" r="7"/><text x="615" y="231" text-anchor="middle">other data types (later)</text></g>
  <g class="n plan"><title>Planned: per protein, do the groups differ?</title><circle cx="765" cy="70" r="7"/><text x="765" y="52" text-anchor="middle">compare groups</text></g>
  <g class="n plan"><title>Planned: a small set of proteins that predicts the outcome, by nested cross-validation</title><circle cx="765" cy="130" r="7"/><text x="765" y="156" text-anchor="middle">find a signature</text></g>
  <g class="n plan"><title>Planned: adjust protein–outcome associations for age, sex and other confounders</title><circle cx="915" cy="70" r="7"/><text x="915" y="52" text-anchor="middle">adjust for confounders</text></g>
  <g class="n plan"><title>Planned: what each protein adds, and how stable the signature is</title><circle cx="915" cy="130" r="7"/><text x="915" y="156" text-anchor="middle">check the signature</text></g>
</svg>
</div>
<p class="sg-legend"><span><i class="dot done"></i>available: click for its reference</span><span><i class="dot"></i>planned</span><span>dotted line: a check that informs the next step</span></p>

## In a few lines

```python
import pandas as pd
import signatomics.olink as so

files = ["plate1.csv", "plate2.csv"]
npx   = pd.concat([so.read_npx(f) for f in files], ignore_index=True)
data  = so.add_sample_info(npx, "samples.xlsx", on="Sample ID")
so.qc_report(data, group="diabetes")            # look
clean = so.apply_qc(data)                       # act
so.bridge_check(clean, group="diabetes")        # look
bridged = so.bridge_plates(clean)               # act
so.plot_bridge_pca(clean, bridged)              # check that it worked
```

Every step prints what it found or changed. Checks never change your data; sample IDs can be
hidden in every report (`show_ids=False`).

[Get started](getting-started.md) · [Reference](reference/index.md) · [GitHub](https://github.com/ElenaStamatelou/signatomics)
{ .sg-links }

<p class="small">Version 0.0.1, in development. An independent project, not affiliated with or
endorsed by Olink or Thermo Fisher Scientific. Olink is a trademark of Thermo Fisher Scientific.</p>
