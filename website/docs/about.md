---
title: "About signatomics: open-source Python tools for Olink data"
description: "signatomics is an independent, open-source family of Python tools for biomarker discovery, starting with Olink (PEA) proteomics."
---

# About

signatomics is a family of open-source Python tools for finding biomarker signatures in omics
data, written by [Elena Stamatelou](https://elenastamatelou.com). The first module,
`signatomics.olink`, is for Olink (PEA) proteomics: Target 96/48, Explore, Explore HT and Reveal.
{ .lede }

## Who it is for

Researchers who analyse Olink data in Python: biomarker and machine-learning studies, cohorts
measured on several plates or projects, and anyone who wants each cleaning step to be visible
and repeatable. If you work in R, Olink's own R package covers much of the cleaning.

## Status

Version 0.1.0, the first release. Preprocessing (data import, quality control, bridging, sample
metadata, the analysis table), statistics (compare groups, correlate, adjusting for covariates)
and the biomarker signature (find it, then check it) work, with 123 automatic tests. `prepare`
does the preprocessing in one call, and built-in synthetic data lets you try everything.
Version 0.x means function names may still change; every change is listed with each release.

Found a bug, or something unclear? [Open an issue on GitHub](https://github.com/ElenaStamatelou/signatomics/issues):
feedback is very welcome.

## Roadmap

<ul class="nodes">
  <li class="done">Read NPX files: long and wide layouts; .csv, .txt, .xlsx and .parquet</li>
  <li class="done">Quality control: a quality check and charts, then your choices applied, with every removal reported</li>
  <li class="done">Bridging plates and projects, with independent checks</li>
  <li class="done">Sample metadata, with a report of every ID without a match</li>
  <li class="done">The usual route in one call (<code>prepare</code>), and built-in example data</li>
  <li class="done">The analysis table: one row per sample, one column per assay (<code>to_wide</code>)</li>
  <li class="done">Compare groups and correlate: one test chosen in advance, effect sizes with confidence intervals, FDR across all assays, volcano plots</li>
  <li class="done">Adjust for confounders: one model per assay, unadjusted and adjusted volcano plots side by side</li>
  <li class="done">Find a biomarker signature: nested cross-validation, how often each protein is chosen, an honest estimate of performance</li>
  <li class="done">Check the signature: what each protein adds, learning and selection curves</li>
  <li>Before/after and repeated measurements</li>
  <li>First release on the Python package index, with a citable DOI</li>
  <li>Later: ELISA and microbiome modules, and combining data types</li>
</ul>

## Principles

- Nothing changes silently: every step that changes data says what it did.
- Look before you act: checks and charts never change data.
- Methods chosen in advance, before seeing the data.

## Cite

Cite the software with its DOI from Zenodo (it will appear here once version 0.1.0 is archived), or
use "Cite this repository" on [GitHub](https://github.com/ElenaStamatelou/signatomics), which gives the
reference in APA and BibTeX. If you use the public PAD000040 data,
also cite the dataset and its paper (see [the PAD000040 example](examples/pad000040.md#1-get-the-data)).

## Licence and trademarks

signatomics is free and open source under the [MIT licence](https://github.com/ElenaStamatelou/signatomics/blob/main/LICENSE):
you may use, change and share it, also commercially, as long as the copyright notice stays with
it. The synthetic example data is free to use for any purpose (CC0). signatomics is an independent
project, not affiliated with or endorsed by Olink or Thermo Fisher Scientific. Olink is a
trademark of Thermo Fisher Scientific.
