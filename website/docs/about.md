---
title: "About: open-source Python tools for Olink data"
description: "signatomics is an independent, open-source family of Python tools for biomarker discovery, starting with Olink (PEA) proteomics."
---

# About

signatomics is a family of open-source Python tools for biomarker discovery from omics data,
written by [Elena Stamatelou](https://elenastamatelou.com). Its first module is for Olink (PEA)
proteomics.
{ .lede }

## Status

**Version 0.0.1, in development.** The Olink extension can read, match, quality-check and bridge
data, and it is tested on public and real study data (70 automatic tests). It is not yet on the
Python package index; see [Getting started](getting-started.md#1-install) for how to install it.
Functions may still change before the first release.

## Roadmap

<ul class="soon">
  <li class="done">Read Olink exports: long and wide layouts; csv, Excel and parquet</li>
  <li class="done">Attach sample information, with a report of what did not match</li>
  <li class="done">Quality report and quality step, with every removal reported</li>
  <li class="done">Bridging plates and projects, with independent checks</li>
  <li class="done">Quality and bridging charts</li>
  <li>One value per person, and comparing groups</li>
  <li>Total-protein correction (for example for oral rinse samples)</li>
  <li>Biomarker signature: nested cross-validation, stability selection, calibration</li>
  <li>Newer formats tested on a second public dataset (Explore HT / Reveal)</li>
  <li>First public release: Python package index and a citable DOI (Zenodo)</li>
  <li>Later: more modules (e.g. ELISA, microbiome) and combining data types</li>
</ul>

## Principles

- **Nothing changes silently.** Every step that changes data says what it did.
- **Look before you decide.** Checks and charts never change data.
- **Honest numbers.** Independent checks, sensitivity checks one argument away, and no hidden
  tuning on the test data.
- **Safe to share.** Sample IDs can be hidden in every report and chart.

## How to cite

A citable version (DOI) will come with the first release. Until then, please cite the GitHub
repository. When you use the example data, cite the PAD000040 dataset and its paper (see
[Getting started](getting-started.md#2-get-example-data)).

## Licence

To be decided before the first release; an open-source licence (MIT) is planned.

## Trademarks

signatomics is an independent project, not affiliated with or endorsed by Olink or Thermo Fisher
Scientific. Olink is a trademark of Thermo Fisher Scientific.
