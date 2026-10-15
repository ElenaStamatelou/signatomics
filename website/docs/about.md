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

Version 0.0.1, in development: reading, sample matching, quality control and plate bridging
work, with 71 automatic tests. Not yet on the Python package index; see
[Get started](getting-started.md#1-install). Function names may still change before the first
release.

## Roadmap

<ul class="nodes">
  <li class="done">Read Olink exports: long and wide layouts; csv, Excel and parquet</li>
  <li class="done">Attach sample information, with a report of what did not match</li>
  <li class="done">Quality report and quality step, with every removal reported</li>
  <li class="done">Bridging plates and projects, with independent checks</li>
  <li>One row per person, with an optional total-protein correction for samples whose dilution varies (saliva, urine)</li>
  <li>Compare groups: one test chosen in advance, effect sizes with confidence intervals, FDR across all proteins</li>
  <li>Find a biomarker signature: nested cross-validation, stability selection, calibration</li>
  <li>Check the signature: what each protein adds, learning and selection curves</li>
  <li>Adjust for confounders; before/after and repeated measurements</li>
  <li>First release on the Python package index, with a citable DOI</li>
  <li>Later: ELISA and microbiome modules, and combining data types</li>
</ul>

## Principles

- Nothing changes silently: every step that changes data says what it did.
- Look before you act: checks and charts never change data.
- Methods chosen in advance, not after seeing the data.

## Cite

A DOI will come with the first release. Until then, cite the
[GitHub repository](https://github.com/ElenaStamatelou/signatomics). If you use the example data,
also cite the PAD000040 dataset (see [Get started](getting-started.md#2-get-example-data)).

## Licence and trademarks

An open-source licence (MIT) is planned for the first release. signatomics is an independent
project, not affiliated with or endorsed by Olink or Thermo Fisher Scientific. Olink is a
trademark of Thermo Fisher Scientific.
