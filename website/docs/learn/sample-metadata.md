---
title: "Sample metadata: add it to Olink data in Python"
description: "Join your sample metadata (demographics, clinical data, lab values) to Olink NPX data by sample ID with add_sample_info, with a report of every ID that did not match."
hide:
  - toc
---

<p class="eyebrow">Preprocessing · last step</p>

# Sample metadata

<ol class="sg-route" aria-label="Preprocessing route"><li class="done"><span class="dot"></span><a href="../data-import/">Data import</a></li><li class="done"><span class="dot"></span><a href="../quality-control/">Quality control</a></li><li class="done opt"><span class="dot"></span><a href="../bridging/">Bridging</a></li><li class="done here"><span class="dot"></span>Sample metadata</li></ol>

**Use this when** you have a table with information about each sample (demographics, clinical
data, lab values, visit, ...) and want it next to the NPX values for the analysis. Olink calls
this table the *sample manifest*.
{ .lede }

Here the cleaned and bridged example from the previous steps gets its sample metadata:

```python
samples = so.example_samples()          # your own: "my_sample_data.xlsx" or a .csv file
data = so.add_sample_info(combined, samples, on="sample_id")
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

The NPX data and your sample metadata are joined by **sample ID**: `on` names the column in your
table that holds the same IDs as the `SampleID` column of the NPX data. Every column of your
table is added, for example `group`, `age` and `sex`. P20 is missing from the NPX data because
quality control removed it (QC status FAIL).

??? "My IDs are written differently in the two tables. Do they still match?"

    Spaces, underscores, hyphens and capitals are ignored when matching: `"308 BIS"` matches
    `"308_bis"`. Leading zeros are kept (`"058"` stays `"058"`).

??? "What happens to a sample that is missing from my sample metadata?"

    It stays in the data with empty sample metadata, and the report lists its ID under
    "missing from your sample metadata". You decide what to do with it.

??? "Why is this the last step?"

    Quality control and bridging need only the NPX data. Your sample metadata is for the analysis:
    next comes the analysis table, with one row per sample and your sample metadata beside the NPX
    values. The checks that need groups can look them up earlier: see
    [Quality control](quality-control.md).

??? question "Check yourself: why did 59 of 60 samples match?"

    The example has 60 people, and quality control removed P20 (QC status FAIL) before this
    step. The report shows it under "missing from the NPX data".

Every option: [`add_sample_info`](../functions/samples.md#add_sample_info) in Functions.

<div class="sg-lessonnav"><a href="../bridging/">&larr; Bridging</a><span>Next: the analysis table (coming)</span></div>

