---
title: "Add sample metadata (clinical data) to Olink NPX data in Python"
description: "add_sample_info merges your sample metadata (group, age, sex, clinical data) with the NPX data on the sample IDs, and lists every sample without a match."
hide:
  - toc
---

<p class="eyebrow">Preprocessing · step 4</p>

# Sample metadata

<ol class="sg-route" aria-label="Preprocessing route"><li class="done"><span class="dot"></span><a href="../data-import/">Data import</a></li><li class="done"><span class="dot"></span><a href="../quality-control/">Quality control</a></li><li class="done opt"><span class="dot"></span><a href="../bridging/">Bridging</a></li><li class="done here"><span class="dot"></span>Sample metadata</li><li class="done"><span class="dot"></span><a href="../analysis-table/">Analysis table</a></li></ol>

`add_sample_info` puts what you know about each sample (group, age, sex, clinical data) next to
its NPX values. Olink calls this table the *sample manifest*.

## 1. Run it

```python
samples = so.example_samples()          # your own: "my_sample_metadata.xlsx" or a .csv file
data = so.add_sample_info(combined, samples, on="sample_id")
```

`combined` is the cleaned and bridged data from the previous step. `on` names the column of your
table with the sample IDs.

## 2. What you get

<div class="sg-example" markdown>
<p class="tag">Output</p>

```text
--8<-- "outputs/add_sample_info_combined.txt"
```

</div>

## 3. How to read it

- **matched**: every sample in the NPX data found its row in your table. Every column of your
  table is added.
- **not in NPX data**: P07 and P20 are in your table only: quality control removed them.
- IDs match even when written a little differently: spaces, underscores, hyphens and capitals are
  ignored, so `"308 BIS"` matches `"308_bis"`.
- A sample without a row in your table stays in the data, and the report lists it.

??? "How it works"

    --8<-- "sample-metadata.html"

Every option and what the report lists: see [`add_sample_info`](../functions/samples.md#add_sample_info).

<div class="sg-lessonnav"><a href="../bridging/">&larr; Bridging</a><a href="../analysis-table/">Analysis table: ready for statistics &rarr;</a></div>
