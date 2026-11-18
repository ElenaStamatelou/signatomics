---
title: "Bridge Olink plates in Python: bridge_check and bridge_plates"
description: "Put Olink plates on one scale with bridge samples: bridge_check shows the plate shift per protein, bridge_plates adjusts each protein by the median difference of its bridge samples."
hide:
  - toc
---

<p class="eyebrow">Preprocessing · step 3 (optional)</p>

# Bridging

<ol class="sg-route" aria-label="Preprocessing route"><li class="done"><span class="dot"></span><a href="../data-import/">Data import</a></li><li class="done"><span class="dot"></span><a href="../quality-control/">Quality control</a></li><li class="done here opt"><span class="dot"></span>Bridging</li><li class="done"><span class="dot"></span><a href="../sample-metadata/">Sample metadata</a></li><li class="done"><span class="dot"></span><a href="../analysis-table/">Analysis table</a></li></ol>

`bridge_plates` puts samples measured on several plates on one scale, using the **bridge samples**:
the samples measured on more than one plate. With one plate, skip to
[Sample metadata](sample-metadata.md).

## 1. Run it

The example comes with a second plate. Read both, clean them, check and bridge:

```python
import pandas as pd

data = pd.concat([so.read_npx(so.example_file(1)),
                  so.read_npx(so.example_file(2))], ignore_index=True)
clean = so.apply_qc(data, drop_samples={"P07": "Plate 1", "P33": "Plate 2"})
so.bridge_check(clean)
combined = so.bridge_plates(clean)
```

P33 is broken on Plate 2 only, so only that measurement is removed.

## 2. What you get

<div class="sg-example" markdown>
<p class="tag">Output of bridge_check</p>

```text
--8<-- "outputs/bridge_check.txt"
```

</div>

<div class="sg-example" markdown>
<p class="tag">Output of bridge_plates</p>

```text
--8<-- "outputs/bridge_plates_drops.txt"
```

</div>

<figure class="sg-chart" markdown>
<div class="sg-plot" data-src="../../assets/figures/example_bridge_pca.json"><img src="../../assets/figures/example_bridge_pca.svg" alt="PCA of the bridge samples before and after bridging, coloured by plate" loading="lazy"></div>
<figcaption markdown>`so.plot_bridge_pca(clean, combined)`. Each line joins one bridge sample's two measurements. Before bridging they lie apart; after bridging the plates overlap.</figcaption>
</figure>

## 3. How to read it

- **Plate 2**: reads higher than Plate 1 by 0.51 NPX (median), by a different amount for each
  protein.
- **shift differently**: P33 moves unlike the other bridge samples. That is why it was removed
  on Plate 2.
- **adjusted**: each protein gets its own adjustment factor, the median of its bridge samples'
  differences, as in Olink's bridge normalization.
- **originals**: the values before bridging stay in the columns `NPX_unbridged` and
  `LOD_unbridged`.

??? "How it works"

    --8<-- "bridging.html"

Choose the reference plate or the bridge samples, bridge projects, save the adjustment factors,
and the other charts: see all options of [`bridge_check` and `bridge_plates`](../functions/bridging.md).

<div class="sg-lessonnav"><a href="../quality-control/">&larr; Quality control</a><a href="../sample-metadata/">Sample metadata: add your clinical data &rarr;</a></div>
