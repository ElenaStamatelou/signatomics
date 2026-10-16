---
title: "Bridging Olink plates in Python"
description: "Put several Olink plates or projects on one scale with bridge samples: bridge_check looks at the plate difference, bridge_plates applies the adjustment factors."
hide:
  - toc
---

<p class="eyebrow">Preprocessing · when your samples are on several plates</p>

# Bridging

<ol class="sg-route" aria-label="Preprocessing route"><li class="done"><span class="dot"></span><a href="../data-import/">Data import</a></li><li class="done"><span class="dot"></span><a href="../quality-control/">Quality control</a></li><li class="done here opt"><span class="dot"></span>Bridging</li><li class="done"><span class="dot"></span><a href="../sample-metadata/">Sample metadata</a></li></ol>

**Use this when** your samples were measured on several plates or in several projects, and
some samples were measured on both: the **bridge samples**. With one plate, go straight to
[Sample metadata](sample-metadata.md).
{ .lede }

--8<-- "bridging.html"

The example comes with a second plate. Read both, apply quality control, then bridge:

```python
import pandas as pd

data = pd.concat([so.read_npx(so.example_file(1)),
                  so.read_npx(so.example_file(2))], ignore_index=True)
clean = so.apply_qc(data)

so.bridge_check(clean)                  # look: how different are the plates?
combined = so.bridge_plates(clean)      # act: put them on one scale
```

<div class="sg-example" markdown>
<p class="tag">On the synthetic example files: the check</p>

```text
Bridging check
  reference: Plate 1
  Plate 2: reads +0.51 NPX compared with Plate 1 (median over assays; per assay -1.05 to +1.67), measured on 12 bridge samples
    assays shifted by more than 0.1 NPX: 83 of 92 - bridging needed
    noise: the bridge samples' differences spread 0.47 NPX per assay (typical IQR); bridging removes the shift; this noise stays
    bridge samples that shift differently from the others (a mix-up or a degraded tube?): P33
```

</div>

<figure class="sg-chart" markdown>
<div class="sg-plot" data-src="../../assets/figures/example_bridge_bland_altman.json"><img src="../../assets/figures/example_bridge_bland_altman.svg" alt="Bland-Altman plots of the bridge samples for three proteins" loading="lazy"></div>
<figcaption markdown>`so.plot_bridge_bland_altman(clean)`: how much the plates differ, before bridging, for three proteins chosen from the data (the largest shift up, about none, the largest shift down). Each point is one bridge sample: the mean of its two measurements (left to right) and the difference, Plate 2 minus Plate 1 (bottom to top). The solid line is the median difference: the shift that bridging removes, different for each protein. Between the dotted lines lie the usual differences (median ± 1.96 robust SD). The point far below in each panel is P33, the broken tube.</figcaption>
</figure>


Plate 2 reads higher than Plate 1, by a different amount for each assay, so each assay gets its
own correction.

<div class="sg-example" markdown>
<p class="tag">On the synthetic example files: the bridging</p>

```text
Bridging report
  reference (kept as it is): Plate 1
  Plate 2: 92 assays adjusted by their adjustment factor (median -0.51 NPX, from -1.67 to +1.05), each from 12 bridge samples
  NPX and LOD shifted together, so 'below the LOD' stays the same; the original values are in NPX_unbridged and LOD_unbridged
  check: plot_bridge_pca(before, after) - each bridge sample's two measurements should now lie together
```

</div>

The **adjustment factor** of an assay is the median of its bridge samples' differences between
the two plates, as in Olink's bridge normalization. The median keeps one broken tube, like P33,
from pulling it.

<figure class="sg-chart" markdown>
<div class="sg-plot" data-src="../../assets/figures/example_bridge_pca.json"><img src="../../assets/figures/example_bridge_pca.svg" alt="PCA of the bridge samples before and after bridging, coloured by plate" loading="lazy"></div>
<figcaption markdown>`so.plot_bridge_pca(clean, combined)`: check that it worked. Each line joins one bridge sample's two measurements. After bridging they lie close together and the plates overlap. The point far from the rest is P33, the broken tube.</figcaption>
</figure>


??? "Which plate is the reference, and can I choose it?"

    By default the plate with the most samples: its values stay as they are and the other
    plates move to its scale. To choose another:

    ```python
    combined = so.bridge_plates(clean, reference="Plate 2")
    ```

??? "How many bridge samples do I need, and can I use only some?"

    Olink recommends 8 to 16 for Target 96, and more for Explore and Reveal. To use only some,
    for example after removing a broken one:

    ```python
    combined = so.bridge_plates(clean, bridges=["P29", "P30", "P31"])
    ```

??? "My data comes from two projects. Can I bridge projects?"

    Yes: name the column that holds the project.

    ```python
    combined = so.bridge_plates(both_projects, batch="Project")
    ```

??? "Can I save the adjustment factors and use them again later?"

    ```python
    factors = so.bridge_check(clean).factors
    factors.to_csv("adjustment_factors.csv", index=False)
    combined = so.bridge_plates(clean, factors=pd.read_csv("adjustment_factors.csv"))
    ```

??? "Is the shift the same in every group?"

    ```python
    so.bridge_check(clean, group="group", samples=so.example_samples(), on="sample_id")
    ```

    The report then shows the shift per group. If one group shifts differently, bridging could
    add or remove a group difference.

??? question "Check yourself: after bridging, bridge_check shows shifts of about 0. Does that prove it worked?"

    The adjustment factors were computed from these same bridge samples, so their shifts are
    about 0 by construction. The independent check is the PCA, which looks at all assays
    together.

<div class="sg-lessonnav"><a href="../quality-control/">&larr; Quality control</a><a href="../sample-metadata/">Sample metadata &rarr;</a></div>

