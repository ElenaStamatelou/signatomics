---
title: "Check a biomarker signature: protein contributions, learning and selection curves"
description: "Check a signature from Olink data: what each protein adds (protein_contributions), whether more people would help (learning_curve), and how many proteins are enough (selection_curve)."
hide:
  - toc
---

<p class="eyebrow">Machine learning · step 2</p>

# Check the signature

<ol class="sg-route" aria-label="Machine learning route"><li class="done"><span class="dot"></span><a href="../find-signature/">Find a signature</a></li><li class="done here"><span class="dot"></span>Check the signature</li></ol>

Three checks on a signature from [Find a signature](find-signature.md): what each protein adds,
whether more people would help, and how many proteins are enough.

## 1. Run it

```python
from signatomics.signature import (protein_contributions, learning_curve, selection_curve,
                                   plot_contributions, plot_curve)

plot_contributions(protein_contributions(result))
plot_curve(learning_curve(result))
plot_curve(selection_curve(result))
```

Each check reruns the folds of `find_signature`: one to two minutes on the example.

## 2. What you get

<div class="sg-example" markdown>
<p class="tag">Output</p>

```text
--8<-- "outputs/protein_contributions.txt"
```

```text
--8<-- "outputs/learning_curve.txt"
```

```text
--8<-- "outputs/selection_curve.txt"
```

</div>

<figure class="sg-chart" markdown>
<div class="sg-plot" data-src="../../assets/figures/example_contributions.json"><img src="../../assets/figures/example_contributions.svg" alt="Each protein's AUC alone and the panel's AUC without it" loading="lazy"></div>
<figcaption markdown>What each protein adds: its AUC alone (orange) and the panel's AUC without it (blue). Dashed: the full panel.</figcaption>
</figure>

<figure class="sg-chart" markdown>
<div class="sg-plot" data-src="../../assets/figures/example_learning_curve.json"><img src="../../assets/figures/example_learning_curve.svg" alt="Learning curve: the AUC rises with the number of training people, then levels off" loading="lazy"></div>
<figcaption markdown>Learning curve: the AUC against the number of training people, with its 95% interval.</figcaption>
</figure>

<figure class="sg-chart" markdown>
<div class="sg-plot" data-src="../../assets/figures/example_selection_curve.json"><img src="../../assets/figures/example_selection_curve.svg" alt="Selection curve: the AUC peaks at three proteins" loading="lazy"></div>
<figcaption markdown>Selection curve: the AUC against the number of proteins kept, with its 95% interval.</figcaption>
</figure>

## 3. How to read it

- **What each protein adds**: IL-17C and OSM are strong alone, and the panel misses them most.
  TWEAK is weak alone and the panel hardly misses it.
- **Would more people help?** The learning curve is flat at the right (last step +0.007): more
  people would add little.
- **How many proteins are enough?** The selection curve peaks at 3 proteins.
- **Together**: IL-17C, OSM and IL6 carry the signature.

??? "How it works"

    --8<-- "check-signature.html"

The table behind each check and every option: see [`protein_contributions`, `learning_curve` and
`selection_curve`](../functions/signature.md).

<div class="sg-lessonnav"><a href="../find-signature/">&larr; Find a signature</a><a href="../../examples/">Examples: real Olink data &rarr;</a></div>
