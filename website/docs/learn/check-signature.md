---
title: "Check a biomarker signature: protein contributions, learning and selection curves"
description: "Check a signature from Olink data: what each protein adds (protein_contributions), whether more people would help (learning_curve), and how many proteins are enough (selection_curve)."
hide:
  - toc
---

<p class="eyebrow">Machine learning · step 2</p>

# Check the signature

<ol class="sg-route" aria-label="Machine learning route"><li class="done"><span class="dot"></span><a href="../find-signature/">Find a signature</a></li><li class="done here"><span class="dot"></span>Check the signature</li></ol>

**Use this when** you have a signature from [Find a signature](find-signature.md) and want to
know three things: what each protein adds, whether more people would help, and how many
proteins are enough.
{ .lede }

Each check takes the result of `find_signature` and reruns its folds, so each takes a little
while (one to two minutes on the example).

### 1. What each protein adds

```python
from signatomics.signature import protein_contributions, plot_contributions

check = protein_contributions(result)
plot_contributions(check)
```

<div class="sg-example" markdown>
<p class="tag">On the synthetic example files</p>

```text
Protein contributions
  AUC of the full panel (4 proteins, model 'linear'): 1.00
  the panel is fixed on all people, so these values compare proteins and sit above the honest estimate in .performance
  largest drop when left out: IL-17C (+0.041); strongest alone: IL-17C (0.93)
```

```text
protein  alone  without  drop
 IL-17C   0.93     0.96  0.04
    OSM   0.90     0.98  0.02
  TWEAK   0.66     0.99  0.01
    IL6   0.83     0.99  0.01
```

</div>

`alone` is the AUC of one protein by itself; `without` is the AUC of the panel without it; `drop`
is how much the panel loses. IL-17C and OSM are strong alone, and the panel misses them most.
TWEAK is weak alone (0.66) and the panel hardly misses it: it adds little. The panel was fixed
on all people, so these numbers compare proteins with each other; the honest estimate stays the
one from `find_signature`.

<figure class="sg-chart" markdown>
<div class="sg-plot" data-src="../../assets/figures/example_contributions.json"><img src="../../assets/figures/example_contributions.svg" alt="Each protein's AUC alone and the panel's AUC without it" loading="lazy"></div>
<figcaption markdown>`plot_contributions(check)`: each protein alone, and the panel without it.</figcaption>
</figure>

### 2. Would more people help?

```python
from signatomics.signature import learning_curve, plot_curve

plot_curve(learning_curve(result))
```

<div class="sg-example" markdown>
<p class="tag">On the synthetic example files</p>

```text
Learning curve
  AUC with 14 training people: 0.76; with 46: 0.98
  change over the last step: +0.007
```

</div>

The training part of every fold is cut to 30%, 50%, ... of its people, and the proteins are
chosen again each time. The curve rises fast and is flat by about 30 people: on this data, more
people would add little. A curve that still rises at the right edge says a larger study would
give a better signature.

<figure class="sg-chart" markdown>
<div class="sg-plot" data-src="../../assets/figures/example_learning_curve.json"><img src="../../assets/figures/example_learning_curve.svg" alt="Learning curve: the AUC rises with the number of training people, then levels off" loading="lazy"></div>
<figcaption markdown>`plot_curve(learning_curve(result))`: the AUC against the number of people in each training part, with its 95% interval.</figcaption>
</figure>

### 3. How many proteins are enough?

```python
from signatomics.signature import selection_curve

plot_curve(selection_curve(result))
```

<div class="sg-example" markdown>
<p class="tag">On the synthetic example files</p>

```text
Selection curve
  AUC for 1 protein(s): 0.89; for 20: 0.91
  highest: 0.99 with 3 proteins (the intervals overlap; read the shape of the curve)
```

</div>

In every fold the proteins are ranked on the training people and only the strongest 1, 2, 3, ...
are kept. Two or three proteins reach an AUC near 0.99; more proteins add noise and the AUC
slowly falls. Read the shape of the curve: the intervals overlap, so the exact top is chance.

<figure class="sg-chart" markdown>
<div class="sg-plot" data-src="../../assets/figures/example_selection_curve.json"><img src="../../assets/figures/example_selection_curve.svg" alt="Selection curve: the AUC peaks at two to three proteins" loading="lazy"></div>
<figcaption markdown>`plot_curve(selection_curve(result))`: the AUC against the number of proteins kept, with its 95% interval.</figcaption>
</figure>

Together: three proteins, IL-17C, OSM and IL6, carry the signature. TWEAK passed the cut in
half of the folds, adds little, and the selection curve shows that three proteins are enough.

??? question "Check yourself: why does the full panel reach 1.00 here, above the 0.98 of find_signature?"

    The checks fix the four signature proteins, which were found with all people. Those
    people then also score the panel, so it looks better than it is. `find_signature` chose the
    proteins again in every fold: its 0.98 is the honest number.

Every option: [`protein_contributions`, `learning_curve` and `selection_curve`](../functions/signature.md) in Functions.

<div class="sg-lessonnav"><a href="../find-signature/">&larr; Find a signature</a><a href="../../examples/">Examples &rarr;</a></div>
