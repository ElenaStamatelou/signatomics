"""Statistical analysis on the analysis table (one row per sample, one column per protein).

Make the table with ``signatomics.olink.to_matrix``. Every function tests all proteins with one
method chosen up front, applies the FDR (Benjamini-Hochberg) over all tested proteins, and
reports effect sizes with confidence intervals next to the p-values.

- ``compare_groups``  two groups, three or more, or groups in a natural order
                      (Welch by default; Mann-Whitney / Kruskal-Wallis as an option)
- ``correlate``       a continuous variable (Spearman by default; Pearson as an option)

Both return a ``StatsResult``: the summary lines plus the tables behind them.
``plot_volcano`` draws a result: every protein's effect against its p-value.
Helpers: ``fdr_bh`` (Benjamini-Hochberg) and ``holm``.
"""

from ._common import StatsResult, fdr_bh, holm
from .compare import compare_groups
from .correlate import correlate
from .plots import plot_volcano

__all__ = ["compare_groups", "correlate", "plot_volcano", "StatsResult", "fdr_bh", "holm"]
