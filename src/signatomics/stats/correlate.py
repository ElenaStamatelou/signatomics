"""Correlate every protein with a continuous variable (e.g. HbA1c, age, BMI).

Spearman's rank correlation by default: it follows any rising or falling relation and is robust
to a few extreme values. ``method="pearson"`` measures a straight-line relation. The method is
chosen up front, the same for every protein. The FDR (Benjamini-Hochberg) runs over all tested
proteins.

Intervals: Fisher's z transformation; for Spearman with the Bonett-Wright standard error
(Bonett & Wright 2000, Psychometrika 65:23).
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats

from ._common import StatsResult, _fmt_p, _n, _proteins, fdr_bh

METHODS = ("spearman", "pearson")


def _interval(r: float, n: int, method: str, level: float) -> tuple[float, float]:
    if n <= 3 or abs(r) >= 1:
        return np.nan, np.nan
    se = np.sqrt((1 + r ** 2 / 2) / (n - 3)) if method == "spearman" else 1 / np.sqrt(n - 3)
    z, q = np.arctanh(r), stats.norm.ppf(0.5 + level / 2)
    return float(np.tanh(z - q * se)), float(np.tanh(z + q * se))


def correlate(table: pd.DataFrame, variable: str, *, method: str = "spearman", proteins=None,
              alpha: float = 0.05, level: float = 0.95, verbose: bool = True) -> StatsResult:
    """Correlate NPX with a continuous variable, protein by protein, with an FDR over the panel.

    Parameters
    ----------
    table : the analysis table from ``to_wide``, or any table with ``proteins=`` given.
    variable : the column with the continuous variable, e.g. ``"hba1c"``. Samples without a
        value are left out and counted in the report.
    method : ``"spearman"`` (default) or ``"pearson"``.
    proteins : the protein columns; default ``table.attrs["proteins"]`` from ``to_wide``.
    alpha : the FDR level used in the report (default 0.05).
    level : confidence level of the intervals (default 0.95).
    verbose : print the report.

    Returns a ``StatsResult`` with the summary ``lines`` and the tables ``table`` (per protein:
    ``n``, ``r`` with ``ci_low`` / ``ci_high``, ``p``, ``fdr``; sorted by p) and ``not_tested``.
    """
    if method not in METHODS:
        raise ValueError(f"method must be one of: {', '.join(METHODS)}")
    if variable not in table.columns:
        raise KeyError(f"Column '{variable}' not found in the table.")
    prots = _proteins(table, proteins, exclude=[variable])
    x_all = pd.to_numeric(table[variable], errors="coerce")
    if x_all.notna().sum() < 3:
        raise ValueError(f"Column '{variable}' has fewer than 3 numbers.")

    name = "Spearman's rank correlation" if method == "spearman" else "Pearson's correlation"
    lines = [f"test: {name} with {variable}, {level:.0%} interval"]
    lines.append(f"samples: {int(x_all.notna().sum())} with a value"
                 + (f"; left out: {_n(int(x_all.isna().sum()), 'sample')} without a value"
                    if x_all.isna().any() else ""))

    func = stats.spearmanr if method == "spearman" else stats.pearsonr
    rows, skipped = [], []
    for p in prots:
        y = pd.to_numeric(table[p], errors="coerce")
        ok = x_all.notna() & y.notna()
        n = int(ok.sum())
        if n < 4:
            skipped.append({"protein": p, "reason": "fewer than 4 samples with both values"})
            continue
        if np.ptp(y[ok]) == 0 or np.ptp(x_all[ok]) == 0:
            skipped.append({"protein": p, "reason": "no variation"})
            continue
        res = func(x_all[ok], y[ok])
        r = float(res.statistic)
        lo, hi = _interval(r, n, method, level)
        rows.append({"protein": p, "n": n, "r": r, "ci_low": lo, "ci_high": hi,
                     "p": float(res.pvalue)})

    out = pd.DataFrame(rows)
    if out.empty:
        raise ValueError("No protein could be tested.")
    out["fdr"] = fdr_bh(out["p"])
    out = out.sort_values("p", kind="stable").reset_index(drop=True)
    lines.append(f"tested: {_n(len(out), 'protein')}"
                 + (f"; not tested: {len(skipped)} (see .not_tested)" if skipped else ""))
    sig = out[out["fdr"] < alpha]
    if sig.empty:
        lines.append(f"FDR below {alpha}: none (lowest FDR {_fmt_p(out['fdr'].min())})")
    else:
        shown = [f"{r['protein']} (r {r['r']:+.2f}, FDR {_fmt_p(r['fdr'])})" for _, r in sig.head(8).iterrows()]
        more = f", ... ({len(sig)} in total)" if len(sig) > 8 else ""
        lines.append(f"FDR below {alpha}: {_n(len(sig), 'protein')}: " + ", ".join(shown) + more)
    result = StatsResult(title="Correlation", lines=lines, table=out,
                         not_tested=pd.DataFrame(skipped, columns=["protein", "reason"]))
    if verbose:
        print(result)
    return result
