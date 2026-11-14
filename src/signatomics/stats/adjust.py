"""Group comparison adjusted for covariates: one linear model per protein, the same for all.

``compare_groups(..., covariates=[...])`` fits, for every protein,

    NPX ~ group + covariate_1 + covariate_2 + ...

with robust (HC3) standard errors. You name the covariates; none is dropped or chosen by its
p-value, and every protein gets the same model. Text and categorical covariates become one
indicator per level (the first level is the reference). Samples with a missing covariate are left
out of every protein and counted in the report.

The group term follows the unadjusted choice: two groups give the adjusted difference in NPX
(with interval), three or more groups a robust Wald F-test of "any group differs" plus every pair
of groups, and ``ordered=True`` the adjusted slope per group step.
"""

from __future__ import annotations

from itertools import combinations

import numpy as np
import pandas as pd
from scipy import stats

from ._common import _n


def _design_covariates(data: pd.DataFrame, covariates: list[str]) -> pd.DataFrame:
    """Numeric covariates as they are; text, categorical and True/False ones as indicators."""
    parts = []
    for c in covariates:
        s = data[c]
        if pd.api.types.is_bool_dtype(s):
            parts.append(s.astype(float).rename(c))
        elif pd.api.types.is_numeric_dtype(s):
            parts.append(pd.to_numeric(s).rename(c))
        else:
            parts.append(pd.get_dummies(s.astype("category"), prefix=c, prefix_sep=" = ",
                                        drop_first=True, dtype=float))
    return pd.concat(parts, axis=1)


def _ols_hc3(y: np.ndarray, X: np.ndarray):
    """Coefficients and the HC3 covariance matrix, or None when the model cannot be fitted."""
    n, p = X.shape
    if n <= p or np.ptp(y) == 0:
        return None
    xtx = X.T @ X
    if np.linalg.matrix_rank(xtx) < p:
        return None
    xtx_inv = np.linalg.inv(xtx)
    beta = xtx_inv @ X.T @ y
    e = y - X @ beta
    h = np.einsum("ij,jk,ik->i", X, xtx_inv, X)
    if (h >= 1 - 1e-10).any():
        return None
    meat = X.T @ (X * (e ** 2 / (1 - h) ** 2)[:, None])
    return beta, xtx_inv @ meat @ xtx_inv


def _contrast(beta, cov, c, df, level):
    est = float(c @ beta)
    se = float(np.sqrt(c @ cov @ c))
    if se == 0:
        return None
    q = stats.t.ppf(0.5 + level / 2, df)
    return est, est - q * se, est + q * se, 2 * stats.t.sf(abs(est / se), df)


def adjusted_comparison(data: pd.DataFrame, proteins: list[str], group: str, levels: list,
                        ref, covariates: list[str], *, ordered: bool, level: float):
    """Fit the adjusted model for every protein.

    ``data`` holds only samples with a group and without a missing covariate. Returns the
    per-protein rows, the pairwise rows (three or more groups, unordered) and the skipped
    proteins, in the layout ``compare_groups`` uses.
    """
    k = len(levels)
    g = data[group]
    cov = _design_covariates(data, covariates)
    ncov = cov.shape[1]
    others = [lv for lv in levels if lv != ref]       # one indicator per group but the reference
    if k == 2 or not ordered:
        gcols = np.column_stack([(g == lv).to_numpy(dtype=float) for lv in others])
    else:
        step = {lv: i for i, lv in enumerate(levels)}
        gcols = g.map(step).to_numpy(dtype=float)[:, None]
    ng = gcols.shape[1]
    X = np.column_stack([np.ones(len(data)), gcols, cov.to_numpy(dtype=float)])
    p_all = X.shape[1]

    rows, pairs, skipped = [], [], []
    for prot in proteins:
        y_all = pd.to_numeric(data[prot], errors="coerce").to_numpy(dtype=float)
        ok = ~np.isnan(y_all)
        counts = {lv: int(((g == lv).to_numpy() & ok).sum()) for lv in levels}
        row = {"protein": prot}
        for lv in levels:
            row[f"n_{lv}"] = counts[lv]
        for lv in levels:
            m = (g == lv).to_numpy() & ok
            row[f"mean_{lv}"] = y_all[m].mean() if m.any() else np.nan
        few = [str(lv) for lv in levels if counts[lv] < 2]
        if few:
            skipped.append({"protein": prot, "reason": f"fewer than 2 values in {', '.join(few)}"})
            continue
        fit = _ols_hc3(y_all[ok], X[ok])
        if fit is None:
            skipped.append({"protein": prot, "reason": "model could not be fitted (no variation, "
                                                       "or covariates that cannot be told apart)"})
            continue
        beta, vc = fit
        df = int(ok.sum()) - p_all
        if k == 2 or ordered:
            c = np.zeros(p_all)
            c[1] = 1.0
            res = _contrast(beta, vc, c, df, level)
            if res is None:
                skipped.append({"protein": prot, "reason": "no variation"})
                continue
            key = "difference" if k == 2 else "slope"
            row.update({key: res[0], "ci_low": res[1], "ci_high": res[2], "p": res[3]})
        else:
            idx = list(range(1, 1 + ng))
            b, v = beta[idx], vc[np.ix_(idx, idx)]
            try:
                f = float(b @ np.linalg.solve(v, b)) / ng
            except np.linalg.LinAlgError:
                skipped.append({"protein": prot, "reason": "no variation"})
                continue
            row.update(statistic=f, p=float(stats.f.sf(f, ng, df)))
            col = {lv: 1 + i for i, lv in enumerate(others)}
            for a, b_ in combinations(levels, 2):
                lo, hi = (b_, a) if b_ == ref else (a, b_)     # the reference, else the earlier group, is subtracted
                c = np.zeros(p_all)
                if lo != ref:
                    c[col[lo]] -= 1.0
                c[col[hi]] += 1.0
                r = _contrast(beta, vc, c, df, level)
                pairs.append({"protein": prot, "group": hi, "reference": lo,
                              "difference": r[0] if r else np.nan,
                              "ci_low": r[1] if r else np.nan, "ci_high": r[2] if r else np.nan,
                              "p": r[3] if r else np.nan})
        rows.append(row)
    return rows, pairs, skipped, p_all


def covariate_lines(table: pd.DataFrame, group: str, covariates: list[str], used: int,
                    dropped: int, n_terms: int) -> list[str]:
    """The report lines about the adjustment."""
    lines = [f"adjusted for: {', '.join(covariates)} (one model for every protein, robust HC3 "
             f"standard errors; no covariate is chosen or dropped by its p-value)"]
    if dropped:
        lines.append(f"left out: {_n(dropped, 'sample')} with a missing covariate; "
                     f"{used} samples remain")
    if used / n_terms < 10:
        lines.append(f"WARNING: {used} samples for {n_terms} model terms (fewer than 10 per term): "
                     f"the adjusted estimates are unstable; use fewer covariates")
    return lines
