"""Compare groups, protein by protein: two groups, three or more, or groups in a natural order.

The test is chosen up front, the same for every protein, so all results are comparable:

- ``test="welch"`` (default): Welch's t-test for two groups, Welch's ANOVA for three or more.
  NPX is on a log2 scale, so the difference in mean NPX is the log2 fold change. Welch's tests
  allow each group its own spread, which keeps the error rate right when groups differ in size
  and spread.
- ``test="mann-whitney"``: Mann-Whitney U for two groups, Kruskal-Wallis for three or more. The
  difference reported is the Hodges-Lehmann shift (the median of all differences between a
  person of one group and a person of the other).

With ``ordered=True`` (e.g. normal glucose, prediabetes, diabetes) one trend test replaces the
omnibus test: a straight line through the groups with robust (HC3) standard errors for Welch,
the Jonckheere-Terpstra test for Mann-Whitney.

The FDR (Benjamini-Hochberg) runs over all tested proteins. With three or more groups, every
pair of groups is also compared; those p-values are Holm-adjusted within each protein and are
meant to be read for the proteins whose omnibus (or trend) FDR is below ``alpha``.
"""

from __future__ import annotations

from itertools import combinations

import numpy as np
import pandas as pd
from scipy import stats

from ._common import StatsResult, _fmt_p, _n, _proteins, fdr_bh, holm

TESTS = ("welch", "mann-whitney")


# ---------------------------------------------------------------------------------------------
# The tests (each takes arrays without missing values)
# ---------------------------------------------------------------------------------------------

def _welch_two(ref: np.ndarray, other: np.ndarray, level: float):
    """Difference in means (other - reference), its CI and p (Welch)."""
    n1, n2 = len(ref), len(other)
    v1, v2 = ref.var(ddof=1), other.var(ddof=1)
    se2 = v1 / n1 + v2 / n2
    if se2 == 0:
        return None
    diff = other.mean() - ref.mean()
    se = np.sqrt(se2)
    df = se2 ** 2 / ((v1 / n1) ** 2 / (n1 - 1) + (v2 / n2) ** 2 / (n2 - 1))
    t = diff / se
    q = stats.t.ppf(0.5 + level / 2, df)
    return diff, diff - q * se, diff + q * se, 2 * stats.t.sf(abs(t), df)


def _mann_whitney_two(ref: np.ndarray, other: np.ndarray, level: float):
    """Hodges-Lehmann shift (other - reference), its distribution-free CI, and p (Mann-Whitney)."""
    if np.ptp(np.concatenate([ref, other])) == 0:
        return None
    d = np.sort((other[:, None] - ref[None, :]).ravel())
    n1, n2 = len(ref), len(other)
    z = stats.norm.ppf(0.5 + level / 2)
    k = int(np.floor(n1 * n2 / 2 - z * np.sqrt(n1 * n2 * (n1 + n2 + 1) / 12)))
    k = min(max(k, 0), len(d) - 1)
    p = stats.mannwhitneyu(other, ref, alternative="two-sided").pvalue
    return float(np.median(d)), d[k], d[len(d) - 1 - k], p


def _welch_anova(groups: list[np.ndarray]):
    n = np.array([len(g) for g in groups], dtype=float)
    v = np.array([g.var(ddof=1) for g in groups])
    if (v == 0).any():
        return None
    m = np.array([g.mean() for g in groups])
    k = len(groups)
    w = n / v
    mw = (w * m).sum() / w.sum()
    a = (w * (m - mw) ** 2).sum() / (k - 1)
    lam = ((1 - w / w.sum()) ** 2 / (n - 1)).sum()
    f = a / (1 + 2 * (k - 2) / (k ** 2 - 1) * lam)
    return f, stats.f.sf(f, k - 1, (k ** 2 - 1) / (3 * lam))


def _kruskal(groups: list[np.ndarray]):
    if np.ptp(np.concatenate(groups)) == 0:
        return None
    res = stats.kruskal(*groups)
    return res.statistic, res.pvalue


def _linear_trend(groups: list[np.ndarray], level: float):
    """Slope of NPX over the group steps (0, 1, 2, ...) with HC3 robust standard errors."""
    y = np.concatenate(groups)
    x = np.concatenate([np.full(len(g), i, dtype=float) for i, g in enumerate(groups)])
    if np.ptp(y) == 0:
        return None
    X = np.column_stack([np.ones_like(x), x])
    xtx_inv = np.linalg.inv(X.T @ X)
    beta = xtx_inv @ X.T @ y
    e = y - X @ beta
    h = np.einsum("ij,jk,ik->i", X, xtx_inv, X)
    meat = X.T @ (X * (e ** 2 / (1 - h) ** 2)[:, None])
    se = np.sqrt((xtx_inv @ meat @ xtx_inv)[1, 1])
    if se == 0:
        return None
    df = len(y) - 2
    q = stats.t.ppf(0.5 + level / 2, df)
    slope = beta[1]
    return slope, slope - q * se, slope + q * se, 2 * stats.t.sf(abs(slope / se), df)


def _jonckheere(groups: list[np.ndarray]):
    """Jonckheere-Terpstra trend test (two-sided, normal approximation with tie correction)."""
    y = np.concatenate(groups)
    if np.ptp(y) == 0:
        return None
    jt = 0.0
    for i, j in combinations(range(len(groups)), 2):
        d = groups[j][:, None] - groups[i][None, :]
        jt += (d > 0).sum() + 0.5 * (d == 0).sum()
    n = np.array([len(g) for g in groups], dtype=float)
    big_n = n.sum()
    t = np.unique(y, return_counts=True)[1].astype(float)
    mean = (big_n ** 2 - (n ** 2).sum()) / 4
    var = ((big_n * (big_n - 1) * (2 * big_n + 5) - (n * (n - 1) * (2 * n + 5)).sum()
            - (t * (t - 1) * (2 * t + 5)).sum()) / 72
           + (n * (n - 1) * (n - 2)).sum() * (t * (t - 1) * (t - 2)).sum()
           / (36 * big_n * (big_n - 1) * (big_n - 2))
           + (n * (n - 1)).sum() * (t * (t - 1)).sum() / (8 * big_n * (big_n - 1)))
    z = (jt - mean) / np.sqrt(var)
    return z, 2 * stats.norm.sf(abs(z))


# ---------------------------------------------------------------------------------------------
# compare_groups
# ---------------------------------------------------------------------------------------------

def _levels(values: pd.Series, reference, order) -> list:
    present = list(pd.unique(values))
    if order is not None:
        order = list(order)
        unknown = [g for g in present if g not in order]
        if unknown:
            raise ValueError(f"Groups not in order=: {', '.join(map(str, unknown))}")
        levels = [g for g in order if g in present]
    elif isinstance(values.dtype, pd.CategoricalDtype):
        levels = [g for g in values.cat.categories if g in present]
    else:
        levels = sorted(present, key=str)
    if reference is not None:
        if reference not in levels:
            raise ValueError(f"Reference '{reference}' not found. Groups: {', '.join(map(str, levels))}")
        if order is None:
            levels = [reference] + [g for g in levels if g != reference]
    return levels


def _below_lod_shares(below_lod: pd.DataFrame, groups: pd.Series, proteins, levels) -> pd.DataFrame:
    b = below_lod.reindex(index=groups.index, columns=proteins).astype(float)
    shares = b.groupby(groups).mean().T
    return shares[levels].rename(columns=lambda g: f"below_lod_{g}")


def compare_groups(table: pd.DataFrame, group: str, *, test: str = "welch", reference=None,
                   ordered: bool = False, order=None, proteins=None,
                   below_lod: pd.DataFrame | None = None, alpha: float = 0.05,
                   level: float = 0.95, verbose: bool = True) -> StatsResult:
    """Compare NPX between groups, protein by protein, with an FDR over all tested proteins.

    Parameters
    ----------
    table : the analysis table from ``to_matrix`` (one row per sample, one column per protein,
        your sample metadata included), or any table with ``proteins=`` given.
    group : the column that holds the groups, e.g. ``"diabetes"``. Samples without a group are
        left out and counted in the report.
    test : ``"welch"`` (default) or ``"mann-whitney"``; see the module text. The same test is
        used for every protein.
    reference : the group the others are compared with; a difference is "other minus reference".
        Default: the first group in ``order``, else the first in alphabetical order.
    ordered : the groups have a natural order: test for a trend over them (needs ``order`` or an
        ordered categorical column).
    order : the groups in their order, e.g. ``["normal", "prediabetes", "diabetes"]``.
    proteins : the protein columns; default ``table.attrs["proteins"]`` from ``to_matrix``.
    below_lod : the table from ``to_matrix(npx, values="below_lod")``; adds the share of values
        below the detection limit per group (a protein can differ by being detectable in one
        group only).
    alpha : the FDR level used in the report (default 0.05).
    level : confidence level of the intervals (default 0.95).
    verbose : print the report.

    Returns a ``StatsResult`` with the summary ``lines`` and the tables ``table`` (one row per
    protein, sorted by p), ``pairwise`` (three or more groups: one row per protein and pair) and
    ``not_tested`` (proteins left out, with the reason).

    Two groups, per protein: ``n_<group>``, ``mean_<group>`` (``median_<group>`` for
    Mann-Whitney), ``difference`` with ``ci_low`` / ``ci_high``, ``p`` and ``fdr``. Three or more
    groups: ``n_`` and ``mean_`` per group, ``statistic``, ``p`` and ``fdr`` of the omnibus or
    trend test; for a Welch trend also the ``slope`` (NPX per group step) with its interval.
    """
    if test not in TESTS:
        raise ValueError(f"test must be one of: {', '.join(TESTS)}")
    if group not in table.columns:
        raise KeyError(f"Column '{group}' not found in the table.")
    if ordered and order is None and not (isinstance(table[group].dtype, pd.CategoricalDtype)
                                          and table[group].cat.ordered):
        raise ValueError("ordered=True needs the order of the groups: pass order=[...].")
    prots = _proteins(table, proteins, exclude=[group])

    lines = []
    has_group = table[group].notna()
    groups = table.loc[has_group, group]
    levels = _levels(groups, reference, order)
    ref = levels[0] if reference is None else reference
    k = len(levels)
    if k == 2:
        levels = [ref] + [g for g in levels if g != ref]
    if k < 2:
        raise ValueError(f"Column '{group}' has fewer than two groups: there is nothing to compare.")
    sizes = groups.value_counts()
    too_small = [g for g in levels if sizes[g] < 2]
    if too_small:
        raise ValueError(f"Groups with fewer than 2 samples: {', '.join(map(str, too_small))}. "
                         f"Leave them out or merge them with another group.")

    mann = test == "mann-whitney"
    if k == 2:
        name = "Mann-Whitney U test" if mann else "Welch's t-test"
        what = ("Hodges-Lehmann shift in NPX" if mann else "difference in mean NPX (= log2 fold change)")
        lines.append(f"test: {name}; effect: {what}, {levels[1]} minus {ref}, "
                     f"with {level:.0%} interval")
    elif ordered:
        lines.append("test: " + ("Jonckheere-Terpstra trend test" if mann else
                                 "linear trend over the group steps (robust HC3 standard errors)")
                     + f"; order: {' < '.join(map(str, levels))}")
    else:
        lines.append("test: " + ("Kruskal-Wallis" if mann else "Welch's ANOVA")
                     + "; then every pair of groups, Holm-adjusted within each protein")
    lines.append("groups: " + ", ".join(f"{g} {sizes[g]}" for g in levels)
                 + (f" (reference: {ref})" if k == 2 or reference is not None else ""))
    if (~has_group).any():
        lines.append(f"left out: {_n(int((~has_group).sum()), 'sample')} without a value in '{group}'")

    center = "median" if mann else "mean"
    rows, pairs, skipped = [], [], []
    data = table.loc[has_group]
    pair_levels = list(combinations(levels, 2))
    for p in prots:
        y = pd.to_numeric(data[p], errors="coerce")
        vals = [y[groups == g].dropna().to_numpy(dtype=float) for g in levels]
        row = {"protein": p}
        for g, v in zip(levels, vals):
            row[f"n_{g}"] = len(v)
        for g, v in zip(levels, vals):
            row[f"{center}_{g}"] = getattr(np, center)(v) if len(v) else np.nan
        few = [str(g) for g, v in zip(levels, vals) if len(v) < 2]
        if few:
            skipped.append({"protein": p, "reason": f"fewer than 2 values in {', '.join(few)}"})
            continue
        if k == 2:
            res = (_mann_whitney_two if mann else _welch_two)(vals[0], vals[1], level)
            if res is None:
                skipped.append({"protein": p, "reason": "no variation in either group"})
                continue
            row.update(difference=res[0], ci_low=res[1], ci_high=res[2], p=res[3])
        else:
            if ordered and not mann:
                res = _linear_trend(vals, level)
                if res is not None:
                    row.update(slope=res[0], ci_low=res[1], ci_high=res[2], p=res[3])
            else:
                res = (_jonckheere if ordered else _kruskal if mann else _welch_anova)(vals)
                if res is not None:
                    row.update(statistic=res[0], p=res[1])
            if res is None:
                reason = "a group without variation" if not mann and not ordered else "no variation"
                skipped.append({"protein": p, "reason": reason})
                continue
            idx = {g: i for i, g in enumerate(levels)}
            for a, b in pair_levels:
                r = (_mann_whitney_two if mann else _welch_two)(vals[idx[a]], vals[idx[b]], level)
                pairs.append({"protein": p, "group": b, "reference": a,
                              "difference": r[0] if r else np.nan,
                              "ci_low": r[1] if r else np.nan, "ci_high": r[2] if r else np.nan,
                              "p": r[3] if r else np.nan})
        rows.append(row)

    out = pd.DataFrame(rows)
    if out.empty:
        raise ValueError("No protein could be tested; see the reasons with verbose output.")
    out["fdr"] = fdr_bh(out["p"])
    if below_lod is not None:
        shares = _below_lod_shares(below_lod, groups, list(out["protein"]), levels)
        out = out.merge(shares, left_on="protein", right_index=True, how="left")
    out = out.sort_values("p", kind="stable").reset_index(drop=True)
    result = StatsResult(title="Group comparison", lines=lines, table=out,
                         not_tested=pd.DataFrame(skipped, columns=["protein", "reason"]))

    if pairs:
        pw = pd.DataFrame(pairs)
        pw["p_holm"] = pw.groupby("protein")["p"].transform(lambda s: pd.Series(holm(s), index=s.index))
        hit = out.set_index("protein")["fdr"] < alpha
        pw["omnibus_significant"] = pw["protein"].map(hit).fillna(False).astype(bool)
        result.pairwise = pw

    lines.append(f"tested: {_n(len(out), 'protein')}"
                 + (f"; not tested: {len(skipped)} (see .not_tested)" if skipped else ""))
    sig = out[out["fdr"] < alpha]
    if sig.empty:
        lines.append(f"FDR below {alpha}: none (lowest FDR {_fmt_p(out['fdr'].min())})")
    else:
        effect = "difference" if k == 2 else "slope" if "slope" in out else None
        shown = []
        for _, r in sig.head(8).iterrows():
            extra = f"{r[effect]:+.2f} NPX, " if effect else ""
            shown.append(f"{r['protein']} ({extra}FDR {_fmt_p(r['fdr'])})")
        more = f", ... ({len(sig)} in total)" if len(sig) > 8 else ""
        lines.append(f"FDR below {alpha}: {_n(len(sig), 'protein')}: " + ", ".join(shown) + more)
    if pairs:
        lines.append(f"pairwise: read .pairwise (p_holm) for the proteins with FDR below {alpha}")
    if verbose:
        print(result)
    return result
