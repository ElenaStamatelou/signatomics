"""Helpers shared by the statistics functions: which columns are proteins, FDR, Holm, numbers."""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pandas as pd


def _proteins(table: pd.DataFrame, proteins, exclude: list[str]) -> list[str]:
    """The protein columns: given, or from ``to_matrix`` (``table.attrs["proteins"]``)."""
    if proteins is None:
        proteins = table.attrs.get("proteins")
        if proteins is None:
            raise ValueError("Which columns are proteins? Make the table with to_matrix, or pass "
                             "proteins=[...] (a list of column names).")
    proteins = [p for p in proteins if p not in exclude]
    missing = [p for p in proteins if p not in table.columns]
    if missing:
        raise KeyError(f"Protein columns not in the table: {', '.join(map(str, missing[:10]))}")
    if not proteins:
        raise ValueError("No protein columns to test.")
    return proteins


def fdr_bh(p) -> np.ndarray:
    """Benjamini-Hochberg FDR. Missing p-values stay missing and do not count in the family."""
    p = np.asarray(p, dtype=float)
    out = np.full(p.shape, np.nan)
    ok = ~np.isnan(p)
    m = int(ok.sum())
    if m == 0:
        return out
    pv = p[ok]
    order = np.argsort(pv)
    ranked = pv[order] * m / np.arange(1, m + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    adj = np.empty(m)
    adj[order] = np.minimum(ranked, 1.0)
    out[ok] = adj
    return out


def holm(p) -> np.ndarray:
    """Holm adjustment within one family. Missing p-values stay missing."""
    p = np.asarray(p, dtype=float)
    out = np.full(p.shape, np.nan)
    ok = ~np.isnan(p)
    m = int(ok.sum())
    if m == 0:
        return out
    pv = p[ok]
    order = np.argsort(pv)
    stepped = np.maximum.accumulate(pv[order] * (m - np.arange(m)))
    adj = np.empty(m)
    adj[order] = np.minimum(stepped, 1.0)
    out[ok] = adj
    return out


def _fmt_p(p: float) -> str:
    if pd.isna(p):
        return "-"
    return f"{p:.2g}" if p >= 0.001 else f"{p:.1e}"


def _n(count: int, word: str) -> str:
    return f"{count} {word}{'s' if count != 1 else ''}"


@dataclass
class StatsResult:
    """Result of a statistics function: summary lines plus the tables behind them."""

    title: str = ""
    lines: list[str] = field(default_factory=list)
    table: pd.DataFrame | None = None       # one row per protein
    pairwise: pd.DataFrame | None = None    # 3+ groups: one row per protein and pair of groups
    not_tested: pd.DataFrame | None = None  # proteins left out, with the reason

    def __str__(self) -> str:
        return self.title + "\n  " + "\n  ".join(self.lines)

    def _repr_markdown_(self) -> str:  # nice display in notebooks
        return f"**{self.title}**\n\n" + "\n".join(f"- {line}" for line in self.lines)
