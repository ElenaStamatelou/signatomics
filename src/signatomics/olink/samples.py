"""Attach sample information (clinical data, groups) to Olink data — and report what did not match.

Silently losing samples while matching IDs is one of the most common hidden errors in omics
analyses, so the matching always reports: what matched, what matched only after tidying the IDs,
and what did not match on either side.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from .reader import _tidy_sample_ids


def _read_sample_list(samples) -> pd.DataFrame:
    if isinstance(samples, pd.DataFrame):
        return samples.copy()
    path = Path(samples)
    if path.suffix.lower() in (".xlsx", ".xls"):
        return pd.read_excel(path, dtype={0: str})
    return pd.read_csv(path, sep=None, engine="python", encoding="utf-8-sig")


def _examples(ids, n=10) -> str:
    ids = sorted(map(str, ids))
    return ", ".join(ids[:n]) + (f", ... ({len(ids)} in total)" if len(ids) > n else "")


def add_sample_info(npx: pd.DataFrame, samples, *, on: str, npx_on: str = "SampleID",
                    verbose: bool = True) -> pd.DataFrame:
    """Add the columns of a sample list (e.g. diabetes status, age, sex) to Olink data.

    Parameters
    ----------
    npx : table from ``read_npx``.
    samples : the sample list — a table, or a path to a csv / Excel file.
    on : the column in the sample list that holds the same IDs as the protein file
         (e.g. ``"FIF ID"``).
    npx_on : the ID column in the protein data (default ``"SampleID"``).
    verbose : print the matching report.

    IDs are tidied on both sides before matching (spaces at the ends are removed; spaces,
    underscores, hyphens and capitals do not matter, so ``"308 BIS"`` matches ``"308_bis"``;
    leading zeros are kept), the original IDs are kept, and the tidy ID is added as
    column ``sample``. Sample-list columns whose name already exists in the protein data get the
    ending ``_samplelist``. Protein rows without a match are kept, with empty sample information.

    Returns the protein data with the sample columns added; the matching report is in
    ``df.attrs["sample_report"]``.
    """
    info = _read_sample_list(samples)
    info.columns = [str(c).strip() for c in info.columns]
    if on not in info.columns:
        raise KeyError(f"Column '{on}' not found in the sample list. Columns: {', '.join(info.columns)}")
    if npx_on not in npx.columns:
        raise KeyError(f"Column '{npx_on}' not found in the protein data.")

    report = []
    info = info[info[on].notna()].copy()
    info["sample"] = _tidy_sample_ids(info[on].astype(str).str.replace(r"\.0$", "", regex=True))
    duplicated = info.loc[info["sample"].duplicated(keep=False), "sample"].unique()
    if len(duplicated):
        raise ValueError(
            f"The sample list has IDs that occur more than once (after tidying): "
            f"{_examples(duplicated)}. Each sample may appear only once — please fix the list."
        )

    df = npx.copy()
    df["sample"] = _tidy_sample_ids(df[npx_on])

    protein_ids, list_ids = set(df["sample"]), set(info["sample"])
    raw_protein = set(df[npx_on].astype(str).str.strip())
    raw_list = set(info[on].astype(str).str.strip().str.replace(r"\.0$", "", regex=True))
    matched = protein_ids & list_ids
    exact = raw_protein & raw_list

    report.append(f"matched on: {npx_on} (NPX data) = {on} (your sample metadata)")
    report.append(f"matched: {len(matched)} of {len(protein_ids)} sample IDs in the NPX data"
                  + (f" ({len(matched) - len(exact)} only after tidying the IDs, e.g. \"308 BIS\" = \"308_BIS\")"
                     if len(matched) > len(exact) else ""))
    if protein_ids - list_ids:
        report.append(f"WARNING missing from your sample metadata "
                      f"(kept, without sample metadata): {_examples(protein_ids - list_ids)}")
    if list_ids - protein_ids:
        report.append(f"missing from the NPX data: {_examples(list_ids - protein_ids)}")

    extra = [c for c in info.columns if c not in ("sample",)]
    clashes = [c for c in extra if c in df.columns]
    info = info.rename(columns={c: f"{c}_samplelist" for c in clashes})
    if clashes:
        report.append("columns in both tables (yours renamed with the ending _samplelist): "
                      + ", ".join(clashes))

    out = df.merge(info, on="sample", how="left", validate="many_to_one")
    out.attrs = dict(npx.attrs)
    out.attrs["sample_report"] = report
    if verbose:
        print("Matching report\n  " + "\n  ".join(report))
    return out
