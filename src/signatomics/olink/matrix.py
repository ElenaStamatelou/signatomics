"""The analysis table: one row per sample, one column per assay, your sample metadata in front.

This is the hand-off between preprocessing and analysis: ``signatomics.stats`` (and later the
signature search) work on this table. Run it last, after quality control, bridging and adding
your sample metadata.

What ``to_matrix`` does, and reports:

- control samples and control assays are left out;
- a sample measured more than once (a bridge sample after bridging, or a repeat) gets the
  **average** of its NPX values, so every sample counts once;
- each assay becomes one column, named after the assay; an assay name used by more than one
  OlinkID (e.g. the same assay on several Explore panels) gets the OlinkID added;
- the columns of your sample metadata (everything ``add_sample_info`` added) come first.

The assay columns are listed in ``table.attrs["proteins"]``, so the analysis functions know which
columns hold NPX values and which hold sample metadata.
"""

from __future__ import annotations

import pandas as pd

from .qc import _control_masks, _protein_col
from .reader import COLUMN_VARIANTS, NUMERIC, _tidy_sample_ids

# Columns that describe a measurement or an assay, never a sample's metadata
OLINK_COLUMNS = set(COLUMN_VARIANTS) | set(NUMERIC) | {
    "sample", "NPX_unbridged", "LOD_unbridged", "Bridge_adjustment", "PlateLOD", "Unit",
    "Sample_Type", "ExploreVersion", "Block", "Count", "ExtNPX", "SampleQC", "AssayQC",
    "DataAnalysisRefID", "Project"}


def _assay_names(df: pd.DataFrame, protein: str) -> pd.Series:
    """Column name per row: the assay name, with the OlinkID added where a name is used twice."""
    if protein == "Assay" or "Assay" not in df.columns:
        return df[protein].astype(str)
    pairs = df[["Assay", protein]].drop_duplicates()
    shared = set(pairs.loc[pairs["Assay"].duplicated(keep=False), "Assay"])
    name = df["Assay"].astype(str)
    clash = df["Assay"].isin(shared)
    name[clash] = name[clash] + " (" + df.loc[clash, protein].astype(str) + ")"
    return name


def _sample_columns(df: pd.DataFrame) -> list[str]:
    """Your sample metadata: columns that are not Olink's and hold one value per sample."""
    candidates = [c for c in df.columns if c not in OLINK_COLUMNS and not c.startswith("_")]
    if not candidates:
        return []
    per_sample = df.groupby("_sample")[candidates].nunique(dropna=True)
    return [c for c in candidates if (per_sample[c] <= 1).all()]


def to_matrix(npx: pd.DataFrame, *, values: str = "NPX", verbose: bool = True) -> pd.DataFrame:
    """Turn the long NPX table into the analysis table: one row per sample, one column per assay.

    Parameters
    ----------
    npx : the preprocessed long table (from ``prepare``, or from ``read_npx`` and the steps after
        it). Sample metadata added with ``add_sample_info`` is carried over.
    values : what goes in the assay columns:
        ``"NPX"`` (default) the NPX value, averaged when a sample was measured more than once;
        ``"below_lod"`` True where that NPX value is below the detection limit (LOD averaged the
        same way). Pass this second table to ``compare_groups(below_lod=...)`` to see the share
        below LOD per group.
    verbose : print the report.

    Returns a table with one row per sample (index ``sample``: the tidy sample ID), your sample
    metadata first, then one column per assay. ``table.attrs["proteins"]`` lists the assay
    columns; ``table.attrs["matrix_report"]`` holds the report.
    """
    if values not in ("NPX", "below_lod"):
        raise ValueError("values must be 'NPX' or 'below_lod'.")
    if "NPX" not in npx.columns:
        raise ValueError("The analysis table needs NPX values. This table has no NPX column.")
    if values == "below_lod" and "LOD" not in npx.columns:
        raise ValueError("values='below_lod' needs a LOD column; this table has none.")

    report = []
    ctrl_samples, ctrl_assays = _control_masks(npx)
    df = npx[~(ctrl_samples | ctrl_assays)].copy()
    if ctrl_samples.any() or ctrl_assays.any():
        parts = []
        if ctrl_samples.any():
            parts.append(f"{npx.loc[ctrl_samples, 'SampleID'].nunique()} control samples")
        if ctrl_assays.any():
            parts.append(f"{npx.loc[ctrl_assays, _protein_col(npx)].nunique()} control assays")
        report.append("left out: " + " and ".join(parts))

    df["_sample"] = df["sample"] if "sample" in df.columns else _tidy_sample_ids(df["SampleID"])
    df = df[df["_sample"].notna()]
    protein = _protein_col(df)
    df["_assay"] = _assay_names(df, protein)

    counts = df[df["NPX"].notna()].groupby(["_sample", "_assay"]).size()
    several = counts[counts > 1].index.get_level_values("_sample").unique()
    if len(several):
        report.append(f"averaged: {len(several)} samples measured more than once (bridge samples "
                      f"or repeats) get the average of their NPX values")

    cells = df.groupby(["_sample", "_assay"])[["NPX"] + (["LOD"] if values == "below_lod" else [])].mean()
    if values == "NPX":
        wide = cells["NPX"].unstack("_assay")
    else:
        below = (cells["NPX"] < cells["LOD"]).where(cells["NPX"].notna() & cells["LOD"].notna())
        wide = below.unstack("_assay")
    order = list(pd.unique(df["_assay"]))
    wide = wide[[a for a in order if a in wide.columns]]

    meta_cols = _sample_columns(df)
    if meta_cols:
        meta = df.groupby("_sample")[meta_cols].first()
        clash = [c for c in meta_cols if c in wide.columns]
        if clash:
            raise ValueError(f"Sample metadata columns with the name of an assay: {', '.join(clash)}. "
                             f"Rename them in your sample metadata.")
        table = meta.join(wide, how="right")
    else:
        table = wide
    table.index.name = "sample"
    table.columns.name = None

    missing = int(wide.isna().sum().sum())
    report.append(f"result: {len(table)} samples x {wide.shape[1]} assays"
                  + (f", plus {len(meta_cols)} columns of sample metadata ({', '.join(meta_cols)})"
                     if meta_cols else " (no sample metadata: add it with add_sample_info)"))
    if missing:
        report.append(f"empty cells: {missing} (a sample without a value for that assay)")

    table.attrs["proteins"] = list(wide.columns)
    table.attrs["matrix_report"] = report
    if verbose:
        print("Analysis table\n  " + "\n  ".join(report))
    return table
