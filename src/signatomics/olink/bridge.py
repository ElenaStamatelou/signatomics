"""Bridging: line up plates (or projects) measured separately, using samples measured on both.

The idea: a *bridging sample* is the same tube measured on two plates. Biology is identical, so
any difference between its two values is technical. Per protein, the **median** difference over
all bridging samples is the plate shift; it is added to every sample of the other plate, so that
plate ends up on the scale of the reference plate. The reference plate is never changed.

Two steps, like the quality steps:

- ``bridge_check``   look: how big is the plate difference, per protein? (changes nothing). Its
                     result holds the tables behind the report, including ``factors``: the
                     shift per protein and plate, which you can inspect or save.
- ``bridge_plates``  apply it; the original values are kept next to the bridged ones

Afterwards, run ``bridge_check`` again: the bridging-sample differences should be about 0.

NPX is on a log2 scale, so a *difference* of 1 NPX means a factor 2. Bridging is therefore an
addition on NPX - never a multiplication. Only NPX data can be bridged this way.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pandas as pd

from .qc import _control_masks, _ids, _is_warning, _protein_col
from .reader import _tidy_sample_ids

ADJUSTED_COLUMNS = ("NPX", "LOD", "MaxLOD")   # shifted together, so "below LOD" stays the same


# ---------------------------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------------------------

def _nz(x: float) -> float:
    """Round to 2 decimals without a "-0.00"."""
    return round(float(x), 2) + 0.0


def _sample_key(df: pd.DataFrame) -> pd.Series:
    """Tidy sample ID: the ``sample`` column from ``add_sample_info``, else made from SampleID."""
    return df["sample"] if "sample" in df.columns else _tidy_sample_ids(df["SampleID"])


def _check_input(npx: pd.DataFrame, batch: str) -> None:
    if "NPX" not in npx.columns:
        raise ValueError("Bridging needs NPX values (log2 scale). This table has no NPX column.")
    if batch not in npx.columns:
        raise KeyError(f"Column '{batch}' not found. Pass the column that tells the plates or "
                       f"projects apart, e.g. batch='PlateID' or batch='Project'.")
    if npx[batch].dropna().nunique() < 2:
        raise ValueError(f"Only one {batch} in the data: there is nothing to bridge.")


def _choose_reference(df: pd.DataFrame, batch: str, reference) -> str:
    plates = df[batch].dropna().unique()
    if reference is not None:
        if reference not in set(plates):
            raise ValueError(f"Reference '{reference}' not found in column '{batch}'. "
                             f"Options: {', '.join(map(str, plates))}")
        return reference
    # default: the plate with the most samples (the most stable scale); ties -> first by name
    counts = df.drop_duplicates(["_sample", batch])[batch].value_counts()
    return sorted(counts[counts == counts.max()].index, key=str)[0]


def _bridge_values(npx: pd.DataFrame, batch: str, bridges, exclude_warnings: bool,
                   exclude_below_lod: bool) -> pd.DataFrame:
    """One value per sample, protein and plate (repeats on a plate averaged), study samples only."""
    ctrl_samples, ctrl_assays = _control_masks(npx)
    df = npx[~(ctrl_samples | ctrl_assays)].copy()
    df["_sample"] = _sample_key(df)
    df = df[df["NPX"].notna() & df[batch].notna() & df["_sample"].notna()]
    if exclude_warnings and "QC_Warning" in df.columns:
        df = df[~_is_warning(df["QC_Warning"])]
    if exclude_below_lod and "LOD" in df.columns:
        df = df[df["NPX"] >= df["LOD"]]
    if bridges is not None:
        df = df[df["_sample"].isin(set(_tidy_sample_ids(pd.Series(list(bridges)).astype(str))))]
    protein = _protein_col(df)
    label = "Assay" if "Assay" in df.columns else protein
    keys = list(dict.fromkeys(["_sample", protein, label, batch]))
    return df.groupby(keys, dropna=False)["NPX"].mean().reset_index()


def _bridge_differences(npx: pd.DataFrame, *, batch: str = "PlateID", reference=None,
                        bridges=None, exclude_warnings: bool = True,
                        exclude_below_lod: bool = False) -> pd.DataFrame:
    """Per bridging sample and protein: how much higher the other plate reads than the reference.

    Returns a long table with one row per bridging sample, protein and non-reference plate:
    ``sample``, the protein columns, ``plate``, ``reference``, ``NPX_reference``, ``NPX_plate``
    and ``difference`` (= plate - reference, in NPX). Repeats on the same plate are averaged.
    The arguments are the same as for ``bridge_check``; shown as ``bridge_check(...).differences``.
    """
    _check_input(npx, batch)
    values = _bridge_values(npx, batch, bridges, exclude_warnings, exclude_below_lod)
    ctrl_samples, ctrl_assays = _control_masks(npx)
    work = npx[~(ctrl_samples | ctrl_assays)]
    ref = _choose_reference(work.assign(_sample=_sample_key(work)), batch, reference)
    protein = _protein_col(values)
    label = "Assay" if "Assay" in values.columns else protein
    keys = ["_sample", protein] + ([label] if label != protein else [])
    on_ref = values[values[batch] == ref][keys + ["NPX"]].rename(columns={"NPX": "NPX_reference"})
    others = values[values[batch] != ref].rename(columns={batch: "plate", "NPX": "NPX_plate"})
    out = others.merge(on_ref, on=keys, how="inner")
    out["reference"] = ref
    out["difference"] = out["NPX_plate"] - out["NPX_reference"]
    out = out.rename(columns={"_sample": "sample"})
    cols = ["sample", protein] + ([label] if label != protein else []) + [
        "plate", "reference", "NPX_reference", "NPX_plate", "difference"]
    out = out[cols].reset_index(drop=True)
    out.attrs["reference"] = ref
    return out


# ---------------------------------------------------------------------------------------------
# Look: bridge_check
# ---------------------------------------------------------------------------------------------

@dataclass
class BridgeCheck:
    """Result of ``bridge_check``: summary lines plus the tables behind them."""

    lines: list[str] = field(default_factory=list)
    differences: pd.DataFrame | None = None  # one row per bridging sample, protein and plate
    proteins: pd.DataFrame | None = None     # per protein and plate: shift, spread, bridges
    samples: pd.DataFrame | None = None      # per bridging sample and plate: overall shift
    factors: pd.DataFrame | None = None      # per plate and protein: the adjustment bridge_plates adds
    by_group: pd.DataFrame | None = None     # per protein: shift in each group (if group given)

    def __str__(self) -> str:
        return "Bridging check\n  " + "\n  ".join(self.lines)

    def _repr_markdown_(self) -> str:  # nice display in notebooks
        return "**Bridging check**\n\n" + "\n".join(f"- {line}" for line in self.lines)


def bridge_check(npx: pd.DataFrame, *, batch: str = "PlateID", reference=None, group: str | None = None,
                 samples=None, on: str | None = None,
                 bridges=None, exclude_warnings: bool = True, exclude_below_lod: bool = False,
                 min_bridges: int = 8, tolerance: float = 0.1, show_ids: bool = True,
                 verbose: bool = True) -> BridgeCheck:
    """Look at the difference between plates, using the bridging samples (changes nothing).

    Run it before bridging (how big is the plate difference?) and after (is it gone?).

    Parameters
    ----------
    npx : table from ``read_npx`` (several plates in one table; ``pd.concat`` separate files).
    batch : the column that tells the plates or projects apart (default ``"PlateID"``; use e.g.
        ``"Project"`` when bridging two whole projects whose plate names may overlap).
    reference : the plate the others are compared with; default the plate with most samples.
    samples, on : your sample metadata and its ID column, when ``group`` is not yet a column of
        ``npx`` (only looked up, as in ``qc_report``).
    group : a column of your sample metadata (e.g. ``"diabetes"``). Checks whether the bridging
        samples represent the groups, and whether the plate shift is the same in every group -
        if not, bridging could add or remove a group difference.
    bridges : the sample IDs to use as bridging samples; default every sample measured on the
        reference and on another plate.
    exclude_warnings : do not use bridging measurements with a quality warning (default True).
    exclude_below_lod : use only pairs where both values are above the detection limit
        (default False: below-LOD values are noisy but not biased, and dropping them leaves few
        pairs for low proteins).
    min_bridges : proteins with fewer usable pairs are flagged (``few_bridges`` in ``factors``).
    tolerance : a protein whose median shift is larger than this (in NPX) is counted as
        "still shifted".
    show_ids : list sample IDs in the report; ``False`` gives counts only (for sharing).

    Returns a ``BridgeCheck``: the summary lines, plus the tables ``differences`` (per bridging
    sample, protein and plate), ``proteins``, ``samples``, ``factors`` (per plate and protein:
    the ``adjustment`` that ``bridge_plates`` adds, ``bridges`` used, ``spread``,
    ``few_bridges``) and, with ``group``, ``by_group``.
    """
    from .qc import _with_groups
    npx = _with_groups(npx, samples, on)
    d = _bridge_differences(npx, batch=batch, reference=reference, bridges=bridges,
                            exclude_warnings=exclude_warnings, exclude_below_lod=exclude_below_lod)
    ref = d.attrs["reference"]
    protein = _protein_col(d)
    label = "Assay" if "Assay" in d.columns else protein
    chk = BridgeCheck(differences=d)
    chk.lines.append(f"reference: {ref}")
    if d.empty:
        chk.lines.append("WARNING no bridging samples found: no sample is measured on the reference "
                         "and on another plate (check that the sample IDs are written the same way)")
        if verbose:
            print(chk)
        return chk

    p = (d.groupby(["plate", protein] + ([label] if label != protein else []))["difference"]
         .agg(shift="median", spread=lambda s: s.quantile(0.75) - s.quantile(0.25), bridges="size")
         .reset_index())
    chk.proteins = p
    s = (d.groupby(["plate", "sample"])["difference"].median().rename("shift").reset_index())
    chk.samples = s
    chk.factors = _bridge_factors(npx, d, batch=batch, min_bridges=min_bridges,
                                  exclude_warnings=exclude_warnings,
                                  exclude_below_lod=exclude_below_lod)

    for plate, pp in p.groupby("plate", sort=False):
        n_br = d.loc[d["plate"] == plate, "sample"].nunique()
        shifted = pp[pp["shift"].abs() > tolerance]
        chk.lines.append(
            f"{plate}: reads {_nz(pp['shift'].median()):+.2f} NPX compared with {ref} (median over "
            f"assays; per assay {_nz(pp['shift'].min()):+.2f} to {_nz(pp['shift'].max()):+.2f}), "
            f"measured on {n_br} bridge samples")
        chk.lines.append(
            f"  assays shifted by more than {tolerance:g} NPX: {len(shifted)} of {len(pp)}"
            + (" - bridging needed" if len(shifted) else " - the plates are lined up"))
        chk.lines.append(
            f"  noise: the bridge samples' differences spread {pp['spread'].median():.2f} NPX per "
            f"assay (typical IQR); bridging removes the shift; this noise stays")
        if n_br < min_bridges:
            chk.lines.append(f"  WARNING only {n_br} bridge samples: the shift per assay is "
                             f"uncertain ({min_bridges} or more are usually advised)")
        ss = s[s["plate"] == plate]
        med = ss["shift"].median()
        mad = (ss["shift"] - med).abs().median() * 1.4826
        odd = ss[(ss["shift"] - med).abs() > max(3 * mad, 0.25)] if len(ss) >= 5 else ss.iloc[0:0]
        if len(odd):
            chk.lines.append(f"  bridge samples that shift differently from the others "
                             f"(a mix-up or a degraded tube?): {_ids(odd['sample'], show_ids)}")

    if "Bridge_adjustment" in npx.columns:
        chk.lines.append(
            "note: this table is already bridged with these bridge samples, so shifts of about 0 "
            "are expected by construction. Independent checks: plot_bridge_pca (the same samples "
            "should now overlap) and, with group=..., whether the shift was the same in every group")

    if group is not None:
        if group not in npx.columns:
            raise KeyError(f"Group column '{group}' not found (did you run add_sample_info?).")
        groups = (npx.assign(_sample=_sample_key(npx)).dropna(subset=[group])
                  .drop_duplicates("_sample").set_index("_sample")[group])
        everyone = groups.value_counts()
        in_bridges = groups.reindex(d["sample"].unique()).value_counts()
        parts = [f"{g}: {int(in_bridges.get(g, 0))} of {int(n)}" for g, n in everyone.items()]
        chk.lines.append(f"'{group}' among the bridge samples: " + "; ".join(parts))
        missing = [str(g) for g in everyone.index if in_bridges.get(g, 0) == 0]
        if missing:
            chk.lines.append(f"  WARNING no bridge samples from: {', '.join(missing)} - it cannot "
                             f"be checked whether the plate shift is the same for this group")
        dg = d.assign(_group=d["sample"].map(groups)).dropna(subset=["_group"])
        per_sample = dg.groupby(["plate", "_group", "sample"])["difference"].median()
        overall = per_sample.groupby(["plate", "_group"]).median()
        for plate in overall.index.get_level_values(0).unique():
            o = overall.loc[plate]
            chk.lines.append(f"  overall shift on {plate} per group: "
                             + "; ".join(f"{g} {_nz(v):+.2f}" for g, v in o.items())
                             + f" (largest difference {o.max() - o.min():.2f} NPX)")
        by = (dg.groupby([label, "_group"])["difference"].median().unstack("_group"))
        by["difference"] = by.max(axis=1) - by.min(axis=1)
        chk.by_group = by.sort_values("difference", ascending=False)
        big = chk.by_group[chk.by_group["difference"] > 0.5]
        if len(big):
            chk.lines.append(
                f"  assays whose shift differs by more than 0.5 NPX between groups: "
                f"{', '.join(map(str, big.index[:10]))}" + (f" ... ({len(big)} in total)" if len(big) > 10 else "")
                + " - with few bridge samples per group this is often noise; look at them before "
                  "trusting group differences for these assays")

    if verbose:
        print(chk)
    return chk


# ---------------------------------------------------------------------------------------------
# Compute and apply: _bridge_factors, bridge_plates
# ---------------------------------------------------------------------------------------------

def _bridge_factors(npx: pd.DataFrame, d: pd.DataFrame, *, batch: str, min_bridges: int,
                    exclude_warnings: bool, exclude_below_lod: bool) -> pd.DataFrame:
    """The shift per protein and plate, from the bridging-sample differences ``d``.

    For every protein and every plate other than the reference:
    ``adjustment = -median(NPX on the plate - NPX on the reference)`` over the bridging samples.
    Adding it to the plate's NPX puts that plate on the reference scale. The median is used (not
    the mean) so that one odd bridging sample cannot pull the shift.

    Returns one row per plate and protein: ``plate``, ``reference``, protein columns,
    ``adjustment`` (NPX to add), ``bridges`` (pairs used) and ``spread`` (IQR of the
    differences), ``few_bridges``. Proteins without any pair get an empty adjustment.
    Shown as ``bridge_check(...).factors``.
    """
    ref = d.attrs["reference"]
    ctrl_samples, ctrl_assays = _control_masks(npx)
    study = npx[~(ctrl_samples | ctrl_assays)]
    protein = _protein_col(study)
    label = "Assay" if "Assay" in study.columns else protein
    names = ([protein, label] if label != protein else [protein])
    grid = (study[[batch] + names].dropna(subset=[batch]).drop_duplicates([batch, protein])
            .rename(columns={batch: "plate"}))
    grid = grid[grid["plate"] != ref]
    stats = (d.groupby(["plate", protein])["difference"]
             .agg(shift="median", bridges="size",
                  spread=lambda s: s.quantile(0.75) - s.quantile(0.25))
             .reset_index())
    f = grid.merge(stats, on=["plate", protein], how="left")
    f["bridges"] = f["bridges"].fillna(0).astype(int)
    f["adjustment"] = -f.pop("shift")
    f["reference"] = ref
    f["few_bridges"] = f["bridges"] < min_bridges
    cols = ["plate", "reference"] + names + ["adjustment", "bridges", "spread", "few_bridges"]
    f = f[cols].sort_values(["plate", label]).reset_index(drop=True)
    f.attrs["settings"] = {"batch": batch, "exclude_warnings": exclude_warnings,
                           "exclude_below_lod": exclude_below_lod}
    return f


def bridge_plates(npx: pd.DataFrame, *, factors: pd.DataFrame | None = None, batch: str = "PlateID",
                  reference=None, bridges=None, exclude_warnings: bool = True,
                  exclude_below_lod: bool = False, min_bridges: int = 8,
                  verbose: bool = True) -> pd.DataFrame:
    """Put every plate on the scale of the reference plate, using the bridging samples.

    The adjustment per protein and plate (the median shift of the bridging samples, as in
    ``bridge_check(...).factors``) is added to NPX, and to LOD and MaxLOD as well - so a value
    that was below the detection limit stays below it. The reference plate is not changed. The
    original values are kept in ``NPX_unbridged`` (and ``LOD_unbridged``), the adjustment used
    in ``Bridge_adjustment``.

    Parameters
    ----------
    npx : table from ``read_npx`` with all plates in one table (best after ``apply_qc``).
    factors : a factors table to apply, e.g. ``bridge_check(...).factors`` checked or saved
        earlier; if not given, it is computed here with the arguments below.
    batch, reference, bridges, exclude_warnings, exclude_below_lod, min_bridges : as in
        ``bridge_check``.
    verbose : print the bridging report.

    Returns the bridged table; the report is in ``df.attrs["bridge_report"]``.
    """
    if "Bridge_adjustment" in npx.columns:
        raise ValueError("This table is already bridged (it has a Bridge_adjustment column). "
                         "Bridge the original table instead, so values are not shifted twice.")
    _check_input(npx, batch)
    if factors is None:
        d = _bridge_differences(npx, batch=batch, reference=reference, bridges=bridges,
                                exclude_warnings=exclude_warnings,
                                exclude_below_lod=exclude_below_lod)
        factors = _bridge_factors(npx, d, batch=batch, min_bridges=min_bridges,
                                  exclude_warnings=exclude_warnings,
                                  exclude_below_lod=exclude_below_lod)
    protein = _protein_col(npx)
    if protein not in factors.columns:
        raise ValueError(f"The factors table has no '{protein}' column; make it with bridge_check "
                         f"on data with the same protein columns.")
    ref = factors["reference"].iloc[0] if len(factors) else reference
    label = "Assay" if "Assay" in factors.columns else protein

    df = npx.copy()
    adj = df[[batch, protein]].merge(
        factors[["plate", protein, "adjustment"]].rename(columns={"plate": batch}),
        on=[batch, protein], how="left")["adjustment"].to_numpy()
    on_ref = (df[batch] == ref).to_numpy()
    adj = np.where(on_ref, 0.0, adj)
    df["NPX_unbridged"] = df["NPX"]
    if "LOD" in df.columns:
        df["LOD_unbridged"] = df["LOD"]
    for col in ADJUSTED_COLUMNS:
        if col in df.columns:
            df[col] = df[col] + np.nan_to_num(adj, nan=0.0)
    df["Bridge_adjustment"] = adj

    report = [f"reference (kept as it is): {ref}"]
    for plate, fp in factors.groupby("plate", sort=False):
        ok = fp[fp["adjustment"].notna()]
        if len(ok):
            report.append(
                f"{plate}: {len(ok)} assays adjusted by their adjustment factor (median "
                f"{ok['adjustment'].median():+.2f} NPX, from {ok['adjustment'].min():+.2f} to "
                f"{ok['adjustment'].max():+.2f}), each from {int(ok['bridges'].min())}"
                + (f"-{int(ok['bridges'].max())}" if ok['bridges'].max() != ok['bridges'].min() else "")
                + " bridge samples")
        none = fp[fp["adjustment"].isna()]
        if len(none):
            report.append(f"  WARNING not adjusted (no usable bridge sample): "
                          f"{', '.join(map(str, none[label]))}")
        few = fp[fp["few_bridges"] & fp["adjustment"].notna()]
        if len(few):
            report.append(f"  few bridge samples (< {min_bridges}): {len(few)} assays, e.g. "
                          f"{', '.join(map(str, few[label].head(5)))}")
    unknown = sorted(set(df.loc[~on_ref & np.isnan(adj), batch].dropna().astype(str))
                     - set(factors["plate"].astype(str)))
    if unknown:
        report.append(f"WARNING plates without adjustment factors (left as they are): {', '.join(unknown)}")
    report.append("NPX and LOD shifted together, so 'below the LOD' stays the same; the original "
                  "values are in NPX_unbridged and LOD_unbridged")
    report.append("check: plot_bridge_pca(before, after) - each bridge sample's two measurements "
                  "should now lie together")

    df.attrs = dict(npx.attrs)
    df.attrs["bridge_report"] = report
    if verbose:
        print("Bridging report\n  " + "\n  ".join(report))
    return df
