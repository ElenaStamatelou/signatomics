"""Quality checks for Olink data: look first (``qc_report``), then decide (``apply_qc``).

Every choice that removes data is a visible argument with a sensible default, and everything
that is removed is reported. Nothing is removed silently.

Control samples (negative, plate and sample controls) and control assays (incubation,
extension, detection, amplification controls) are reported separately and left out of all other
counts: they are part of the lab process, not of the study.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import pandas as pd

PASS_WORDS = {"PASS", "PASSED", "OK"}


def _value_col(df: pd.DataFrame) -> str:
    for col in ("NPX", "Quantified_value", "Ct"):
        if col in df.columns:
            return col
    raise KeyError("No value column (NPX, Quantified_value or Ct) found.")


def _protein_col(df: pd.DataFrame) -> str:
    return "OlinkID" if "OlinkID" in df.columns else "Assay"


def _measurement_keys(df: pd.DataFrame) -> list[str]:
    """One measurement = one sample in one well of one plate (per panel for multi-panel data)."""
    return [c for c in ("SampleID", "Index", "WellID", "Panel", "PlateID") if c in df.columns]


def _measure_col(df: pd.DataFrame) -> str | None:
    """Column that tells two measurements of the same sample on the same plate apart."""
    return next((c for c in ("Index", "WellID") if c in df.columns), None)


def _is_warning(qc: pd.Series) -> pd.Series:
    return qc.notna() & ~qc.astype(str).str.strip().str.upper().isin(PASS_WORDS)


def _control_masks(df: pd.DataFrame) -> tuple[pd.Series, pd.Series]:
    """Rows of control samples, and rows of control assays."""
    no = pd.Series(False, index=df.index)
    samples = (df["SampleType"].astype("string").str.upper().ne("SAMPLE") & df["SampleType"].notna()
               if "SampleType" in df.columns else no)
    assays = (df["AssayType"].astype("string").str.lower().ne("assay") & df["AssayType"].notna()
              if "AssayType" in df.columns else no)
    return samples.fillna(False).astype(bool), assays.fillna(False).astype(bool)


def _ids(values, show_ids: bool, n: int = 15) -> str:
    values = [str(v) for v in pd.unique(pd.Series(values).dropna())]
    if not show_ids:
        return f"{len(values)} (IDs hidden)"
    more = f", ... ({len(values)} in total)" if len(values) > n else ""
    return ", ".join(values[:n]) + more


def _measurement_table(df: pd.DataFrame) -> pd.DataFrame:
    value, keys = _value_col(df), _measurement_keys(df)
    agg = {"values": (value, "size"), "missing": (value, lambda s: int(s.isna().sum()))}
    if "QC_Warning" in df.columns:
        agg["qc_warning"] = ("QC_Warning", "first")
    m = df.groupby(keys, dropna=False).agg(**agg).reset_index()
    m["empty"] = m["missing"] == m["values"]
    m["warning"] = _is_warning(m["qc_warning"]) if "qc_warning" in m else False
    per_plate = [c for c in ("SampleID", "Panel", "PlateID") if c in m.columns]
    measure = _measure_col(m)
    m["repeat"] = (m.groupby(per_plate, dropna=False)[measure].transform("nunique") > 1
                   if measure else False)
    m["bridging"] = (m.groupby("SampleID")["PlateID"].transform("nunique") > 1
                     if "PlateID" in m else False)
    return m


def below_lod_by_group(npx: pd.DataFrame, group: str) -> pd.DataFrame:
    """Share (0-1) of measurements below the detection limit, per protein and per group.

    Returns one row per protein, one column per group, plus ``difference`` (largest minus
    smallest share). A large difference means the protein is detectable in one group but not in
    another - that can itself be a biological signal, so such proteins should not be removed
    just because they are often below LOD overall.
    """
    ctrl_samples, ctrl_assays = _control_masks(npx)
    df = npx[~(ctrl_samples | ctrl_assays)]
    value, protein = _value_col(df), _protein_col(df)
    if "LOD" not in df.columns:
        raise ValueError("No LOD column.")
    if group not in df.columns:
        raise KeyError(f"Group column '{group}' not found (did you run add_sample_info?).")
    below = (df[value] < df["LOD"]).where(df[value].notna() & df["LOD"].notna())
    label = "Assay" if "Assay" in df.columns else protein
    t = (df.assign(_below=below).dropna(subset=[group])
         .groupby([label, group])["_below"].mean().unstack(group))
    t["difference"] = t.max(axis=1) - t.min(axis=1)
    return t.sort_values("difference", ascending=False)


@dataclass
class QCReport:
    """Result of ``qc_report``: summary lines plus the tables behind them."""

    lines: list[str] = field(default_factory=list)
    measurements: pd.DataFrame | None = None    # one row per sample measurement (no controls)
    plates: pd.DataFrame | None = None          # one row per plate
    proteins: pd.DataFrame | None = None        # one row per protein (no control assays)
    group_by_plate: pd.DataFrame | None = None  # counts of each group on each plate
    below_lod_by_group: pd.DataFrame | None = None  # share below LOD per protein and group

    def __str__(self) -> str:
        return "Quality report\n  " + "\n  ".join(self.lines)

    def _repr_markdown_(self) -> str:  # nice display in notebooks
        return "**Quality report**\n\n" + "\n".join(f"- {line}" for line in self.lines)


def qc_report(npx: pd.DataFrame, *, group: str | None = None, max_below_lod: float = 0.5,
              show_ids: bool = True, verbose: bool = True) -> QCReport:
    """Look at data quality without changing anything.

    Checks: control samples and assays, quality warnings, empty measurements, samples measured
    more than once on the same plate (repeats), samples on more than one plate (bridging
    samples), proteins mostly below the detection limit, and — if ``group`` is given — whether a
    group is mixed up with plates.

    Parameters
    ----------
    npx : table from ``read_npx`` (optionally with ``add_sample_info``).
    group : a sample-information column to check against the plates, e.g. ``"diabetes"``.
    max_below_lod : share (0-1) above which a protein is reported as "mostly below LOD".
    show_ids : list sample IDs in the report; ``False`` gives counts only (for sharing).
    """
    value, protein = _value_col(npx), _protein_col(npx)
    rep = QCReport()

    ctrl_samples, ctrl_assays = _control_masks(npx)
    if ctrl_samples.any() or ctrl_assays.any():
        ctrl = npx[ctrl_samples]
        n_meas = ctrl.groupby(_measurement_keys(ctrl), dropna=False).ngroups if len(ctrl) else 0
        kinds = (ctrl.drop_duplicates(_measurement_keys(ctrl))["SampleType"].value_counts()
                 if len(ctrl) else pd.Series(dtype=int))
        rep.lines.append(
            f"controls (left out of everything below; apply_qc removes them): {n_meas} control "
            f"measurements" + (f" ({', '.join(f'{n} {k}' for k, n in kinds.items())})" if len(kinds) else "")
            + f", {npx.loc[ctrl_assays, protein].nunique()} control assays")
    df = npx[~(ctrl_samples | ctrl_assays)]

    no_id = df["SampleID"].isna()
    if no_id.any():
        rep.lines.append(f"WARNING {int(no_id.sum())} rows without a sample ID")

    m = _measurement_table(df)
    rep.measurements = m
    rep.lines.append(f"{len(m)} sample measurements of {df['SampleID'].nunique()} sample IDs, "
                     f"{df[protein].nunique()} proteins"
                     + (f", {df['PlateID'].nunique()} plates" if "PlateID" in df else ""))
    rep.lines.append(f"quality warnings: {int(m['warning'].sum())} measurements")
    if m["empty"].any():
        rep.lines.append(f"WARNING completely empty measurements (no values at all): "
                         f"{_ids(m.loc[m['empty'], 'SampleID'], show_ids)}")
    if m["repeat"].any():
        rep.lines.append(f"measured more than once on the same plate (repeats): "
                         f"{_ids(m.loc[m['repeat'], 'SampleID'], show_ids)}")
    if m["bridging"].any():
        rep.lines.append(f"on more than one plate (bridging samples): {m.loc[m['bridging'], 'SampleID'].nunique()}")
    if _measure_col(df) is None:
        rep.lines.append("no Index or WellID column: repeats on the same plate cannot be recognised")

    if "PlateID" in m:
        rep.plates = (m.groupby("PlateID")
                      .agg(measurements=("SampleID", "size"), warnings=("warning", "sum"),
                           empty=("empty", "sum"))
                      .reset_index())

    if "LOD" in df.columns:
        below = (df[value] < df["LOD"]).where(df[value].notna() & df["LOD"].notna())
        p = (df.assign(_below=below)
             .groupby(protein)
             .agg(below_lod=("_below", "mean"), missing=(value, lambda s: s.isna().mean()))
             .reset_index()
             .sort_values("below_lod", ascending=False))
        if "Assay" in df.columns and protein != "Assay":
            p = p.merge(df[[protein, "Assay"]].drop_duplicates(protein), on=protein, how="left")
        rep.proteins = p
        mostly = p[p["below_lod"] > max_below_lod]
        rep.lines.append(f"proteins below the detection limit in more than {max_below_lod:.0%} of "
                         f"measurements: {len(mostly)} of {len(p)}")
    else:
        rep.lines.append("no LOD column: detection-limit check skipped")

    if group:
        if group not in df.columns:
            raise KeyError(f"Group column '{group}' not found (did you run add_sample_info?).")
        if "PlateID" not in df.columns:
            rep.lines.append("no PlateID column: group-vs-plate check skipped")
        else:
            per_sample = df.drop_duplicates(["SampleID", "PlateID"])
            plates_of = per_sample.groupby("SampleID")["PlateID"].agg(
                lambda s: "both/several plates" if s.nunique() > 1 else s.iloc[0])
            # samples without sample information (no match in the sample list) are not a group
            groups = df.drop_duplicates("SampleID").set_index("SampleID")[group].dropna()
            table = pd.crosstab(groups, plates_of.reindex(groups.index), dropna=False)
            rep.group_by_plate = table
            single = [c for c in table.columns if c != "both/several plates"]
            missing = [(g, pl) for g in table.index for pl in single
                       if table.loc[g, pl] == 0 and table[pl].sum() >= 5]
            if missing:
                rep.lines.append(
                    "WARNING group mixed up with plate: " +
                    "; ".join(f"no '{g}' measured only on {pl}" for g, pl in missing) +
                    ". Plate differences can look like group differences — bridge/normalise the "
                    "plates and check before comparing groups.")
            else:
                rep.lines.append(f"'{group}': every group is present on every plate")
        if "LOD" in df.columns:
            by_group = below_lod_by_group(df, group)
            rep.below_lod_by_group = by_group
            big = by_group[by_group["difference"] >= 0.30]
            if len(big):
                rep.lines.append(
                    f"detectability differs between '{group}' groups (>= 30 percentage points "
                    f"below LOD): {', '.join(map(str, big.index[:10]))}"
                    + (f" ... ({len(big)} in total)" if len(big) > 10 else "")
                    + " - keep these: being undetectable in one group can itself be a signal")

    if verbose:
        print(rep)
    return rep


def apply_qc(npx: pd.DataFrame, *, drop_controls: bool = True, drop_warnings: bool = True,
             drop_empty: bool = True, max_below_lod: float | None = None,
             repeats: str = "prefer_pass_then_last", drop_samples=None,
             show_ids: bool = True, verbose: bool = True) -> pd.DataFrame:
    """Remove low-quality data — every choice visible, everything removed reported.

    Parameters
    ----------
    drop_controls : remove control samples (``SampleType`` other than SAMPLE: negative, plate
        and sample controls) and control assays (``AssayType`` other than assay).
    drop_warnings : remove measurements with a quality warning.
    drop_empty : remove measurements without any values.
    max_below_lod : remove proteins below the detection limit in more than this share of the
        measurements (0-1). Default ``None``: keep all proteins - a protein below LOD in one
        group but not in another can be informative, and uninformative proteins do not come out
        of the biomarker search anyway. Look at ``qc_report(..., group=...)`` first.
    repeats : what to do with a sample measured more than once on the same plate:
        ``"prefer_pass_then_last"`` (default: keep a passed measurement, and of those the last —
        usually the rerun), ``"first"``, ``"last"`` or ``"keep_all"``.
        Bridging samples (same sample on *different* plates) are always kept: they are needed
        to line up the plates.
    drop_samples : sample IDs to remove that you found yourself (e.g. an outlier in
        ``plot_sample_qc`` or ``bridge_check``). A list removes every measurement of those
        samples; a dict ``{sample ID: PlateID}`` removes only the measurement on that plate.
        IDs are compared after ``tidy_sample_ids``.
    show_ids : list sample IDs in the report; ``False`` gives counts only (for sharing).

    Returns the cleaned table; the list of removals is in ``df.attrs["qc_applied"]``.
    """
    value, protein = _value_col(npx), _protein_col(npx)
    keys = _measurement_keys(npx)
    df, done = npx.copy(), []

    def drop_measurements(mask_df: pd.DataFrame, why: str):
        nonlocal df
        if mask_df.empty:
            return
        marked = df[keys].merge(mask_df[keys].drop_duplicates().assign(_drop=True),
                                on=keys, how="left")
        df = df[~marked["_drop"].notna().to_numpy()]
        done.append(f"removed {len(mask_df)} measurements ({why}): {_ids(mask_df['SampleID'], show_ids)}")

    if drop_controls:
        ctrl_samples, ctrl_assays = _control_masks(df)
        if ctrl_samples.any():
            ctrl = df[ctrl_samples].drop_duplicates(keys)
            kinds = ctrl["SampleType"].value_counts()
            done.append(f"removed {len(ctrl)} control measurements: "
                        + ", ".join(f"{n} {k}" for k, n in kinds.items()))
        if ctrl_assays.any():
            done.append(f"removed {df.loc[ctrl_assays, protein].nunique()} control assays: "
                        + ", ".join(sorted(df.loc[ctrl_assays, "AssayType"].astype(str).unique())))
        df = df[~(ctrl_samples | ctrl_assays)]

    if drop_samples is not None:
        from .reader import tidy_sample_ids
        tidy = tidy_sample_ids(df["SampleID"])
        if isinstance(drop_samples, dict):
            wanted = {(t, str(p)) for t, p in zip(tidy_sample_ids(pd.Series(list(drop_samples), dtype="string")),
                                                  drop_samples.values())}
            plates = df["PlateID"].astype(str) if "PlateID" in df.columns else pd.Series("", index=df.index)
            hit = pd.Series([(t, p) in wanted for t, p in zip(tidy, plates)], index=df.index)
        else:
            hit = tidy.isin(set(tidy_sample_ids(pd.Series(list(drop_samples), dtype="string"))))
        hit = hit.fillna(False).astype(bool)
        if hit.any():
            drop_measurements(df[hit].drop_duplicates(keys), "chosen by you (drop_samples)")
        else:
            done.append("drop_samples: none of the given samples was found")

    if drop_empty:
        m = _measurement_table(df)
        drop_measurements(m[m["empty"]], "no values")
    if drop_warnings and "QC_Warning" in df.columns:
        m = _measurement_table(df)
        drop_measurements(m[m["warning"]], "quality warning")

    measure = _measure_col(df)
    if repeats != "keep_all" and measure:
        m = _measurement_table(df)
        rep = m[m["repeat"]].copy()
        if not rep.empty:
            per_plate = [c for c in ("SampleID", "Panel", "PlateID") if c in rep.columns]
            if repeats == "prefer_pass_then_last":
                rep = rep.sort_values(per_plate + ["warning", measure],
                                      ascending=[True] * len(per_plate) + [False, True])
                keep = rep.groupby(per_plate).tail(1)
            elif repeats in ("first", "last"):
                rep = rep.sort_values(per_plate + [measure])
                keep = rep.groupby(per_plate).head(1) if repeats == "first" else rep.groupby(per_plate).tail(1)
            else:
                raise ValueError("repeats must be 'prefer_pass_then_last', 'first', 'last' or 'keep_all'")
            drop = rep.merge(keep[keys].assign(_k=True), on=keys, how="left")
            drop_measurements(drop[drop["_k"].isna()], f"repeat on the same plate, rule '{repeats}'")

    if max_below_lod is not None and "LOD" in df.columns:
        below = (df[value] < df["LOD"]).where(df[value].notna() & df["LOD"].notna())
        share = df.assign(_b=below).groupby(protein)["_b"].mean()
        gone = share[share > max_below_lod].index
        if len(gone):
            names = df.loc[df[protein].isin(gone), "Assay" if "Assay" in df else protein].unique()
            df = df[~df[protein].isin(gone)]
            done.append(f"removed {len(gone)} proteins below the detection limit in more than "
                        f"{max_below_lod:.0%} of measurements: {', '.join(map(str, names))}")

    kept = df.groupby(keys, dropna=False).ngroups
    done.append(f"kept: {kept} measurements, {df['SampleID'].nunique()} sample IDs, "
                f"{df[protein].nunique()} proteins")
    df = df.reset_index(drop=True)
    df.attrs = dict(npx.attrs)
    df.attrs["qc_applied"] = done
    if verbose:
        print("Quality steps applied\n  " + "\n  ".join(done))
    return df
