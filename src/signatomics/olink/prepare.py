"""The usual preprocessing in one call: read, quality control, bridging, sample metadata.

``prepare`` runs the same steps you can also run one by one (``read_npx``, ``qc_report``,
``apply_qc``, ``bridge_check``, ``bridge_plates``, ``add_sample_info``), each with its default
choices, and prints one short report. For other choices, run the steps yourself.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from .bridge import bridge_check, bridge_plates
from .qc import apply_qc, qc_report
from .reader import read_npx
from .samples import add_sample_info


def _as_tables(data) -> tuple[list[pd.DataFrame], int]:
    items = data if isinstance(data, (list, tuple)) else [data]
    tables, files = [], 0
    for item in items:
        if isinstance(item, pd.DataFrame):
            tables.append(item)
        elif isinstance(item, (str, Path)):
            tables.append(read_npx(item, verbose=False))
            files += 1
        else:
            raise TypeError("data must be a file path, a table, or a list of them.")
    if not tables:
        raise ValueError("No data given.")
    return tables, files


def _n(count: int, word: str) -> str:
    return f"{count} {word}{'s' if count != 1 else ''}"


def prepare(data, samples=None, *, on: str | None = None, group: str | None = None,
            show_ids: bool = True, verbose: bool = True) -> pd.DataFrame:
    """Preprocess Olink data in one call, with the default choices of each step.

    Parameters
    ----------
    data : an NPX file (path), several NPX files (list of paths, e.g. one per plate), or tables
        from ``read_npx``.
    samples : your sample metadata (a table or a path), optional. Needs ``on``.
    on : the column in your sample metadata that holds the same IDs as ``SampleID`` in the NPX data.
    group : a column of your sample metadata, e.g. ``"diabetes"``: adds the group checks.
    show_ids : list sample IDs in the report; ``False`` gives counts only (for sharing).
    verbose : print the report.

    Steps: read the files → quality control (removes control samples, measurements with QC status
    FAIL, empty measurements, and the extra copy of a sample measured twice on one plate; WARN and
    all assays are kept) → bridging, when there are several plates with bridge samples → your
    sample metadata added.

    Returns the preprocessed table. The short report is in ``df.attrs["prepare_report"]``, every
    step's full report in ``df.attrs["prepare_details"]``. Things that need your decision, such
    as an outlier, are listed under "to look at"; ``prepare`` leaves them to you.
    """
    if samples is not None and on is None:
        raise ValueError("Give on=..., the column in your sample metadata that holds the sample IDs.")
    tables, files = _as_tables(data)
    df = pd.concat(tables, ignore_index=True) if len(tables) > 1 else tables[0].copy()
    protein = "Assay" if "Assay" in df.columns else "OlinkID"
    plates = df["PlateID"].nunique() if "PlateID" in df.columns else 1
    from .reader import _is_control_type
    ctrl = _is_control_type(df["SampleType"]) if "SampleType" in df.columns else df["SampleID"].isna()
    n_ctrl = df.loc[ctrl, "SampleID"].nunique()
    what = (f"{_n(df.loc[~ctrl, 'SampleID'].nunique(), 'sample')}"
            + (f" + {_n(n_ctrl, 'control')}" if n_ctrl else "")
            + f", {_n(df[protein].nunique(), 'assay')}" + (f", {plates} plates" if plates > 1 else ""))
    short = [f"read {_n(files, 'file')}: {what}" if files else f"data: {what}"]
    details, to_check = list(short), []

    # quality control (the group checks look your sample metadata up; nothing is added yet)
    look = qc_report(df, group=group, samples=samples if group else None, on=on,
                     show_ids=show_ids, verbose=False)
    details += ["quality control (found):"] + [f"  {line}" for line in look.lines]
    # what apply_qc removes anyway (empty measurements, controls) is left out of "to look at"
    to_check += [line for line in look.lines
                 if (line.startswith("WARNING")
                     or (line.startswith("outliers:") and not line.startswith("outliers: none"))
                     or "far more often in one" in line)
                 and "empty measurements" not in line and "apply_qc removes" not in line]
    clean = apply_qc(df, show_ids=show_ids, verbose=False)
    removed = clean.attrs["qc_applied"][:-1]
    details += ["quality control (removed):"] + [f"  {line}" for line in removed]
    short += ["quality control: " + line for line in removed] or ["quality control: nothing removed"]

    out = clean
    if plates > 1 and "NPX" in clean.columns:
        check = bridge_check(clean, group=group, samples=samples if group else None, on=on,
                             show_ids=show_ids, verbose=False)
        if check.factors is None:
            short.append("bridging: several plates but no bridge samples, so not bridged")
        else:
            out = bridge_plates(clean, factors=check.factors, verbose=False)
            details += ["bridging (check):"] + [f"  {line}" for line in check.lines]
            details += ["bridging:"] + [f"  {line}" for line in out.attrs["bridge_report"]]
            ref = check.factors["reference"].iloc[0]
            for plate, fp in check.factors.groupby("plate", sort=False):
                short.append(f"bridging: {plate} put on the scale of {ref} "
                             f"({_n(int(fp['adjustment'].notna().sum()), 'assay')} adjusted, "
                             f"{_n(int(fp['bridges'].max()), 'bridge sample')})")
            to_check += [line.strip() for line in check.lines
                         if "shift differently" in line or "differs by more" in line
                         or "WARNING" in line]

    if samples is not None:
        before = dict(out.attrs)
        out = add_sample_info(out, samples, on=on, verbose=False)
        lines = out.attrs["sample_report"]
        details += ["sample metadata:"] + [f"  {line}" for line in lines]
        short += ["sample metadata: " + line for line in lines
                  if line.startswith("matched:") or line.startswith("WARNING")]
        out.attrs = before | dict(out.attrs)

    result = f"result: {_n(out['SampleID'].nunique(), 'sample')}, {_n(out[protein].nunique(), 'assay')}"
    short.append(result)
    details.append(result)
    if to_check:
        short.append("to look at (prepare leaves these to you):")
        short += [f"  {line}" for line in to_check]
    out.attrs = dict(clean.attrs) | dict(out.attrs)
    out.attrs["prepare_report"] = short
    out.attrs["prepare_details"] = details
    if verbose:
        print("Prepared data\n  " + "\n  ".join(short))
    return out
