"""Read Olink NPX exports into one standard table.

No AI and no internet: the reader knows the vocabulary of Olink exports (column names and
their variants between products and software versions) and applies fixed rules. The same file
therefore always gives the same table, and every change it makes is written to a reading report.

Two layouts are recognised automatically:

- **long** (current exports): one row per sample and protein;
- **wide** (older exports, e.g. NPX Manager): a header block (Panel / Assay / Uniprot ID /
  OlinkID), one row per sample with all proteins side by side, plate and QC-warning columns per
  panel, and a bottom block (LOD, Missing Data freq., Normalization, Plate LOD, ...).
  It is converted to the long layout, so everything after reading works the same.
"""

from __future__ import annotations

import csv
from pathlib import Path

import pandas as pd

# ---------------------------------------------------------------------------------------------
# Vocabulary
# ---------------------------------------------------------------------------------------------

# Standard column name -> names seen in Olink exports (compared ignoring case, spaces and "_").
COLUMN_VARIANTS = {
    "SampleID": ["SampleID", "Sample_ID", "Sample ID"],
    "SampleType": ["SampleType", "Sample_Type"],
    "AssayType": ["AssayType", "Assay_Type"],
    "Index": ["Index"],
    "OlinkID": ["OlinkID", "Olink_ID", "OID"],
    "UniProt": ["UniProt", "Uniprot_ID", "UniProtID"],
    "Assay": ["Assay", "Assay_Name"],
    "Panel": ["Panel"],
    "Panel_Version": ["Panel_Version", "Panel_Lot_Nr"],
    "PlateID": ["PlateID", "Plate_ID", "Plate ID"],
    "QC_Warning": ["QC_Warning", "SampleQC", "Sample_QC", "QC Warning"],
    "Assay_Warning": ["Assay_Warning", "AssayQC", "Assay_QC", "Assay warning"],
    "LOD": ["LOD", "LODNPX", "LOD_NPX", "PlateLOD", "Plate_LOD"],
    "MaxLOD": ["MaxLOD", "Max_LOD"],
    "NPX": ["NPX"],
    "MissingFreq": ["MissingFreq", "Missing_Freq"],
    "Normalization": ["Normalization"],
}

REQUIRED_ONE_OF = ["Assay", "OlinkID"]  # at least one protein identifier

NUMERIC = ["NPX", "Ct", "Quantified_value", "LOD", "PlateLOD", "MaxLOD", "MissingFreq",
           "LLOQ", "ULOQ", "Lowest_quantifiable_level", "Adj_factor", "Index",
           "QC Deviation Inc Ctrl", "QC Deviation Det Ctrl"]

# Words that mean "no value" (compared in capitals).
MISSING_WORDS = {"", "NA", "N/A", "#N/A", "NAN", "NONE", "NULL", "NO DATA", "-", "."}

# "True"/"False" in the languages spreadsheet programs write them in.
TRUE_WORDS = {"true", "vero", "wahr", "vrai", "verdadero", "verdadeiro", "waar", "sant", "sand", "prawda"}
FALSE_WORDS = {"false", "falso", "falsch", "faux", "onwaar", "falskt", "falsk", "fałsz", "falsz"}

# --- wide layout ---
WIDE_TOP = {"panel": "Panel", "assay": "Assay", "uniprotid": "UniProt", "olinkid": "OlinkID",
            "unit": "Unit"}
WIDE_DATA_TYPES = {"npxdata": "NPX", "ctdata": "Ct", "quantifieddata": "Quantified_value"}
# bottom-row label -> (standard column, one row per plate?)
WIDE_BOTTOM = {
    "lod": ("LOD", False),
    "maxlod": ("MaxLOD", False),
    "missingdatafreq": ("MissingFreq", False),
    "normalization": ("Normalization", False),
    "lloq": ("LLOQ", False),
    "uloq": ("ULOQ", False),
    "platelod": ("PlateLOD", True),
    "assaywarning": ("Assay_Warning", True),
    "lowestquantifiablelevel": ("Lowest_quantifiable_level", True),
}
WIDE_INTERNAL_CONTROLS = {"incctrl", "incctrl1", "incctrl2", "detctrl", "extctrl"}


def _key(name) -> str:
    return str(name).strip().lower().replace("_", "").replace(" ", "")


def _wkey(value) -> str:
    """Key for labels in wide files ('Missing Data freq.' -> 'missingdatafreq')."""
    if pd.isna(value):
        return ""
    return _key(value).replace(".", "").replace("-", "")


# ---------------------------------------------------------------------------------------------
# Reading the file
# ---------------------------------------------------------------------------------------------

def _detect_separator(path: Path) -> str:
    with open(path, encoding="utf-8-sig", errors="replace") as f:
        sample = f.read(50_000)
    lines = [line for line in sample.splitlines() if line.strip()][:20]
    # The separator that appears most often in the first lines (robust for wide files, whose
    # first rows are mostly empty cells).
    counts = {sep: sum(line.count(sep) for line in lines) for sep in [",", ";", "\t", "|"]}
    best = max(counts, key=counts.get)
    if counts[best] == 0:
        try:
            return csv.Sniffer().sniff(sample, delimiters=",;\t|").delimiter
        except csv.Error:
            return ","
    return best


def _read_grid(path: Path, sep: str | None, report: list[str]) -> pd.DataFrame | None:
    """Read a text or Excel file as a grid of text cells, without assuming a header row."""
    suffix = path.suffix.lower()
    if suffix in (".xlsx", ".xls"):
        report.append(f"file type: {suffix}")
        return pd.read_excel(path, header=None, dtype=str)
    sep = sep or _detect_separator(path)
    report.append(f"file type: {suffix or 'text'}")
    for encoding in ("utf-8-sig", "latin-1"):
        try:
            return pd.read_csv(path, sep=sep, header=None, dtype=str, encoding=encoding,
                               keep_default_na=True, skip_blank_lines=False)
        except UnicodeDecodeError:
            continue
    raise ValueError(f"Could not read {path.name}: unknown text encoding.")


def _is_wide(grid: pd.DataFrame) -> bool:
    first = set(grid.iloc[:25, 0].map(_wkey))
    return {"assay", "olinkid"} <= first


def _grid_to_long_table(grid: pd.DataFrame) -> pd.DataFrame:
    header = [str(c).strip() if pd.notna(c) else f"unnamed_{i}" for i, c in enumerate(grid.iloc[0])]
    df = grid.iloc[1:].reset_index(drop=True)
    df.columns = header
    return df.dropna(how="all")


# ---------------------------------------------------------------------------------------------
# Wide layout -> long layout
# ---------------------------------------------------------------------------------------------

def _wide_to_long(grid: pd.DataFrame, report: list[str]) -> pd.DataFrame:
    grid = grid.reset_index(drop=True)
    labels = grid.iloc[:, 0].map(_wkey)

    top_rows = {}
    for i in range(min(25, len(grid))):
        if labels[i] in WIDE_TOP and WIDE_TOP[labels[i]] not in top_rows:
            top_rows[WIDE_TOP[labels[i]]] = i
    last_top = max(top_rows.values())

    data_type = "NPX"
    for i in range(min(top_rows.values())):
        if labels[i] in WIDE_DATA_TYPES:
            data_type = WIDE_DATA_TYPES[labels[i]]
    value_col = data_type

    # rows below the header block: samples first, then the bottom block
    below = range(last_top + 1, len(grid))
    bottom_rows = [i for i in below if labels[i] in WIDE_BOTTOM]
    first_bottom = min(bottom_rows) if bottom_rows else len(grid)
    sample_rows = [i for i in below if i < first_bottom and labels[i] != ""]

    def top(name, j):
        row = top_rows.get(name)
        if row is None:
            return pd.NA
        value = grid.iat[row, j]
        return str(value).strip() if pd.notna(value) else pd.NA

    # classify the columns
    assay_cols, plate_col, qc_col, dev_cols = [], {}, {}, {}
    n_controls = 0
    for j in range(1, grid.shape[1]):
        label = top("Assay", j)
        if pd.isna(label):
            continue
        key, panel = _wkey(label), top("Panel", j)
        if key == "plateid":
            plate_col[panel] = j
        elif key == "qcwarning":
            qc_col[panel] = j
        elif key == "qcdeviationfrommedian":
            dev_cols[(panel, f"QC Deviation {top('UniProt', j)}")] = j
        elif key in WIDE_INTERNAL_CONTROLS:
            n_controls += 1
        else:
            assay_cols.append(j)

    def panel_column(columns, panel):
        if panel in columns:
            return columns[panel]
        return next(iter(columns.values())) if len(columns) == 1 else None

    samples = grid.iloc[sample_rows].reset_index(drop=True)
    sample_ids = samples[0].astype("string").str.strip()
    index = pd.Series(range(1, len(samples) + 1))

    parts = []
    for j in assay_cols:
        panel = top("Panel", j)
        part = pd.DataFrame({
            "SampleID": sample_ids,
            "Index": index,
            "OlinkID": top("OlinkID", j),
            "UniProt": top("UniProt", j),
            "Assay": top("Assay", j),
            "Panel": panel,
            value_col: samples[j].values,
            "_col": j,
        })
        if "Unit" in top_rows:
            part["Unit"] = top("Unit", j)
        pc, qc = panel_column(plate_col, panel), panel_column(qc_col, panel)
        if pc is not None:
            part["PlateID"] = samples[pc].astype("string").str.strip().values
        if qc is not None:
            part["QC_Warning"] = samples[qc].astype("string").str.strip().values
        for (dev_panel, name), dj in dev_cols.items():
            if dev_panel == panel:
                part[name] = samples[dj].values
        parts.append(part)
    long = pd.concat(parts, ignore_index=True)

    # bottom block: per protein, or per protein and plate
    used, skipped = [], []
    by_std: dict[str, list[int]] = {}
    for i in bottom_rows:
        by_std.setdefault(WIDE_BOTTOM[labels[i]][0], []).append(i)
    for std, rows in by_std.items():
        plate_specific = WIDE_BOTTOM[_wkey(grid.iat[rows[0], 0])][1] or len(rows) > 1
        if not plate_specific:
            values = {j: grid.iat[rows[0], j] for j in assay_cols}
            long[std] = long["_col"].map(values)
            used.append(str(grid.iat[rows[0], 0]).strip())
            continue
        if "PlateID" not in long.columns:
            skipped.append(f"{std} (no plate column)")
            continue
        mapping, unlinked = {}, set()
        for j in assay_cols:
            panel = top("Panel", j)
            pc = panel_column(plate_col, panel)
            if pc is None:
                unlinked.add(panel)
                continue
            plates_in_order = list(dict.fromkeys(samples[pc].dropna().astype(str).str.strip()))
            row_plates = [grid.iat[r, pc] for r in rows]
            if all(pd.notna(p) for p in row_plates):
                pairs = zip(rows, [str(p).strip() for p in row_plates])   # plate written in the row
            elif len(plates_in_order) <= len(rows) and all(
                    pd.isna(grid.iat[r, j]) for r in rows[len(plates_in_order):]):
                pairs = zip(rows, plates_in_order)        # rows in the order the plates appear;
            else:                                         # extra rows (other panels) are empty
                unlinked.add(panel)
                continue
            for r, plate in pairs:
                mapping[(j, plate)] = grid.iat[r, j]
        if not mapping:
            skipped.append(f"{std} ({len(rows)} rows, could not link them to the plates)")
            continue
        long[std] = [mapping.get((j, p)) for j, p in zip(long["_col"], long["PlateID"])]
        used.append(f"{grid.iat[rows[0], 0]} (per plate)")
        if unlinked:
            skipped.append(f"{std} for panel(s) {', '.join(sorted(map(str, unlinked)))} "
                           f"(could not link the rows to the plates)")

    long = long.drop(columns="_col")
    long = long[long[value_col].notna() | long["SampleID"].notna()]

    panels = sorted(long["Panel"].dropna().unique())
    report.append(
        f"layout: wide format (older Olink export, one row per sample), turned into long format: "
        f"{len(samples)} samples x {len(assay_cols)} assays; values: {data_type}"
        + (f"; panels: {', '.join(panels)}" if panels else "")
    )
    if n_controls:
        report.append(f"wide format: {n_controls} internal control assay columns left out")
    if used:
        report.append("wide format: summary rows at the bottom used: " + ", ".join(used))
    if skipped:
        report.append("wide format: summary rows at the bottom skipped: " + ", ".join(skipped))
    if data_type != "NPX":
        report.append(f"note: values are {data_type} (the steps expect NPX); they are in column '{value_col}'")
    return long


# ---------------------------------------------------------------------------------------------
# Cleaning (same for both layouts)
# ---------------------------------------------------------------------------------------------

def _standardise_columns(df: pd.DataFrame, column_map: dict | None, report: list[str]) -> pd.DataFrame:
    df = df.rename(columns=lambda c: str(c).strip())
    if column_map:
        df = df.rename(columns=column_map)
        report.append(f"columns renamed by you: {column_map}")
    lookup = {_key(v): std for std, variants in COLUMN_VARIANTS.items() for v in variants}
    renames = {}
    for col in df.columns:
        std = lookup.get(_key(col))
        if std and std != col and std not in df.columns and std not in renames.values():
            renames[col] = std
    if renames:
        df = df.rename(columns=renames)
        report.append("columns renamed to the names signatomics uses: "
                      + ", ".join(f"{a} -> {b}" for a, b in renames.items()))
    return df


def _value_column(df: pd.DataFrame) -> str:
    for col in ("NPX", "Quantified_value", "Ct"):
        if col in df.columns:
            return col
    return "NPX"


# words that hint at a column's role, for the suggestion when a needed column is missing
HINTS = {
    "SampleID": ("sample", "patient", "subject", "participant", "donor", "id"),
    "Assay": ("assay", "protein", "target", "analyte", "gene", "marker"),
    "NPX": ("npx", "value", "expression", "level", "log2"),
    "LOD": ("lod", "detection", "limit"),
    "PlateID": ("plate",),
    "QC_Warning": ("qc", "warning", "flag"),
}

# optional columns and what they are needed for (listed in the report when missing)
OPTIONAL = {
    "LOD": "the below-LOD checks",
    "PlateID": "bridging and the plate checks",
    "QC_Warning": "the QC status",
}


def _guess(df: pd.DataFrame, needed: list[str]) -> dict:
    """For each needed column, the file column whose name hints at it (first match)."""
    free = [c for c in df.columns if _key(c) not in {_key(s) for s in ("SampleID", "Assay", "OlinkID", "NPX")}]
    guesses = {}
    for std in needed:
        for col in free:
            if col not in guesses and any(h in str(col).lower() for h in HINTS[std]):
                guesses[col] = std
                break
    return guesses


def _check_required(df: pd.DataFrame, path: Path) -> None:
    value = _value_column(df)
    needed = [c for c in ("SampleID",) if c not in df.columns]
    if not any(c in df.columns for c in REQUIRED_ONE_OF):
        needed.append("Assay")
    if value not in df.columns:
        needed.append("NPX")
    if needed:
        names = {"SampleID": "SampleID (the sample ID)", "Assay": "Assay or OlinkID (the protein)",
                 "NPX": "NPX (the value)"}
        guesses = _guess(df, needed)
        hint = (f"  Did you mean: read_npx(\"{path.name}\", column_map={guesses})"
                if guesses else
                f"  Tell the reader which column is which, e.g. "
                f"read_npx(\"{path.name}\", column_map={{'Patient': 'SampleID'}})")
        raise ValueError(
            f"{path.name}: the columns signatomics needs are missing.\n"
            f"  Needed: {', '.join(names[c] for c in needed)}\n"
            f"  Your file has: {', '.join(map(str, df.columns))}\n"
            f"{hint}"
        )


def _to_number(s: pd.Series) -> tuple[pd.Series, list[str]]:
    """Text to numbers: handles decimal commas ("0,5"), percentages ("100%") and "No Data"."""
    notes = []
    text = s.astype("string").str.strip()
    if text.str.contains(",", regex=False).any() and not text.str.contains(".", regex=False).any():
        text = text.str.replace(",", ".", regex=False)
        notes.append("decimal commas")
    is_pct = text.str.endswith("%", na=False)
    if is_pct.any():
        notes.append("percentages (\"12%\") read as fractions (0.12)")
    values = pd.to_numeric(text.str.rstrip("%"), errors="coerce").astype("float64")
    values = values.where(~is_pct.fillna(False).astype(bool), values / 100)
    marked_missing = text.str.upper().isin(MISSING_WORDS).fillna(False)
    if marked_missing.any():
        words = sorted(set(text[marked_missing]) - {""})
        if words:
            n = int(marked_missing.sum())
            notes.append(f"{n} cell{'s' if n != 1 else ''} said {'/'.join(words)}, read as missing values")
    bad = (s.notna() & values.isna() & ~marked_missing).fillna(False)
    if bad.any():
        n = int(bad.sum())
        notes.append(f"{n} cell{'s' if n != 1 else ''} could not be read as numbers, read as missing values")
    return values, notes


def _convert_types(df: pd.DataFrame, report: list[str]) -> pd.DataFrame:
    df["SampleID"] = df["SampleID"].astype("string").str.strip()
    for col in NUMERIC:
        if col in df.columns and not pd.api.types.is_numeric_dtype(df[col]):
            df[col], notes = _to_number(df[col])
            if notes:
                label = "MissingFreq (missing frequency)" if col == "MissingFreq" else col
                report.append(f"{label}: " + "; ".join(notes))
    if "Index" in df.columns and pd.api.types.is_float_dtype(df["Index"]):
        whole = df["Index"].dropna()
        if (whole == whole.round()).all():
            df["Index"] = df["Index"].astype("Int64")
    for col in df.columns:
        if col in NUMERIC or col == "SampleID" or pd.api.types.is_numeric_dtype(df[col]):
            continue
        values = df[col].dropna().astype(str).str.strip().str.lower()
        if values.empty or not set(values) <= TRUE_WORDS | FALSE_WORDS:
            continue
        words = sorted(set(df[col].dropna().astype(str).str.strip()))
        df[col] = df[col].astype(str).str.strip().str.lower().map(
            lambda v: True if v in TRUE_WORDS else (False if v in FALSE_WORDS else pd.NA)
        ).astype("boolean")
        if set(words) - {"True", "False", "TRUE", "FALSE", "true", "false"}:
            report.append(f"{col}: {'/'.join(words)} read as True/False")
    return df


def _tidy_sample_ids(ids: pd.Series) -> pd.Series:
    """Make sample IDs comparable between files (used wherever sample IDs are matched).

    Removes spaces at the ends, uses capitals, and treats spaces, underscores and hyphens the
    same: ``"308 BIS"``, ``"308_bis"`` and ``"308-BIS"`` all become ``"308_BIS"``.
    Leading zeros are kept (``"058"`` stays ``"058"``).
    """
    return (
        ids.astype("string")
        .str.strip()
        .str.upper()
        .str.replace(r"[\s_\-]+", "_", regex=True)
    )


# ---------------------------------------------------------------------------------------------
# Public functions
# ---------------------------------------------------------------------------------------------

def read_npx(path, *, sep: str | None = None, column_map: dict | None = None,
             verbose: bool = True) -> pd.DataFrame:
    """Read an Olink export (csv/txt/tsv, Excel or parquet; long or wide) into a standard table.

    One row per sample and protein. Older *wide* exports (NPX Manager style) are recognised
    automatically and converted. Column names are made standard (e.g. ``PlateLOD`` -> ``LOD``,
    ``SampleQC`` -> ``QC_Warning``), sample IDs are kept as text, numbers written as text
    (decimal commas, percentages, "No Data") become numbers, and true/false words in other
    languages (``VERO``/``FALSO``, ``WAHR``/``FALSCH``, ...) become True/False.

    Parameters
    ----------
    path : file to read.
    sep : column separator for text files; detected automatically when not given.
    column_map : your own renames for unusual files, e.g. ``{"Patient": "SampleID"}``.
    verbose : print the reading report.

    Returns
    -------
    pandas.DataFrame. The reading report (what was changed) is in ``df.attrs["read_report"]``.
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"File not found: {path}")
    report: list[str] = [f"file: {path.name}"]

    if path.suffix.lower() == ".parquet":
        report.append("file type: .parquet")
        df = pd.read_parquet(path)
    else:
        grid = _read_grid(path, sep, report)
        if _is_wide(grid):
            df = _wide_to_long(grid, report)
        else:
            report.append("layout: long format (one row per sample and assay)")
            df = _grid_to_long_table(grid)

    df = _standardise_columns(df, column_map, report)
    _check_required(df, path)
    df = _convert_types(df, report)
    absent = [c for c in OPTIONAL if c not in df.columns]
    if absent:
        report.append("optional columns missing from the file: "
                      + "; ".join(f"{c} (for {OPTIONAL[c]})" for c in absent))
        guesses = _guess(df, absent)
        if guesses:
            report.append(f"  if your file has them under another name: read_npx(..., column_map={guesses})")

    for col in ("QC_Warning", "Assay_Warning"):
        if col in df.columns and df[col].dtype != bool:
            df[col] = df[col].mask(df[col].astype("string").str.strip().str.upper().isin(MISSING_WORDS))
    if "SampleType" in df.columns:
        is_ctrl = df["SampleType"].notna() & (df["SampleType"].astype("string").str.upper() != "SAMPLE")
        keys = [c for c in ("SampleID", "Index", "WellID", "Panel", "PlateID") if c in df.columns]
        n_ctrl = df[is_ctrl.fillna(False).astype(bool)].groupby(keys, dropna=False).ngroups
        if n_ctrl:
            report.append(f"control samples: {n_ctrl} (SampleType other than SAMPLE); apply_qc removes them")
    if "AssayType" in df.columns:
        n_ctrl = df.loc[df["AssayType"].astype("string").str.lower() != "assay", "OlinkID"].nunique()
        if n_ctrl:
            report.append(f"internal control assays: {n_ctrl} (AssayType other than assay); apply_qc removes them")

    runs = df.drop_duplicates([c for c in ("SampleID", "Index", "PlateID") if c in df.columns])
    protein_col = "Assay" if "Assay" in df.columns else "OlinkID"
    n_ids, n_assays, n_runs = df["SampleID"].nunique(), df[protein_col].nunique(), len(runs)
    n_plates = df["PlateID"].nunique() if "PlateID" in df.columns else None
    plates = f", on {n_plates} plates" if n_plates and n_plates > 1 else ""
    def n(k, word):
        return f"{k} {word}" + ("" if k == 1 else "s")
    shape = (f"{n(n_runs, 'sample')} x {n(n_assays, 'assay')}" if n_runs == n_ids
             else f"{n(n_runs, 'measurement')} x {n(n_assays, 'assay')} ({n_ids} samples, some measured twice)")
    report.append(f"result: {n(len(df), 'row')} = {shape}{plates}"
                  if len(df) == n_runs * n_assays
                  else f"result: {n(len(df), 'row')}: {shape}{plates}")

    df = df.reset_index(drop=True)
    df.attrs["read_report"] = report
    if verbose:
        print("Reading report\n  " + "\n  ".join(report))
    return df
