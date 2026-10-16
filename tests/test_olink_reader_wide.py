"""Tests for reading older *wide* Olink exports (NPX Manager style).

The layout (studied from Olink's documentation of the format, re-created here with made-up
values):

    <project>        <software version>
    NPX Data
    Panel            Panel A   Panel A   Panel B   Panel A    Panel A     Panel B    Panel B
    Assay            IL8       IL2       CCL3      Plate ID   QC Warning  Plate ID   QC Warning
    Uniprot ID       P10145    P60568    P10147
    OlinkID          OID001    OID002    OID101
    (empty)
    S1               4.8       0.2       7.1       plateA1    Pass        plateB1    Pass
    ...
    (empty)
    LOD              0.1       0.6       0.3
    Missing Data freq.  0      0.5       0
"""

import numpy as np
import pandas as pd
import pytest

from signatomics.olink import read_npx

ASSAYS = [  # panel, assay, uniprot, olinkid
    ("Panel A", "IL8", "P10145", "OID001"),
    ("Panel A", "IL2", "P60568", "OID002"),
    ("Panel B", "CCL3", "P10147", "OID101"),
]
SAMPLES = [  # id, values per assay, plate A, qc A, plate B, qc B
    ("S1", ["4.8", "0.2", "7.1"], "plateA1", "Pass", "plateB1", "Pass"),
    ("S2", ["5.0", "0.9", "6.9"], "plateA1", "Warning", "plateB1", "Pass"),
    ("058", ["4.1", "NaN", "7.4"], "plateA2", "Pass", "plateB1", "Pass"),
    ("S1", ["4.7", "0.3", "7.0"], "plateA2", "Pass", "plateB1", "Pass"),  # repeat of S1
]


def wide_grid(data_type="NPX Data", bottom=None, controls=False):
    """Build a wide Olink export as a grid of text cells."""
    extra = [("Panel A", "Plate ID"), ("Panel A", "QC Warning"),
             ("Panel B", "Plate ID"), ("Panel B", "QC Warning")]
    if controls:
        extra.insert(0, ("Panel A", "Inc Ctrl 1"))
    n_cols = 1 + len(ASSAYS) + len(extra)
    blank = [""] * n_cols
    rows = [["my_project", "NPX Manager 2.1"] + [""] * (n_cols - 2), [data_type] + blank[1:]]
    rows.append(["Panel"] + [a[0] for a in ASSAYS] + [e[0] for e in extra])
    rows.append(["Assay"] + [a[1] for a in ASSAYS] + [e[1] for e in extra])
    rows.append(["Uniprot ID"] + [a[2] for a in ASSAYS] + [""] * len(extra))
    rows.append(["OlinkID"] + [a[3] for a in ASSAYS] + [""] * len(extra))
    rows.append(blank)
    for sid, values, pa, qa, pb, qb in SAMPLES:
        ctrl = ["11.2"] if controls else []
        rows.append([sid] + values + ctrl + [pa, qa, pb, qb])
    rows.append(blank)
    for label, values in (bottom or [("LOD", ["0.1", "0.6", "0.3"]),
                                     ("Missing Data freq.", ["0", "0.5", "0"])]):
        rows.append([label] + values + [""] * len(extra))
    return rows


def write_csv(tmp_path, rows, sep=";", name="wide.csv"):
    path = tmp_path / name
    path.write_text("\n".join(sep.join(r) for r in rows) + "\n", encoding="utf-8")
    return path


def value(df, sid_index, olinkid, col="NPX"):
    return df.loc[(df["Index"] == sid_index) & (df["OlinkID"] == olinkid), col].iloc[0]


def test_wide_is_recognised_and_converted(tmp_path):
    df = read_npx(write_csv(tmp_path, wide_grid()), verbose=False)
    assert len(df) == len(SAMPLES) * len(ASSAYS)               # one row per sample and protein
    assert any("layout: wide" in line for line in df.attrs["read_report"])
    assert value(df, 1, "OID001") == 4.8
    assert value(df, 3, "OID101") == 7.4
    assert np.isnan(value(df, 3, "OID002"))                    # "NaN" cell stays empty


def test_each_panel_gets_its_own_plate_and_qc(tmp_path):
    df = read_npx(write_csv(tmp_path, wide_grid()), verbose=False)
    assert value(df, 2, "OID001", "QC_Warning") == "Warning"   # Panel A column
    assert value(df, 2, "OID101", "QC_Warning") == "Pass"      # Panel B column
    assert value(df, 3, "OID001", "PlateID") == "plateA2"
    assert value(df, 3, "OID101", "PlateID") == "plateB1"


def test_repeated_sample_keeps_separate_measurements(tmp_path):
    df = read_npx(write_csv(tmp_path, wide_grid()), verbose=False)
    s1 = df[(df["SampleID"] == "S1") & (df["OlinkID"] == "OID001")]
    assert sorted(s1["NPX"]) == [4.7, 4.8]
    assert s1["Index"].nunique() == 2


def test_bottom_rows_become_columns(tmp_path):
    df = read_npx(write_csv(tmp_path, wide_grid()), verbose=False)
    assert value(df, 1, "OID002", "LOD") == 0.6
    assert value(df, 4, "OID002", "MissingFreq") == 0.5


def test_plate_lod_rows_are_linked_to_plates(tmp_path):
    # Two "Plate LOD" rows: Panel A has two plates (in order of appearance), Panel B only one
    # (its second row is empty).
    bottom = [("Plate LOD", ["0.1", "0.6", "0.3"]), ("Plate LOD", ["0.2", "0.7", ""])]
    df = read_npx(write_csv(tmp_path, wide_grid(bottom=bottom)), verbose=False)
    assert value(df, 1, "OID002", "LOD") == 0.6               # Panel A, plateA1
    assert value(df, 3, "OID002", "LOD") == 0.7               # Panel A, plateA2
    assert value(df, 3, "OID101", "LOD") == 0.3               # Panel B, plateB1


def test_plate_lod_that_cannot_be_linked_is_reported(tmp_path):
    # Three rows but Panel A has two plates and the third row is not empty -> ambiguous
    bottom = [("Plate LOD", ["0.1", "0.6", "0.3"]), ("Plate LOD", ["0.2", "0.7", ""]),
              ("Plate LOD", ["0.3", "0.8", ""])]
    df = read_npx(write_csv(tmp_path, wide_grid(bottom=bottom)), verbose=False)
    assert any("skipped" in line and "Panel A" in line for line in df.attrs["read_report"])


def test_internal_control_columns_are_left_out(tmp_path):
    df = read_npx(write_csv(tmp_path, wide_grid(controls=True)), verbose=False)
    assert set(df["Assay"]) == {"IL8", "IL2", "CCL3"}
    assert any("internal control assay columns left out" in line for line in df.attrs["read_report"])


def test_sample_ids_stay_text(tmp_path):
    df = read_npx(write_csv(tmp_path, wide_grid()), verbose=False)
    assert "058" in set(df["SampleID"])


def test_quantified_data_goes_to_its_own_column(tmp_path):
    df = read_npx(write_csv(tmp_path, wide_grid(data_type="Quantified Data")), verbose=False)
    assert "Quantified_value" in df.columns and "NPX" not in df.columns
    assert any("expect NPX" in line for line in df.attrs["read_report"])


@pytest.mark.parametrize("fmt", ["csv_semicolon", "csv_comma", "excel"])
def test_same_table_from_every_file_type(tmp_path, fmt):
    rows = wide_grid()
    if fmt == "excel":
        path = tmp_path / "wide.xlsx"
        pd.DataFrame(rows).replace("", None).to_excel(path, header=False, index=False)
    else:
        path = write_csv(tmp_path, rows, sep=";" if fmt == "csv_semicolon" else ",")
    df = read_npx(path, verbose=False)
    assert value(df, 4, "OID001") == 4.7
    assert value(df, 2, "OID001", "QC_Warning") == "Warning"


def test_round_trip_long_to_wide_and_back(tmp_path):
    """Write a known long table in the wide layout, read it back: nothing may change."""
    rng = np.random.default_rng(0)
    proteins = [("Panel A", f"P{i}", f"U{i}", f"OID{i:03d}") for i in range(12)]
    ids = [f"{n:03d}" for n in range(20)]
    truth = pd.DataFrame(
        [(sid, k + 1, oid, round(float(rng.normal(5, 2)), 6), "plate1" if k < 10 else "plate2")
         for k, sid in enumerate(ids) for (_, _, _, oid) in proteins],
        columns=["SampleID", "Index", "OlinkID", "NPX", "PlateID"],
    )
    header = ["Panel", "Assay", "Uniprot ID", "OlinkID"]
    rows = [["proj", "v1"] + [""] * (len(proteins)), ["NPX Data"] + [""] * (len(proteins) + 1)]
    for h, pos in zip(header, range(4)):
        rows.append([h] + [p[pos] for p in proteins] + ["Plate ID" if h == "Assay" else ("Panel A" if h == "Panel" else "")])
    rows.append([""] * (len(proteins) + 2))
    for k, sid in enumerate(ids):
        vals = truth[truth["Index"] == k + 1].sort_values("OlinkID")["NPX"].astype(str).tolist()
        rows.append([sid] + vals + ["plate1" if k < 10 else "plate2"])
    df = read_npx(write_csv(tmp_path, rows), verbose=False)
    merged = truth.merge(df, on=["Index", "OlinkID"], suffixes=("", "_read"))
    assert len(merged) == len(truth)
    assert (merged["SampleID"] == merged["SampleID_read"]).all()
    assert np.allclose(merged["NPX"], merged["NPX_read"])
    assert (merged["PlateID"] == merged["PlateID_read"]).all()
