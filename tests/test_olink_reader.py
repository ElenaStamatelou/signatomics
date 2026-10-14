"""Tests for signatomics.olink.read_npx and the sample-ID tidying it relies on.

Each test writes a tiny made-up file with one known quirk (found in real data, see
analyses/01_PAD000040_first_look.ipynb) and checks that the reader handles it.
"""

import pandas as pd
import pytest

from signatomics.olink import read_npx

HEADER = "SampleID,Index,OlinkID,Assay,PlateID,QC_Warning,LOD,NPX,MissingFreq,AssayFlag"
ROWS = [
    "058,1,OID00471,IL8,plate1,Pass,0.10,4.81,0,FALSO",
    "058,1,OID00472,IL2,plate1,Pass,0.65,0.20,100%,VERO",
    "308_BIS,2,OID00471,IL8,plate1,Warning,0.10,,0.5,FALSO",
    "308_BIS,2,OID00472,IL2,plate1,Warning,0.65,,100%,VERO",
]


def write(tmp_path, text, name="npx.csv"):
    path = tmp_path / name
    path.write_text(text, encoding="utf-8")
    return path


def test_reads_comma_file_and_reports(tmp_path):
    df = read_npx(write(tmp_path, "\n".join([HEADER, *ROWS])), verbose=False)
    assert len(df) == 4
    assert df["SampleID"].nunique() == 2
    report = " ".join(df.attrs["read_report"])
    assert "comma" in report and "2 sample IDs" in report and "2 proteins" in report


def test_sample_ids_stay_text_with_leading_zeros(tmp_path):
    df = read_npx(write(tmp_path, "\n".join([HEADER, *ROWS])), verbose=False)
    assert set(df["SampleID"]) == {"058", "308_BIS"}


def test_percentages_and_decimals_become_numbers(tmp_path):
    df = read_npx(write(tmp_path, "\n".join([HEADER, *ROWS])), verbose=False)
    assert df["MissingFreq"].tolist() == [0.0, 1.0, 0.5, 1.0]
    assert df["NPX"].isna().sum() == 2          # empty values stay empty


def test_italian_true_false_become_booleans(tmp_path):
    df = read_npx(write(tmp_path, "\n".join([HEADER, *ROWS])), verbose=False)
    assert df["AssayFlag"].tolist() == [False, True, False, True]
    assert any("VERO" in line for line in df.attrs["read_report"])


def test_column_with_only_whole_percentages(tmp_path):
    # Real case: the separate plate files of PAD000040 contain only "0%" ... "100%"
    text = "SampleID,Assay,NPX,MissingFreq\nA1,IL8,4.8,0%\nA1,IL2,0.2,100%\nA1,IL4,1.1,\n"
    df = read_npx(write(tmp_path, text), verbose=False)
    assert df["MissingFreq"].tolist()[:2] == [0.0, 1.0]
    assert pd.isna(df["MissingFreq"].iloc[2])


def test_no_data_words_become_empty_and_are_reported(tmp_path):
    # Real case: PAD000040 plate 2 writes "No Data" for two empty measurements
    # ("#N/A", "NA" and similar are already made empty by pandas itself while reading)
    text = "SampleID,Assay,NPX\nA1,IL8,4.8\nA2,IL8,No Data\nA3,IL8,#N/A\nA4,IL8,oops\n"
    df = read_npx(write(tmp_path, text), verbose=False)
    assert df["NPX"].isna().tolist() == [False, True, True, True]
    report = " ".join(df.attrs["read_report"])
    assert "1 values written as No Data" in report
    assert "1 values not readable as numbers" in report


def test_semicolon_file_with_decimal_commas(tmp_path):
    text = "SampleID;Assay;NPX;LOD\nA1;IL8;4,81;0,10\nA2;IL8;5,02;0,10\n"
    df = read_npx(write(tmp_path, text), verbose=False)
    assert df["NPX"].tolist() == [4.81, 5.02]
    assert "semicolon" in " ".join(df.attrs["read_report"])


def test_column_name_variants_are_standardised(tmp_path):
    text = "Sample_ID,Assay,PlateLOD,SampleQC,NPX\nA1,IL8,0.1,PASS,4.8\n"
    df = read_npx(write(tmp_path, text), verbose=False)
    assert {"SampleID", "LOD", "QC_Warning"} <= set(df.columns)


def test_own_column_names_via_column_map(tmp_path):
    text = "Patient,Assay,NPX\nA1,IL8,4.8\n"
    df = read_npx(write(tmp_path, text), column_map={"Patient": "SampleID"}, verbose=False)
    assert df["SampleID"].tolist() == ["A1"]


def test_clear_message_when_not_an_npx_file(tmp_path):
    text = "Patient,Glucose\nA1,5.6\n"
    with pytest.raises(ValueError, match="Missing: SampleID, NPX"):
        read_npx(write(tmp_path, text), verbose=False)


def test_parquet_and_excel_give_the_same_table(tmp_path):
    csv_df = read_npx(write(tmp_path, "\n".join([HEADER, *ROWS])), verbose=False)
    raw = pd.read_csv(tmp_path / "npx.csv", dtype=str)
    raw.to_parquet(tmp_path / "npx.parquet")
    raw.to_excel(tmp_path / "npx.xlsx", index=False)
    for name in ("npx.parquet", "npx.xlsx"):
        other = read_npx(tmp_path / name, verbose=False)
        pd.testing.assert_series_equal(other["NPX"], csv_df["NPX"])
        assert other["SampleID"].tolist() == csv_df["SampleID"].tolist()


def test_missing_file_gives_clear_error(tmp_path):
    with pytest.raises(FileNotFoundError):
        read_npx(tmp_path / "does_not_exist.csv")


def test_tidy_sample_ids():
    from signatomics.olink.reader import _tidy_sample_ids   # internal: used wherever IDs are matched
    ids = pd.Series(["308 BIS", " 308_bis", "308-BIS", "058", "a  b"])
    assert _tidy_sample_ids(ids).tolist() == ["308_BIS", "308_BIS", "308_BIS", "058", "A_B"]
