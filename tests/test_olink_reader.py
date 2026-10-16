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
    assert "file type: .csv" in report and "2 samples x 2 assays" in report


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
    assert "1 cell said No Data" in report
    assert "1 cell could not be read as numbers" in report


def test_semicolon_file_with_decimal_commas(tmp_path):
    text = "SampleID;Assay;NPX;LOD\nA1;IL8;4,81;0,10\nA2;IL8;5,02;0,10\n"
    df = read_npx(write(tmp_path, text), verbose=False)
    assert df["NPX"].tolist() == [4.81, 5.02]
    assert "decimal commas" in " ".join(df.attrs["read_report"])


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
    with pytest.raises(ValueError, match="Needed: SampleID .*Assay or OlinkID .*NPX") as err:
        read_npx(write(tmp_path, text), verbose=False)
    assert "'Patient': 'SampleID'" in str(err.value)          # a suggestion to copy


def test_sample_id_column_holding_something_else(tmp_path):
    rows = [f"{i},A{i % 2 + 1},{assay},4.{i}" for i, assay in enumerate(["IL8", "IL6", "IL10"] * 2)]
    text = "SampleID,SampleName,Assay,NPX\n" + "\n".join(rows) + "\n"
    path = write(tmp_path, text)
    report = read_npx(path, verbose=False).attrs["read_report"]
    assert any("SampleID may hold something else" in line and "SampleName" in line for line in report)
    fixed = read_npx(path, column_map={"SampleName": "SampleID"}, verbose=False)
    assert sorted(fixed["SampleID"].unique()) == ["A1", "A2"] and "SampleID_in_file" in fixed.columns


def test_damaged_numbers_and_group_sample_type(tmp_path):
    text = ("SampleID,Sample_Type,Assay,NPX\nCONTROL_SAMPLE,AD,IL8,1.205.603\n"
            "A1,AD,IL8,912.855\nA2,HC,IL8,4.8\n")
    df = read_npx(write(tmp_path, text), verbose=False)
    report = " ".join(df.attrs["read_report"])
    assert "thousands separator" in report and "Sample_Type kept as it is" in report
    assert "Sample_Type" in df.columns and "SampleType" not in df.columns


def test_missing_optional_columns_are_reported_with_a_suggestion(tmp_path):
    text = "SampleID,Assay,NPX,Plate\nA1,IL8,4.8,p1\n"
    report = read_npx(write(tmp_path, text), verbose=False).attrs["read_report"]
    assert any(line.startswith("optional columns missing from the file: LOD") for line in report)
    assert any("{'Plate': 'PlateID'}" in line for line in report)


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
