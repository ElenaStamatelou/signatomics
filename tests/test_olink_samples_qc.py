"""Tests for add_sample_info, qc_report and apply_qc on a small made-up dataset.

The made-up data copies the situations found in PAD000040: an ID written differently in the
two files, a completely empty measurement, a sample rerun on the same plate (first with a
warning), a bridging sample on both plates, a protein mostly below the detection limit, and a
group that was never measured on one of the plates.
"""

import numpy as np
import pandas as pd
import pytest

from signatomics.olink import add_sample_info, apply_qc, qc_report

PROTEINS = {"P1": 0.0, "P2": 0.0, "LOW": 9.0}  # LOW: LOD so high that most values are below it


def measurement(sid, index, plate, qc="Pass", empty=False, low_above=False):
    rows = []
    for p, lod in PROTEINS.items():
        npx = np.nan if empty else (10.0 if (p == "LOW" and low_above) else 5.0)
        rows.append(dict(SampleID=sid, Index=index, Assay=p, OlinkID=f"OID_{p}", PlateID=plate,
                         QC_Warning=qc, LOD=lod, NPX=npx))
    return rows


@pytest.fixture
def npx():
    rows = []
    rows += measurement("A1", 1, "plate1", low_above=True)
    rows += measurement("A2", 2, "plate1")
    rows += measurement("A3", 3, "plate1", qc="Warning")            # rerun below
    rows += measurement("A3", 4, "plate1")                          # the rerun (passes)
    rows += measurement("A4", 5, "plate1", qc="Warning", empty=True)
    rows += measurement("A5", 6, "plate1")
    rows += measurement("A6", 7, "plate1")
    rows += measurement("B1", 1, "plate2")
    rows += measurement("B2", 2, "plate2")
    rows += measurement("308_BIS", 3, "plate2")
    rows += measurement("A1", 4, "plate2")                          # bridging sample
    return pd.DataFrame(rows)


@pytest.fixture
def sample_list():
    return pd.DataFrame({
        "FIF ID": ["A1", "A2", "A3", "A4", "A5", "A6", "B1", "B2", "308 BIS", "C9"],
        "group": ["no", "no", "no", "yes", "no", "no", "pre", "pre", "yes", "no"],
        "PlateID": ["P1"] * 10,                 # clashes with the protein data's PlateID
    })


# --- add_sample_info -----------------------------------------------------------------------

def test_matching_after_tidying_and_report(npx, sample_list):
    df = add_sample_info(npx, sample_list, on="FIF ID", verbose=False)
    assert df["group"].notna().all()                                 # every protein row matched
    report = " ".join(df.attrs["sample_report"])
    assert "1 only after tidying" in report                          # 308_BIS vs 308 BIS
    assert "C9" in report                                            # in list, no protein data


def test_unmatched_protein_samples_are_kept_and_reported(npx, sample_list):
    df = add_sample_info(npx, sample_list[sample_list["FIF ID"] != "B2"], on="FIF ID", verbose=False)
    assert df.loc[df["SampleID"] == "B2", "group"].isna().all()
    assert any("NOT in the sample list" in line and "B2" in line for line in df.attrs["sample_report"])


def test_clashing_columns_get_a_suffix(npx, sample_list):
    df = add_sample_info(npx, sample_list, on="FIF ID", verbose=False)
    assert "PlateID_samplelist" in df.columns
    assert set(df["PlateID"]) == {"plate1", "plate2"}               # original column untouched


def test_duplicate_ids_in_sample_list_are_refused(npx, sample_list):
    dup = pd.concat([sample_list, sample_list.iloc[[0]]])
    with pytest.raises(ValueError, match="more than once"):
        add_sample_info(npx, dup, on="FIF ID", verbose=False)


def test_sample_list_from_file(tmp_path, npx, sample_list):
    path = tmp_path / "samples.xlsx"
    sample_list.to_excel(path, index=False)
    df = add_sample_info(npx, path, on="FIF ID", verbose=False)
    assert df["group"].notna().all()


# --- qc_report -----------------------------------------------------------------------------

def test_qc_report_finds_everything(npx, sample_list):
    df = add_sample_info(npx, sample_list, on="FIF ID", verbose=False)
    rep = qc_report(df, group="group", verbose=False)
    text = str(rep)
    assert "quality warnings: 2 measurements" in text
    assert "empty measurements (no values at all): A4" in text
    assert "(repeats): A3" in text
    assert "(bridging samples): 1" in text
    assert "1 of 3" in text                                          # LOW mostly below LOD
    assert "no 'pre' measured only on plate1" in text                # group mixed up with plate


def test_qc_report_changes_nothing(npx):
    before = npx.copy()
    qc_report(npx, verbose=False)
    pd.testing.assert_frame_equal(npx, before)


# --- apply_qc ------------------------------------------------------------------------------

def test_apply_qc_defaults(npx):
    df = apply_qc(npx, verbose=False)
    m = df.drop_duplicates(["SampleID", "Index", "PlateID"])
    assert "A4" not in set(m["SampleID"])                            # empty -> removed
    assert m.loc[m["SampleID"] == "A3", "Index"].tolist() == [4]     # warned run removed, rerun kept
    assert (m["SampleID"] == "A1").sum() == 2                        # bridging sample kept twice
    assert "LOW" in set(df["Assay"])                                 # all proteins kept by default
    assert any("kept:" in line for line in df.attrs["qc_applied"])


def test_apply_qc_protein_threshold_only_when_asked(npx):
    df = apply_qc(npx, max_below_lod=0.5, verbose=False)
    assert "LOW" not in set(df["Assay"])
    assert any("removed 1 proteins" in line for line in df.attrs["qc_applied"])


def test_detectability_difference_between_groups_is_reported():
    from signatomics.olink import below_lod_by_group
    # protein X: below LOD in every "healthy" sample, above LOD in every "disease" sample
    rows = []
    for k in range(20):
        g = "healthy" if k < 10 else "disease"
        rows.append(dict(SampleID=f"S{k}", Index=k, Assay="X", PlateID="p1", LOD=1.0,
                         NPX=0.5 if g == "healthy" else 2.0, group=g))
        rows.append(dict(SampleID=f"S{k}", Index=k, Assay="Y", PlateID="p1", LOD=1.0,
                         NPX=3.0, group=g))
    df = pd.DataFrame(rows)
    table = below_lod_by_group(df, "group")
    assert table.loc["X", "healthy"] == 1.0 and table.loc["X", "disease"] == 0.0
    assert table.loc["X", "difference"] == 1.0 and table.loc["Y", "difference"] == 0.0
    text = str(qc_report(df, group="group", verbose=False))
    assert "detectability differs" in text and "X" in text


def test_apply_qc_choices_are_respected(npx):
    df = apply_qc(npx, drop_warnings=False, drop_empty=False, max_below_lod=None,
                  repeats="keep_all", verbose=False)
    assert len(df) == len(npx)                                        # nothing removed


def test_repeat_rule_prefers_passed_measurement(npx):
    # keep warnings, but resolve the repeat: the passed rerun (Index 4) must win
    df = apply_qc(npx, drop_warnings=False, drop_empty=False, max_below_lod=None, verbose=False)
    a3 = df.loc[df["SampleID"] == "A3", "Index"].unique().tolist()
    assert a3 == [4]


def test_apply_qc_does_not_change_the_input(npx):
    before = npx.copy()
    apply_qc(npx, verbose=False)
    pd.testing.assert_frame_equal(npx, before)


# --- control samples and control assays (newer exports) ----------------------------------

def test_controls_are_reported_and_removed(tmp_path):
    from signatomics.olink import read_npx
    text = (
        "SampleID;SampleType;Assay;OlinkID;AssayType;PlateID;QC_Warning;LODNPX;NPX\n"
        "S1;SAMPLE;IL8;OID1;assay;p1;PASS;0,1;4,8\n"
        "S1;SAMPLE;Inc Ctrl 1;OID9;inc_ctrl_1;p1;PASS;0,1;9,0\n"
        "NEG;NEGATIVE_CONTROL;IL8;OID1;assay;p1;PASS;0,1;0,0\n"
        "NEG;NEGATIVE_CONTROL;Inc Ctrl 1;OID9;inc_ctrl_1;p1;PASS;0,1;9,1\n"
    )
    path = tmp_path / "ext.csv"
    path.write_text(text, encoding="utf-8")
    df = read_npx(path, verbose=False)
    report = " ".join(df.attrs["read_report"])
    assert "1 control measurements" in report and "1 control assays" in report
    assert "LODNPX -> LOD" in report                                  # newer LOD column name
    clean = apply_qc(df, verbose=False)
    assert set(clean["SampleID"]) == {"S1"} and set(clean["Assay"]) == {"IL8"}
    kept = apply_qc(df, drop_controls=False, max_below_lod=None, verbose=False)
    assert len(kept) == 4


def test_controls_without_id_and_repeats_by_well(tmp_path):
    from signatomics.olink import read_npx
    # newer exports: no Index, but WellID; a control without SampleID; S1 measured twice
    text = (
        "SampleID;SampleType;WellID;Assay;PlateID;QC_Warning;NPX\n"
        "S1;SAMPLE;A1;IL8;p1;WARN;4.8\n"
        "S1;SAMPLE;B1;IL8;p1;PASS;4.9\n"
        "S2;SAMPLE;C1;IL8;p1;PASS;5.0\n"
        ";SAMPLE_CONTROL;H1;IL8;p1;PASS;6.0\n"
        ";SAMPLE_CONTROL;H2;IL8;p1;PASS;6.1\n"
    )
    path = tmp_path / "new.csv"
    path.write_text(text, encoding="utf-8")
    df = read_npx(path, verbose=False)
    assert any("2 control measurements" in line for line in df.attrs["read_report"])
    rep = qc_report(df, verbose=False)
    text_rep = str(rep)
    assert "2 control measurements" in text_rep
    assert "3 sample measurements of 2 sample IDs" in text_rep          # controls left out
    assert "(repeats): S1" in text_rep                                    # recognised by WellID
    clean = apply_qc(df, drop_warnings=False, verbose=False)
    assert clean.loc[clean["SampleID"] == "S1", "WellID"].tolist() == ["B1"]   # passed one kept


def test_ids_can_be_hidden_in_reports(npx):
    rep = qc_report(npx, show_ids=False, verbose=False)
    assert "A4" not in str(rep) and "IDs hidden" in str(rep)
    df = apply_qc(npx, show_ids=False, verbose=False)
    assert "A3" not in " ".join(df.attrs["qc_applied"])


def test_samples_without_group_are_not_a_group(npx, sample_list):
    df = add_sample_info(npx, sample_list[sample_list["FIF ID"] != "B2"], on="FIF ID", verbose=False)
    rep = qc_report(df, group="group", verbose=False)
    assert "nan" not in str(rep)
    assert rep.group_by_plate.index.notna().all()
