"""Tests for the made-up example data and for prepare (the usual preprocessing in one call)."""

import pandas as pd
import pytest

from signatomics.olink import (add_sample_info, apply_qc, bridge_plates, example_data,
                               example_file, example_samples, prepare, qc_report, read_npx)


def test_example_data_is_the_same_every_time():
    a, b = example_data(), example_data()
    pd.testing.assert_frame_equal(a, b)
    assert a["SampleID"].nunique() == 60 and a["Assay"].nunique() == 24
    assert set(a["PlateID"]) == {"Plate 1", "Plate 2"}


def test_example_data_has_what_the_checks_should_find():
    rep = qc_report(example_data(), verbose=False)
    m = rep.measurements
    assert m["warning"].sum() == 2 and m["fail"].sum() == 1
    assert not m["repeat"].any()
    assert m.loc[m["bridging"], "SampleID"].nunique() == 12
    assert sorted(rep.samples.loc[rep.samples["outlier"], "SampleID"]) == ["P07", "P33"]
    assert (rep.proteins["below_lod"] > 0.5).sum() == 3


def test_example_files_are_the_two_plates():
    one = read_npx(example_file(), verbose=False)
    two = read_npx(example_file(2), verbose=False)
    assert set(one["PlateID"]) == {"Plate 1"} and set(two["PlateID"]) == {"Plate 2"}
    both = pd.concat([one, two], ignore_index=True)
    pd.testing.assert_frame_equal(both, example_data(), check_dtype=False, check_like=True)
    assert one.attrs["read_report"][-1] == "result: 960 rows = 40 samples x 24 assays"
    with pytest.raises(ValueError):
        example_file(3)


def test_example_samples_match_every_sample():
    data = add_sample_info(example_data(), example_samples(), on="sample_id", verbose=False)
    assert data["group"].notna().all()
    assert set(data["group"]) == {"case", "control"}


def test_group_checks_look_up_the_sample_data_without_adding_it():
    data = example_data()
    rep = qc_report(data, group="group", samples=example_samples(), on="sample_id", verbose=False)
    assert rep.below_lod_by_group is not None and rep.group_by_plate is not None
    assert "group" not in data.columns


def test_prepare_equals_the_steps_with_their_defaults():
    data, samples = example_data(), example_samples()
    one_call = prepare(data, samples, on="sample_id", verbose=False)
    steps = add_sample_info(bridge_plates(apply_qc(data, verbose=False), verbose=False),
                            samples, on="sample_id", verbose=False)
    pd.testing.assert_frame_equal(one_call.reset_index(drop=True), steps.reset_index(drop=True))


def test_prepare_report_is_short_and_lists_what_to_look_at():
    out = prepare(example_data(), verbose=False)
    report = out.attrs["prepare_report"]
    assert report[0] == "data: 60 samples, 24 assays, 2 plates"
    assert any(line.startswith("bridging: Plate 2 put on the scale of Plate 1") for line in report)
    assert any(line.strip().startswith("outliers") and "P33" in line for line in report)
    assert len(out.attrs["prepare_details"]) > len(report)
    hidden = " ".join(prepare(example_data(), show_ids=False, verbose=False).attrs["prepare_report"])
    assert "P33" not in hidden and "P20" not in hidden


def test_prepare_reads_files():
    out = prepare([example_file(1), example_file(2)], verbose=False)
    assert out.attrs["prepare_report"][0].startswith("read 2 files")
    assert "Bridge_adjustment" in out.columns


def test_prepare_one_plate_is_not_bridged():
    out = prepare(example_file(), verbose=False)
    assert "Bridge_adjustment" not in out.columns
    assert not any(line.startswith("bridging") for line in out.attrs["prepare_report"])


def test_prepare_needs_on_with_samples():
    with pytest.raises(ValueError, match="on="):
        prepare(example_data(), example_samples())
