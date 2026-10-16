"""Tests for bridging (bridge_check with its differences/factors tables, bridge_plates), the
bridging charts, and apply_qc(drop_samples=...), on made-up data with a known plate shift.

The made-up data: 3 proteins measured on two plates; plate 2 reads every protein higher by a
known amount (P1 +1.0, P2 -0.5, P3 +2.0). 12 bridging samples are on both plates, one of them
with a quality warning on plate 2 and wildly wrong values there (it must not be used), and one
written differently on the two plates ("B 01" vs "B_01").
"""

import numpy as np
import pandas as pd
import pytest

from signatomics.olink import (apply_qc, bridge_check, bridge_plates, plot_bridge_pca,
                               plot_bridge_scatter, plot_bridge_shift)

SHIFT = {"P1": 1.0, "P2": -0.5, "P3": 2.0}       # plate 2 minus plate 1
LOD = {"P1": 1.0, "P2": 0.5, "P3": 2.0}


def make_data(seed=1, n_bridges=12, plates=("plate1", "plate2")):
    rng = np.random.default_rng(seed)
    rows = []
    truth = {f"S{i:02d}": rng.normal(5, 1, size=3) for i in range(40)}
    bridges = [f"B_{i:02d}" for i in range(n_bridges)]
    truth |= {b: rng.normal(5, 1, size=3) for b in bridges}

    def add(sid, plate, index, qc="Pass", junk=False):
        for k, p in enumerate(SHIFT):
            shift = 0.0 if plate == plates[0] else SHIFT[p] * (1 + plates.index(plate) - 1)
            value = truth[sid.replace(" ", "_")][k] + shift + rng.normal(0, 0.05)
            rows.append(dict(SampleID=sid, Index=index, OlinkID=f"OID_{p}", Assay=p, PlateID=plate,
                             QC_Warning=qc, NPX=-8.0 if junk else value,
                             LOD=LOD[p] + shift, MaxLOD=LOD[p] + shift + 0.1))

    for i in range(40):
        add(f"S{i:02d}", plates[0] if i < 25 else plates[1], i)
    for j, b in enumerate(bridges):
        add(b, plates[0], 100 + j)
        for pl in plates[1:]:
            if j == 0:
                add("B 00", pl, 200 + j)                     # same sample, written differently
            elif j == 1:
                add(b, pl, 200 + j, qc="Warning", junk=True) # broken measurement, with a warning
            else:
                add(b, pl, 200 + j)
    return pd.DataFrame(rows)


@pytest.fixture
def npx():
    return make_data()


def test_differences_one_row_per_bridge_and_protein(npx):
    d = bridge_check(npx, verbose=False).differences
    assert d.attrs["reference"] == "plate1"                  # the plate with most samples
    assert d["sample"].nunique() == 11                       # 12 bridges minus the warned one
    assert set(d["plate"]) == {"plate2"}
    assert np.allclose(d.groupby("Assay")["difference"].median().reindex(list(SHIFT)),
                       list(SHIFT.values()), atol=0.1)


def test_warned_bridge_is_used_only_when_asked(npx):
    with_warned = bridge_check(npx, exclude_warnings=False, verbose=False).differences
    assert with_warned["sample"].nunique() == 12


def test_factors_recover_the_known_shift(npx):
    f = bridge_check(npx, verbose=False).factors
    got = f.set_index("Assay")["adjustment"].reindex(list(SHIFT))
    assert np.allclose(got, [-v for v in SHIFT.values()], atol=0.1)
    assert (f["bridges"] == 11).all()
    assert not f["few_bridges"].any()
    assert (f["reference"] == "plate1").all()


def test_bridging_lines_up_the_plates(npx):
    b = bridge_plates(npx, verbose=False)
    chk = bridge_check(b, verbose=False)
    assert chk.proteins["shift"].abs().max() < 0.05
    assert "already bridged" in " ".join(chk.lines)


def test_reference_unchanged_and_originals_kept(npx):
    b = bridge_plates(npx, verbose=False)
    ref = b["PlateID"] == "plate1"
    assert np.allclose(b.loc[ref, "NPX"], npx.loc[ref, "NPX"])
    assert (b.loc[ref, "Bridge_adjustment"] == 0).all()
    assert np.allclose(b["NPX_unbridged"], npx["NPX"])
    assert np.allclose(b["LOD_unbridged"], npx["LOD"])


def test_lod_shifted_with_npx_so_below_lod_is_unchanged(npx):
    b = bridge_plates(npx, verbose=False)
    assert ((b["NPX"] < b["LOD"]) == (npx["NPX"] < npx["LOD"])).all()
    assert np.allclose(b["MaxLOD"] - b["LOD"], 0.1)
    # after bridging, plate 2's LOD is back on the plate 1 scale
    assert np.allclose(b.groupby("Assay")["LOD"].max() - b.groupby("Assay")["LOD"].min(), 0, atol=0.1)


def test_chosen_reference_and_report(npx):
    b = bridge_plates(npx, reference="plate2", verbose=False)
    assert (b.loc[b["PlateID"] == "plate2", "Bridge_adjustment"] == 0).all()
    assert "plate2" in b.attrs["bridge_report"][0]
    with pytest.raises(ValueError, match="not found"):
        bridge_plates(npx, reference="plate9", verbose=False)


def test_refuses_to_bridge_twice(npx):
    b = bridge_plates(npx, verbose=False)
    with pytest.raises(ValueError, match="already bridged"):
        bridge_plates(b, verbose=False)


def test_saved_factors_can_be_applied(npx):
    f = bridge_check(npx, verbose=False).factors
    a = bridge_plates(npx, factors=f, verbose=False)
    b = bridge_plates(npx, verbose=False)
    assert np.allclose(a["NPX"], b["NPX"])


def test_only_chosen_bridges(npx):
    f = bridge_check(npx, bridges=["B_02", "B_03", "B_04"], verbose=False).factors
    assert (f["bridges"] == 3).all() and f["few_bridges"].all()
    report = " ".join(bridge_plates(npx, factors=f, verbose=False).attrs["bridge_report"])
    assert "few bridge samples" in report


def test_no_bridges_is_reported_not_crashed():
    df = make_data()
    df = df[~df["SampleID"].str.startswith("B")]
    chk = bridge_check(df, verbose=False)
    assert "no bridging samples" in " ".join(chk.lines)
    b = bridge_plates(df, verbose=False)
    assert "WARNING not adjusted" in " ".join(b.attrs["bridge_report"])
    assert np.allclose(b["NPX"], df["NPX"])                  # nothing shifted without bridges


def test_three_plates():
    df = make_data(plates=("p1", "p2", "p3"))
    b = bridge_plates(df, reference="p1", verbose=False)
    f = bridge_check(df, reference="p1", verbose=False).factors
    assert set(f["plate"]) == {"p2", "p3"}
    assert bridge_check(b, reference="p1", verbose=False).proteins["shift"].abs().max() < 0.05


def test_needs_npx_and_two_plates(npx):
    with pytest.raises(ValueError, match="NPX"):
        bridge_plates(npx.drop(columns="NPX").assign(Ct=1.0), verbose=False)
    with pytest.raises(ValueError, match="nothing to bridge"):
        bridge_plates(npx[npx["PlateID"] == "plate1"], verbose=False)
    with pytest.raises(KeyError, match="Project"):
        bridge_plates(npx, batch="Project", verbose=False)


def test_group_check_and_odd_bridge(npx):
    groups = {sid: ("a" if i % 2 else "b") for i, sid in enumerate(sorted(npx["SampleID"].unique()))}
    df = npx.assign(group=npx["SampleID"].map(groups))
    df.loc[(df["SampleID"] == "B_05") & (df["PlateID"] == "plate2"), "NPX"] -= 3   # odd bridge
    chk = bridge_check(df, group="group", show_ids=True, verbose=False)
    text = " ".join(chk.lines)
    assert "among the bridge samples" in text
    assert chk.by_group is not None and "difference" in chk.by_group.columns
    assert "B_05" in text                                    # flagged as shifting differently
    hidden = " ".join(bridge_check(df, group="group", show_ids=False, verbose=False).lines)
    assert "B_05" not in hidden


# --- charts --------------------------------------------------------------------------------

def test_shift_chart(npx):
    fig = plot_bridge_shift(npx)
    medians = [t for t in fig.data if "median" in t.name]
    assert len(medians) == 1 and len(medians[0].x) == 3
    assert "sample" not in fig.data[0].hovertemplate               # IDs hidden by default


def test_pca_chart_before_after(npx):
    b = bridge_plates(npx, verbose=False)
    fig = plot_bridge_pca(npx, b)
    assert "before bridging" in fig.layout.title.text and "after bridging" in fig.layout.title.text


def test_scatter_chart(npx):
    fig = plot_bridge_scatter(npx, ["P1", "P3"])
    assert len([t for t in fig.data if t.mode == "markers"]) == 2
    with pytest.raises(ValueError, match="Not found"):
        plot_bridge_scatter(npx, ["nope"])


# --- apply_qc(drop_samples=...) ------------------------------------------------------------

def test_drop_samples_everywhere_or_on_one_plate(npx):
    gone = apply_qc(npx, drop_samples=["b-03"], verbose=False)
    assert "B_03" not in set(gone["SampleID"])
    one = apply_qc(npx, drop_samples={"B_03": "plate2"}, verbose=False)
    assert set(one.loc[one["SampleID"] == "B_03", "PlateID"]) == {"plate1"}
    assert "chosen by you" in " ".join(one.attrs["qc_applied"])
    none = apply_qc(npx, drop_samples=["ZZZ"], verbose=False)
    assert "none of the given samples" in " ".join(none.attrs["qc_applied"])
