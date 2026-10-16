"""Tests for the quality charts: they build, show what they should, and hide IDs by default."""

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import pytest

from signatomics.olink import (plot_below_lod, plot_group_by_plate, plot_npx_strip,
                               plot_sample_qc, qc_report)


@pytest.fixture
def npx():
    rng = np.random.default_rng(1)
    rows = []
    for k in range(40):
        sid, plate = f"S{k:02d}", "plate1" if k < 20 else "plate2"
        group = ["no", "pre", "t2d"][k % 3] if plate == "plate2" else ["no", "t2d"][k % 2]
        for p, (mean, lod) in {"P1": (5, 0.0), "P2": (3, 0.5), "LOW": (0, 1.0)}.items():
            rows.append(dict(SampleID=sid, Index=k, Assay=p, OlinkID=f"OID_{p}", PlateID=plate,
                             Panel="Panel A", LOD=lod, NPX=mean + rng.normal(0, 0.5), group=group))
    df = pd.DataFrame(rows)
    df.loc[df["SampleID"] == "S07", "NPX"] += 6          # one clearly abnormal sample
    return df


def test_below_lod_chart_marks_nothing_by_default(npx):
    fig = plot_below_lod(npx)
    assert isinstance(fig, go.Figure)
    assert not any("removed" in (t.name or "") for t in fig.data)


def test_below_lod_chart_with_threshold(npx):
    fig = plot_below_lod(npx, max_below_lod=0.5)
    removed = [t for t in fig.data if "removed" in t.name][0]
    assert "LOW" in list(removed.y)                         # mostly below LOD -> red bar


def test_below_lod_chart_per_group(npx):
    fig = plot_below_lod(npx, group="group")
    assert {t.name for t in fig.data} == {"no", "pre", "t2d"}


def test_strip_chart_marks_below_lod_hollow_and_hides_ids(npx):
    fig = plot_npx_strip(npx, group="group")
    names = [t.name for t in fig.data]
    assert any("(below LOD)" in n for n in names)
    assert "detection limit (LOD)" in names
    hollow = [t for t in fig.data if t.name and "(below LOD)" in t.name][0]
    assert hollow.marker.color == "rgba(0,0,0,0)"
    assert not any("sample:" in (t.hovertemplate or "") for t in fig.data)


def test_strip_chart_can_show_ids(npx):
    fig = plot_npx_strip(npx, group="group", show_ids=True)
    assert any("sample:" in (t.hovertemplate or "") for t in fig.data)


def test_strip_chart_refuses_too_many_groups(npx):
    npx = npx.assign(group=npx["SampleID"])
    with pytest.raises(ValueError, match="groups"):
        plot_npx_strip(npx, group="group")


def test_group_by_plate_chart(npx):
    fig = plot_group_by_plate(npx, "group")
    assert {t.name for t in fig.data} == {"plate1", "plate2"}


def test_sample_qc_finds_the_abnormal_sample():
    # realistic size: 60 samples x 20 proteins; one sample shifted on every protein
    rng = np.random.default_rng(2)
    rows = [dict(SampleID=f"S{k:02d}", Index=k, Assay=f"P{j}", PlateID="plate1",
                 NPX=j * 0.5 + rng.normal(0, 0.4) + (4 if k == 7 else 0))
            for k in range(60) for j in range(20)]
    df = pd.DataFrame(rows)
    rep = qc_report(df, outlier_sd=3, verbose=False)
    table = rep.samples
    assert table.loc[table["outlier"], "SampleID"].tolist() == ["S07"]
    assert "outliers: 1 of 60 samples (S07)" in str(rep)
    fig = plot_sample_qc(df, sd=3)
    assert any(t.name and "outlier" in t.name for t in fig.data)


def test_bland_altman_shows_the_plate_shift_and_its_removal():
    import plotly.graph_objects as go

    from signatomics.olink import apply_qc, bridge_plates, example_data, plot_bridge_bland_altman
    clean = apply_qc(example_data(), verbose=False)
    before = plot_bridge_bland_altman(clean)
    after = plot_bridge_bland_altman(bridge_plates(clean, verbose=False))
    assert isinstance(before, go.Figure) and len(before.data) == 1

    def bias(fig):
        return float([a.text for a in fig.layout.annotations if a.text.startswith("median")][0].split()[1])
    assert bias(before) > 0.3 and abs(bias(after)) < 0.1
    assert len(plot_bridge_bland_altman(clean, ["IL6", "HGF"]).data) == 2
