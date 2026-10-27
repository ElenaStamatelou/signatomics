"""Tests for the analysis table (to_wide) and signatomics.stats (compare_groups, correlate)."""

import numpy as np
import pandas as pd
import pytest
from scipy import stats

from signatomics.olink import example_file, example_samples, prepare, to_wide
from signatomics.stats import compare_groups, correlate, fdr_bh, holm
from signatomics.stats.compare import _jonckheere, _welch_anova, _welch_two


@pytest.fixture(scope="module")
def prepared():
    return prepare([example_file(1), example_file(2)], example_samples(), on="sample_id",
                   verbose=False)


@pytest.fixture(scope="module")
def table(prepared):
    return to_wide(prepared, verbose=False)


def _toy(groups: dict, proteins: dict) -> pd.DataFrame:
    t = pd.DataFrame({"group": [g for g, n in groups.items() for _ in range(n)]})
    for name, values in proteins.items():
        t[name] = values
    t.attrs["proteins"] = list(proteins)
    return t


# --- to_wide --------------------------------------------------------------------------------

def test_to_wide_one_row_per_sample(prepared, table):
    assert table.shape[0] == 59 and len(table.attrs["proteins"]) == 92
    assert list(table.columns[:4]) == ["sample_id", "group", "age", "sex"]
    assert "PlateID" not in table.columns and table.index.name == "sample"


def test_to_wide_averages_bridge_samples(prepared, table):
    p34 = prepared[(prepared["sample"] == "P34") & (prepared["Assay"] == "IL6")]
    assert len(p34) == 2
    assert table.loc["P34", "IL6"] == pytest.approx(p34["NPX"].mean())


def test_to_wide_below_lod(prepared):
    b = to_wide(prepared, values="below_lod", verbose=False)
    assert b["IL4"].mean() > 0.75 and b["IL6"].mean() < 0.2
    with pytest.raises(ValueError):
        to_wide(prepared, values="LOD")


def test_to_wide_names_shared_assay_by_olinkid():
    rows = [{"SampleID": s, "OlinkID": o, "Assay": "IL6", "Panel": p, "NPX": v}
            for s, v in (("A", 1.0), ("B", 2.0)) for o, p in (("OID1", "Inflammation"), ("OID2", "Oncology"))]
    t = to_wide(pd.DataFrame(rows), verbose=False)
    assert t.attrs["proteins"] == ["IL6 (OID1)", "IL6 (OID2)"]


# --- the tests themselves ---------------------------------------------------------------------

def test_welch_matches_scipy():
    rng = np.random.default_rng(0)
    a, b = rng.normal(0, 1, 20), rng.normal(0.5, 3, 60)
    d, lo, hi, p = _welch_two(a, b, 0.95)
    ref = stats.ttest_ind(b, a, equal_var=False)
    ci = ref.confidence_interval(0.95)
    assert np.allclose([d, lo, hi, p], [b.mean() - a.mean(), ci.low, ci.high, ref.pvalue])
    f, p_anova = _welch_anova([a, b])            # two groups: Welch ANOVA = Welch t-test
    assert p_anova == pytest.approx(ref.pvalue) and f == pytest.approx(ref.statistic ** 2)


def test_jonckheere_finds_a_trend_and_not_its_absence():
    rng = np.random.default_rng(1)
    up = [rng.normal(m, 1, 30) for m in (0, 0.7, 1.4)]
    assert _jonckheere(up)[1] < 0.001 and _jonckheere(up)[0] > 0
    assert _jonckheere(up[::-1])[0] < 0


def test_fdr_and_holm_skip_missing_p_values():
    p = [0.01, 0.04, np.nan, 0.03, 0.2, 0.001]
    assert np.allclose(fdr_bh(p), [0.025, 0.05, np.nan, 0.05, 0.2, 0.005], equal_nan=True)
    assert np.allclose(holm(p), [0.04, 0.09, np.nan, 0.09, 0.2, 0.005], equal_nan=True)


# --- compare_groups ---------------------------------------------------------------------------

def test_compare_groups_finds_the_built_in_differences(table):
    r = compare_groups(table, "group", reference="control", verbose=False)
    top = set(r.table.loc[r.table["fdr"] < 0.05, "protein"])
    assert top <= {"IL6", "CXCL10", "HGF", "OSM", "IL-17C"} and {"IL6", "OSM", "IL-17C"} <= top
    il6 = r.table.set_index("protein").loc["IL6"]
    assert il6["difference"] == pytest.approx(il6["mean_case"] - il6["mean_control"])
    assert il6["ci_low"] < il6["difference"] < il6["ci_high"]


def test_reference_sets_the_direction(table):
    a = compare_groups(table, "group", reference="control", verbose=False).table.set_index("protein")
    b = compare_groups(table, "group", reference="case", verbose=False).table.set_index("protein")
    assert a.loc["IL6", "difference"] == pytest.approx(-b.loc["IL6", "difference"])
    assert a.loc["IL6", "p"] == pytest.approx(b.loc["IL6", "p"])


def test_mann_whitney_and_below_lod(prepared, table):
    b = to_wide(prepared, values="below_lod", verbose=False)
    r = compare_groups(table, "group", reference="control", test="mann-whitney", below_lod=b,
                       verbose=False)
    t = r.table.set_index("protein")
    assert {"median_case", "below_lod_case", "below_lod_control"} <= set(t.columns)
    assert t.loc["IL-17C", "below_lod_control"] > t.loc["IL-17C", "below_lod_case"]


def test_untestable_proteins_are_reported_and_left_out_of_the_fdr():
    t = _toy({"a": 5, "b": 5}, {"flat": [1.0] * 10, "few": [1.0, np.nan, np.nan, np.nan, np.nan] + [2.0] * 5,
                                 "real": list(range(10))})
    r = compare_groups(t, "group", reference="a", verbose=False)
    assert list(r.table["protein"]) == ["real"] and not r.table["fdr"].isna().any()
    assert dict(zip(r.not_tested["protein"], r.not_tested["reason"])) == {
        "flat": "no variation in either group", "few": "fewer than 2 values in a"}


def test_samples_without_group_are_dropped_and_counted():
    t = _toy({"a": 4, "b": 4}, {"x": [1, 2, 3, 4, 5, 6, 7, 8]})
    t.loc[0, "group"] = None
    r = compare_groups(t, "group", reference="a", verbose=False)
    assert r.table.loc[0, "n_a"] == 3
    assert any("1 sample without a value" in line for line in r.lines)


def test_three_groups_omnibus_pairwise_and_trend():
    rng = np.random.default_rng(2)
    g = {"low": 20, "mid": 20, "high": 20}
    t = _toy(g, {"up": np.concatenate([rng.normal(m, 1, 20) for m in (0, 1, 2)]),
                 "none": rng.normal(0, 1, 60)})
    r = compare_groups(t, "group", reference="low", verbose=False)
    assert r.table.iloc[0]["protein"] == "up" and r.table.iloc[0]["p"] < 1e-6
    pw = r.pairwise
    assert len(pw) == 6
    assert set(pw.loc[(pw["protein"] == "up") & (pw["reference"] == "low"), "group"]) == {"mid", "high"}
    assert pw.loc[pw["protein"] == "up", "omnibus_significant"].all()
    assert (pw["p_holm"] >= pw["p"]).all()
    tr = compare_groups(t, "group", reference="low", ordered=True, order=["low", "mid", "high"],
                        verbose=False)
    up = tr.table.set_index("protein").loc["up"]
    assert up["ci_low"] < 1 < up["ci_high"]     # true slope: 1 NPX per step


def test_input_errors():
    t = _toy({"a": 3, "b": 1}, {"x": [1.0, 2, 3, 4]})
    with pytest.raises(ValueError, match="fewer than 2 samples"):
        compare_groups(t, "group", verbose=False)
    with pytest.raises(ValueError, match="order"):
        compare_groups(t, "group", ordered=True, verbose=False)
    with pytest.raises(ValueError, match="test must be"):
        compare_groups(t, "group", test="student", verbose=False)
    plain = pd.DataFrame({"group": ["a", "b"], "x": [1, 2]})
    with pytest.raises(ValueError, match="proteins"):
        compare_groups(plain, "group", verbose=False)


# --- correlate --------------------------------------------------------------------------------

def test_correlate_matches_scipy_and_reports_intervals():
    rng = np.random.default_rng(3)
    x = rng.normal(0, 1, 50)
    t = pd.DataFrame({"hba1c": x, "p1": x + rng.normal(0, 1, 50), "p2": rng.normal(0, 1, 50)})
    t.attrs["proteins"] = ["p1", "p2"]
    for method, fn in (("spearman", stats.spearmanr), ("pearson", stats.pearsonr)):
        r = correlate(t, "hba1c", method=method, verbose=False).table.set_index("protein")
        ref = fn(t["hba1c"], t["p1"])
        assert r.loc["p1", "r"] == pytest.approx(ref.statistic)
        assert r.loc["p1", "p"] == pytest.approx(ref.pvalue)
        assert r.loc["p1", "ci_low"] < r.loc["p1", "r"] < r.loc["p1", "ci_high"]


def test_correlate_on_the_analysis_table(table):
    r = correlate(table, "age", verbose=False)
    assert len(r.table) == 92 and r.table["n"].eq(59).all()


def test_volcano_for_two_groups_and_correlation(table):
    import plotly.graph_objects as go

    from signatomics.stats import plot_volcano
    fig = plot_volcano(compare_groups(table, "group", reference="control", verbose=False))
    assert isinstance(fig, go.Figure)
    assert {"IL-17C", "OSM", "IL6"} <= {a.text for a in fig.layout.annotations}
    assert isinstance(plot_volcano(correlate(table, "age", verbose=False)), go.Figure)
    with pytest.raises(ValueError, match="effect"):
        three = _toy({"a": 5, "b": 5, "c": 5}, {"x": list(range(15))})
        plot_volcano(compare_groups(three, "group", reference="a", verbose=False))


def test_the_reference_is_always_required():
    t = _toy({"control": 4, "case": 4}, {"x": [1.0, 2, 3, 4, 6, 7, 8, 9]})
    with pytest.raises(ValueError, match="reference='case' or reference='control'"):
        compare_groups(t, "group", verbose=False)
    with pytest.raises(ValueError, match="Which group is the reference"):
        compare_groups(t, "group", order=["control", "case"], verbose=False)
    assert compare_groups(t, "group", reference="control", verbose=False).table.loc[0, "difference"] == 5.0
    three = _toy({"low": 4, "mid": 4, "high": 4}, {"x": [1.0, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]})
    with pytest.raises(ValueError, match="reference='high', reference='low' or reference='mid'"):
        compare_groups(three, "group", verbose=False)
    with pytest.raises(ValueError, match="Which group is the reference"):
        compare_groups(three, "group", ordered=True, order=["low", "mid", "high"], verbose=False)


def test_pairs_subtract_the_reference_then_the_earlier_group():
    t = _toy({"low": 4, "mid": 4, "high": 4}, {"x": [1.0, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]})
    pw = compare_groups(t, "group", reference="mid", ordered=True, order=["low", "mid", "high"],
                        verbose=False).pairwise
    pairs = set(zip(pw["group"], pw["reference"]))
    assert pairs == {("low", "mid"), ("high", "mid"), ("high", "low")}
    assert pw.set_index("group").loc["low", "difference"] == -4.0
