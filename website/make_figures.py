"""Make the charts for the documentation website: PAD000040 (public) and the built-in example data.

Run from the repository root, after downloading PAD000040 into data/PAD000040/ (see
website/docs/guides/example-analysis.md):

    python website/make_figures.py

For every chart it writes two files to website/docs/assets/figures/:
  <name>.json  the chart as data (values rounded): drawn interactively by assets/javascripts/charts.js
  <name>.svg   a static picture: shown until the chart is drawn, and if the chart library cannot load
Every chart is made by the package itself; sample IDs are hidden (the charts' default).
Needs kaleido for the static pictures (pip install kaleido). Add --example-only to make only the
charts of the built-in example (a few minutes; PAD000040 not needed).
"""

import base64
import json
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

import signatomics.olink as so
from signatomics.stats import compare_groups, correlate, plot_volcano, plot_volcano_pair
from signatomics.signature import (find_signature, learning_curve, plot_contributions, plot_curve,
                                   plot_selection_frequency, plot_signature_performance,
                                   protein_contributions, selection_curve)

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "PAD000040"
OUT = ROOT / "website" / "docs" / "assets" / "figures"
DOWNLOADS = ROOT / "website" / "docs" / "assets" / "data"


# the lessons' quality-control decision: P07 is broken; P33 only on Plate 2
EXAMPLE_DROPS = {"P07": "Plate 1", "P33": "Plate 2"}


def _round(x, digits=3):
    if isinstance(x, float):
        return round(x, digits)
    if isinstance(x, list):
        return [_round(v, digits) for v in x]
    if isinstance(x, dict) and "bdata" in x and "dtype" in x:
        # Plotly packs numbers as base64 (8 bytes each): unpack to short rounded numbers
        values = np.frombuffer(base64.b64decode(x["bdata"]), dtype=x["dtype"])
        if "shape" in x:
            values = values.reshape([int(s) for s in str(x["shape"]).split(",")])
        return _round(values.tolist(), digits)
    if isinstance(x, dict):
        return {k: _round(v, digits) for k, v in x.items()}
    return x


# kaleido (the static pictures) runs a hidden browser that can hang when it closes. Each picture
# is therefore made in its own process with a time limit; a hung one is killed with its browser
# and tried again, so one bad picture never stops the run.
_RENDER = "import sys, plotly.io as pio; pio.read_json(sys.argv[1]).write_image(sys.argv[2])"
FAILED = []


def _picture(fig, path, tries=3, limit=120):
    with tempfile.TemporaryDirectory() as tmp:
        spec = Path(tmp) / "figure.json"
        fig.write_json(spec)                      # the full figure, style included
        for attempt in range(1, tries + 1):
            proc = subprocess.Popen([sys.executable, "-c", _RENDER, str(spec), str(path)])
            try:
                if proc.wait(timeout=limit) == 0:
                    return
            except subprocess.TimeoutExpired:
                # kill the process and the browser it started (Windows: the whole tree)
                if os.name == "nt":
                    subprocess.run(["taskkill", "/T", "/F", "/PID", str(proc.pid)], capture_output=True)
                else:
                    proc.kill()
            print(f"  picture {path.name}: try {attempt} of {tries} failed")
    FAILED.append(path.name)


def save(name, fig, static=True):
    if static:
        _picture(fig, OUT / f"{name}.svg")
    spec = _round(json.loads(fig.to_json()), digits=2)
    spec.get("layout", {}).pop("template", None)          # the charts set their own style
    (OUT / f"{name}.json").write_text(json.dumps(spec, separators=(",", ":")), encoding="utf-8")
    kb = (OUT / f"{name}.json").stat().st_size / 1024
    print(f"{name}: {kb:.0f} KB data")


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    write_downloads()

    # --- the built-in example (made-up): used in Start and the guides ---
    plate1 = so.read_npx(so.example_file(1), verbose=False)
    both = pd.concat([plate1, so.read_npx(so.example_file(2), verbose=False)], ignore_index=True)
    both_clean = so.apply_qc(both, verbose=False)
    save("example_sample_qc", so.plot_sample_qc(plate1))
    save("example_npx_strip", so.plot_npx_strip(plate1))
    save("example_npx_strip_groups", so.plot_npx_strip(plate1, group="group", samples=so.example_samples(),
                                                       on="sample_id"))
    save("example_bridge_bland_altman", so.plot_bridge_bland_altman(both_clean))
    save("example_bridge_scatter", so.plot_bridge_scatter(both_clean, ["IL6", "CCL4", "HGF"]))
    save("example_bridge_pca", so.plot_bridge_pca(both_clean, so.bridge_plates(both_clean, verbose=False)))

    # statistical analysis on the analysis table: the lessons remove P07 and P33 (Plate 2 only)
    # in quality control; the "outliers" chart shows what happens when they stay
    files = [so.example_file(1), so.example_file(2)]
    table = so.to_wide(so.prepare(files, so.example_samples(), on="sample_id", verbose=False),
                       verbose=False)
    clean_table = so.to_wide(so.prepare(files, so.example_samples(), on="sample_id",
                                        drop_samples=EXAMPLE_DROPS, verbose=False), verbose=False)
    save("example_volcano_outliers", plot_volcano(compare_groups(table, "group", reference="control",
                                                                 verbose=False)))
    unadjusted = compare_groups(clean_table, "group", reference="control", verbose=False)
    save("example_volcano", plot_volcano(unadjusted))
    save("example_volcano_pair", plot_volcano_pair(
        unadjusted, compare_groups(clean_table, "group", reference="control",
                                   covariates=["age", "sex"], verbose=False)))
    save("example_volcano_age", plot_volcano(correlate(clean_table, "age", verbose=False)))
    signature_figures(clean_table, "group", "case", prefix="example_")
    if "--example-only" in sys.argv:
        return

    # --- PAD000040 (public, CC0): the full example and the reference ---
    plates = pd.concat([so.read_npx(DATA / "inflam_NASH_dt1.raw.csv", verbose=False),
                        so.read_npx(DATA / "inflam_PRA_dt2.raw.csv", verbose=False)],
                       ignore_index=True)
    manifest = dict(samples=DATA / "2025.05.08_sample_manifest.xlsx", on="FIF ID")
    report = so.qc_report(plates, group="diabetes", **manifest, verbose=False)
    clean = so.apply_qc(plates, verbose=False)
    check = so.bridge_check(clean, verbose=False)
    bridged = so.bridge_plates(clean, factors=check.factors, verbose=False)
    data = so.add_sample_info(bridged, manifest["samples"], on=manifest["on"], verbose=False)

    # three proteins whose plate shift differs: the largest up, about none, the largest down
    f = check.factors.dropna(subset=["adjustment"]).sort_values("adjustment")
    near_zero = f.iloc[(f["adjustment"].abs()).argsort()[:1]]
    scatter_proteins = [f["Assay"].iloc[0], near_zero["Assay"].iloc[0], f["Assay"].iloc[-1]]
    # the shift chart for all 92 proteins is very tall: show 24, spread evenly over the range
    pick = f.iloc[[round(i * (len(f) - 1) / 23) for i in range(24)]]["OlinkID"]

    save("group_by_plate", so.plot_group_by_plate(plates, "diabetes", **manifest))
    save("sample_qc", so.plot_sample_qc(plates))
    save("below_lod_by_group", so.plot_below_lod(plates, group="diabetes", **manifest))
    save("npx_strip", so.plot_npx_strip(data, group="diabetes",
                                        proteins=list(report.below_lod_by_group.index[:6])))
    save("bridge_shift", so.plot_bridge_shift(clean[clean["OlinkID"].isin(set(pick))]))
    save("bridge_bland_altman", so.plot_bridge_bland_altman(clean))
    save("bridge_scatter", so.plot_bridge_scatter(clean, scatter_proteins))
    save("bridge_pca", so.plot_bridge_pca(clean, bridged))
    print("scatter proteins:", scatter_proteins)
    stats_figures(plates, manifest)


def stats_figures(plates, manifest):
    """The statistics charts of the PAD000040 example (sample 249 removed, plates bridged)."""
    data = so.add_sample_info(plates, manifest["samples"], on=manifest["on"], verbose=False)
    clean = so.apply_qc(data, drop_samples={"249": "Inflammation PRAESIIDIUM Giraudi 04_03_2025"},
                        show_ids=False, verbose=False)
    table = so.to_wide(so.bridge_plates(clean, verbose=False), verbose=False)
    order = ["No diabetes", "Prediabetes", "T2DM"]
    trend = compare_groups(table, "diabetes", reference="No diabetes", ordered=True, order=order,
                           verbose=False)
    save("volcano_trend", plot_volcano(trend))
    two = table[table["diabetes"].isin(["No diabetes", "T2DM"])]
    save("volcano_pair", plot_volcano_pair(
        compare_groups(two, "diabetes", reference="No diabetes", verbose=False),
        compare_groups(two, "diabetes", reference="No diabetes", covariates=["sex"], verbose=False)))
    signature_figures(two, "diabetes", "T2DM", prefix="", covariates=["sex"], checks=False)


def signature_figures(table, outcome, positive, *, prefix, covariates=None, checks=True):
    """The machine learning charts: the signature and (for the lessons) its checks."""
    result = find_signature(table, outcome, positive=positive, covariates=covariates, verbose=False)
    save(f"{prefix}signature_performance", plot_signature_performance(result))
    save(f"{prefix}selection_frequency", plot_selection_frequency(result))
    if checks and result.signature:
        save(f"{prefix}contributions", plot_contributions(protein_contributions(result, verbose=False)))
        save(f"{prefix}learning_curve", plot_curve(learning_curve(result, verbose=False)))
        save(f"{prefix}selection_curve", plot_curve(selection_curve(result, verbose=False)))


def write_downloads():
    """The synthetic example data as files to download from the website (Learn > Example data)."""
    DOWNLOADS.mkdir(parents=True, exist_ok=True)
    for plate in (1, 2):
        shutil.copyfile(so.example_file(plate), DOWNLOADS / f"signatomics_example_plate{plate}.csv")
    so.example_samples().to_excel(DOWNLOADS / "signatomics_example_samples.xlsx", index=False,
                                  sheet_name="sample metadata")
    table = so.to_wide(so.prepare([so.example_file(1), so.example_file(2)], so.example_samples(),
                                  on="sample_id", drop_samples=EXAMPLE_DROPS, verbose=False),
                       verbose=False)
    table.round(4).to_excel(DOWNLOADS / "signatomics_example_analysis_table.xlsx", index=False,
                            sheet_name="analysis table")
    print("downloads:", ", ".join(sorted(f.name for f in DOWNLOADS.iterdir())))


if __name__ == "__main__":
    main()
    if FAILED:
        print("pictures not made (the chart data is there; rerun to retry):", ", ".join(FAILED))
