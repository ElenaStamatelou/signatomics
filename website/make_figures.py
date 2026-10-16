"""Make the charts for the documentation website: PAD000040 (public) and the built-in example data.

Run from the repository root, after downloading PAD000040 into data/PAD000040/ (see
website/docs/guides/example-analysis.md):

    python website/make_figures.py

For every chart it writes two files to website/docs/assets/figures/:
  <name>.json  the chart as data (values rounded): drawn interactively by assets/javascripts/charts.js
  <name>.svg   a static picture: shown until the chart is drawn, and if the chart library cannot load
Every chart is made by the package itself; sample IDs are hidden (the charts' default).
Needs kaleido for the static pictures (pip install kaleido).
"""

import json
from pathlib import Path

import pandas as pd

import signatomics.olink as so

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "PAD000040"
OUT = ROOT / "website" / "docs" / "assets" / "figures"


def _round(x, digits=3):
    if isinstance(x, float):
        return round(x, digits)
    if isinstance(x, list):
        return [_round(v, digits) for v in x]
    if isinstance(x, dict):
        return {k: _round(v, digits) for k, v in x.items()}
    return x


def save(name, fig):
    fig.write_image(OUT / f"{name}.svg")
    spec = _round(json.loads(fig.to_json()))
    spec.get("layout", {}).pop("template", None)          # the charts set their own style
    (OUT / f"{name}.json").write_text(json.dumps(spec, separators=(",", ":")), encoding="utf-8")
    kb = (OUT / f"{name}.json").stat().st_size / 1024
    print(f"{name}: {kb:.0f} KB data")


def main():
    OUT.mkdir(parents=True, exist_ok=True)

    # --- the built-in example (made-up): used in Start and the guides ---
    plate1 = so.read_npx(so.example_file(1), verbose=False)
    both = pd.concat([plate1, so.read_npx(so.example_file(2), verbose=False)], ignore_index=True)
    both_clean = so.apply_qc(both, verbose=False)
    save("example_sample_qc", so.plot_sample_qc(plate1))
    save("example_npx_strip", so.plot_npx_strip(plate1))
    save("example_npx_strip_groups", so.plot_npx_strip(plate1, group="group", samples=so.example_samples(),
                                                       on="sample_id"))
    save("example_bridge_scatter", so.plot_bridge_scatter(both_clean, ["IL6", "CCL4", "HGF"]))
    save("example_bridge_pca", so.plot_bridge_pca(both_clean, so.bridge_plates(both_clean, verbose=False)))

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
    save("bridge_scatter", so.plot_bridge_scatter(clean, scatter_proteins))
    save("bridge_pca", so.plot_bridge_pca(clean, bridged))
    print("scatter proteins:", scatter_proteins)


if __name__ == "__main__":
    main()
