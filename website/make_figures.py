"""Make the example figures for the documentation website, from the public dataset PAD000040.

Run from the repository root, after downloading PAD000040 into data/PAD000040/ (see
website/docs/getting-started.md):

    python website/make_figures.py

Writes SVG files to website/docs/assets/figures/. Every chart is made by the package itself,
on the same steps as the Get started page; sample IDs are hidden (the charts' default).
Needs kaleido for the static export (pip install kaleido).
"""

from pathlib import Path

import pandas as pd

import signatomics.olink as so

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "PAD000040"
OUT = ROOT / "website" / "docs" / "assets" / "figures"


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    plates = pd.concat([so.read_npx(DATA / "inflam_NASH_dt1.raw.csv", verbose=False),
                        so.read_npx(DATA / "inflam_PRA_dt2.raw.csv", verbose=False)],
                       ignore_index=True)
    data = so.add_sample_info(plates, DATA / "2025.05.08_sample_manifest.xlsx", on="FIF ID",
                              verbose=False)
    report = so.qc_report(data, group="diabetes", verbose=False)
    clean = so.apply_qc(data, verbose=False)
    check = so.bridge_check(clean, verbose=False)
    bridged = so.bridge_plates(clean, factors=check.factors, verbose=False)

    # three proteins whose plate shift differs: the largest up, about none, the largest down
    f = check.factors.dropna(subset=["adjustment"]).sort_values("adjustment")
    near_zero = f.iloc[(f["adjustment"].abs()).argsort()[:1]]
    scatter_proteins = [f["Assay"].iloc[0], near_zero["Assay"].iloc[0], f["Assay"].iloc[-1]]

    # proteins whose detectability differs most between the diabetes groups, for the strip plot
    strip_proteins = list(report.below_lod_by_group.index[:6])

    # the shift chart for all 92 proteins is very tall: show 24, spread evenly over the range
    pick = f.iloc[[round(i * (len(f) - 1) / 23) for i in range(24)]]["OlinkID"]
    shift_subset = clean[clean["OlinkID"].isin(set(pick))]

    figures = {
        "group_by_plate": so.plot_group_by_plate(data, "diabetes"),
        "sample_qc": so.plot_sample_qc(data),
        "below_lod_by_group": so.plot_below_lod(data, group="diabetes"),
        "npx_strip": so.plot_npx_strip(bridged, group="diabetes", proteins=strip_proteins),
        "bridge_shift": so.plot_bridge_shift(shift_subset),
        "bridge_scatter": so.plot_bridge_scatter(clean, scatter_proteins),
        "bridge_pca": so.plot_bridge_pca(clean, bridged),
    }
    for name, fig in figures.items():
        fig.write_image(OUT / f"{name}.svg")
        print(f"{name}.svg  ({fig.layout.width} x {fig.layout.height})")
    print("scatter proteins:", scatter_proteins)
    print("strip proteins:", strip_proteins)


if __name__ == "__main__":
    main()
