"""Write the printed output of every example on the website, from the real code.

Each output goes to website/snippets/outputs/<name>.txt, and the pages include it inside their
output box with  --8<-- "outputs/<name>.txt". Run this after a change to the library, then
rebuild the site: the boxes can never show an older version of a report.

Run from the repository root:

    python website/make_outputs.py                  # everything (PAD000040 needs data/PAD000040/)
    python website/make_outputs.py --example-only   # only the built-in synthetic example
"""

from __future__ import annotations

import contextlib
import io
import sys
import warnings
from pathlib import Path

import pandas as pd

import signatomics.olink as so
from signatomics.signature import find_signature, learning_curve, protein_contributions, selection_curve
from signatomics.stats import compare_groups, correlate

warnings.filterwarnings("ignore")
ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "website" / "snippets" / "outputs"
DATA = ROOT / "data" / "PAD000040"
DROPS = {"P07": "Plate 1", "P33": "Plate 2"}          # the lessons' decisions in quality control


def save(name: str, run) -> None:
    """Run ``run()`` and save what it prints, plus the text of what it returns (if not None)."""
    buffer = io.StringIO()
    with contextlib.redirect_stdout(buffer):
        value = run()
    if value is not None:
        buffer.write(str(value) + "\n")
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / f"{name}.txt").write_text(buffer.getvalue().rstrip() + "\n", encoding="utf-8")
    print("  ", name)


def printed(func, *args, **kwargs):
    """Call a function for the report it prints; return nothing."""
    func(*args, **kwargs)


def quiet(func, *args, **kwargs):
    """Call a function without its printed report; return its result."""
    with contextlib.redirect_stdout(io.StringIO()):
        return func(*args, **kwargs)


def example():
    print("synthetic example")
    plate1 = quiet(so.read_npx, so.example_file())
    samples1 = so.example_samples(1)
    save("read_plate1", lambda: printed(so.read_npx, so.example_file()))
    save("qc_report_plate1", lambda: so.qc_report(plate1, verbose=False))
    save("apply_qc_plate1", lambda: printed(so.apply_qc, plate1))
    save("apply_qc_plate1_keep_fail", lambda: printed(so.apply_qc, plate1, drop_fail=False))
    save("apply_qc_plate1_drop_p07", lambda: printed(so.apply_qc, plate1, drop_samples=["P07"]))
    save("example_samples_head", lambda: samples1.head())
    save("add_sample_info_plate1", lambda: printed(so.add_sample_info, plate1, samples1, on="sample_id"))

    # Start: one plate, all the way to a first group comparison (Elena, 17 Nov 2026)
    start_clean = quiet(so.apply_qc, plate1, drop_samples=["P07"])
    save("start_add_sample_info", lambda: printed(so.add_sample_info, start_clean, samples1, on="sample_id"))
    start_table = quiet(so.to_wide, quiet(so.add_sample_info, start_clean, samples1, on="sample_id"))
    save("start_to_wide", lambda: printed(so.to_wide, quiet(so.add_sample_info, start_clean, samples1,
                                                                 on="sample_id")))
    save("start_compare_groups", lambda: printed(compare_groups, start_table, "group", reference="control"))

    both = pd.concat([quiet(so.read_npx, so.example_file(1)), quiet(so.read_npx, so.example_file(2))],
                     ignore_index=True)
    clean = quiet(so.apply_qc, both)
    save("bridge_check", lambda: so.bridge_check(clean, verbose=False))
    save("bridge_plates", lambda: printed(so.bridge_plates, clean))
    clean_drops = quiet(so.apply_qc, both, drop_samples=DROPS)
    combined = quiet(so.bridge_plates, clean_drops)
    save("bridge_plates_drops", lambda: printed(so.bridge_plates, clean_drops))
    save("add_sample_info_combined",
         lambda: printed(so.add_sample_info, combined, so.example_samples(), on="sample_id"))
    with_samples = quiet(so.add_sample_info, combined, so.example_samples(), on="sample_id")
    save("to_wide", lambda: printed(so.to_wide, with_samples))
    table = quiet(so.to_wide, with_samples)

    save("analysis_table_head", lambda: table.iloc[:3, :7].round(2))
    clean_all = quiet(so.apply_qc, both)
    no_drops = quiet(so.add_sample_info, quiet(so.bridge_plates, clean_all), so.example_samples(),
                     on="sample_id")
    save("to_wide_no_drops", lambda: printed(so.to_wide, no_drops))

    save("compare_groups", lambda: printed(compare_groups, table, "group", reference="control"))
    save("correlate", lambda: printed(correlate, table, "age"))
    corr = quiet(correlate, table, "age")
    cols = ["protein", "n", "r", "ci_low", "ci_high", "p", "fdr"]
    save("correlate_table", lambda: corr.table[cols].head(3).to_string(
        index=False, float_format=lambda v: f"{v:.3g}"))

    result = quiet(find_signature, table, "group", positive="case")
    save("find_signature", lambda: result)
    check = quiet(protein_contributions, result)
    save("protein_contributions", lambda: check)
    save("protein_contributions_table",
         lambda: check.table[["protein", "alone", "without", "drop"]].round(2).to_string(index=False))
    save("learning_curve", lambda: printed(learning_curve, result))
    save("selection_curve", lambda: printed(selection_curve, result))


def pad000040():
    print("PAD000040")
    dt1, dt2 = DATA / "inflam_NASH_dt1.raw.csv", DATA / "inflam_PRA_dt2.raw.csv"
    manifest = DATA / "2025.05.08_sample_manifest.xlsx"
    save("pad_read_plate2", lambda: printed(so.read_npx, dt2))
    plates = pd.concat([quiet(so.read_npx, dt1), quiet(so.read_npx, dt2)], ignore_index=True)
    save("pad_qc_report", lambda: so.qc_report(plates, group="diabetes", samples=manifest, on="FIF ID",
                                               show_ids=False, verbose=False))
    save("pad_apply_qc", lambda: printed(so.apply_qc, plates, show_ids=False))
    clean = quiet(so.apply_qc, plates, show_ids=False)
    save("pad_bridge_check", lambda: so.bridge_check(clean, group="diabetes", samples=manifest,
                                                     on="FIF ID", show_ids=False, verbose=False))
    save("pad_bridge_plates", lambda: printed(so.bridge_plates, clean))
    bridged = quiet(so.bridge_plates, clean)
    save("pad_add_sample_info", lambda: printed(so.add_sample_info, bridged, manifest, on="FIF ID"))
    save("pad_add_sample_info_plates", lambda: printed(so.add_sample_info, plates, manifest, on="FIF ID"))
    data = quiet(so.add_sample_info, plates, manifest, on="FIF ID")
    clean249 = quiet(so.apply_qc, data, drop_samples={"249": "Inflammation PRAESIIDIUM Giraudi 04_03_2025"})
    bridged249 = quiet(so.bridge_plates, clean249)
    save("pad_to_wide", lambda: printed(so.to_wide, bridged249))
    table = quiet(so.to_wide, bridged249)
    order = ["No diabetes", "Prediabetes", "T2DM"]
    save("pad_trend", lambda: printed(compare_groups, table, "diabetes", reference="No diabetes",
                                      ordered=True, order=order))
    two = table[table["diabetes"].isin(["No diabetes", "T2DM"])]
    save("pad_find_signature", lambda: printed(find_signature, two, "diabetes", positive="T2DM",
                                               covariates=["sex"]))


if __name__ == "__main__":
    example()
    if "--example-only" not in sys.argv:
        pad000040()
