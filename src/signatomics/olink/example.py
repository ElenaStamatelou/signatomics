"""Synthetic example data, so every step can be tried without an NPX file of your own.

The values are synthetic (no person, no study) and the same every time. They look like a small
Olink Target 96 study: 60 people (30 cases, 30 controls), 92 assays (named after the
Target 96 Inflammation panel), measured on two plates, one NPX file per plate: P01-P40 on
plate 1, P29-P60 on plate 2. People P29-P40 are on both plates: the bridge samples. Plate 2
reads higher than plate 1, by a different amount per assay. Each plate also holds the 8
controls of a Target 96 plate (3 negative, 3 plate and 2 sample controls), marked in
``SampleType`` as in Olink's exports.

Built in on purpose, so the checks have something to find:

- Plate 1: P05 with QC status WARN, P20 with FAIL, P07 far too low on every assay (an outlier);
- Plate 2: P41 with QC status WARN, P33 far too low on every assay (a broken bridge sample);
- ten assays mostly below the LOD, and a few around it;
- IL-17C detectable mostly in cases;
- four assays higher in cases (IL6, CXCL10, HGF, OSM).
"""

from __future__ import annotations

from importlib.resources import files

import numpy as np
import pandas as pd

PROTEINS = [
    "IL8", "VEGFA", "CD8A", "MCP-3", "GDNF", "CDCP1", "CD244", "IL7", "OPG", "LAP TGF-beta-1",
    "uPA", "IL6", "IL-17C", "MCP-1", "IL-17A", "CXCL11", "AXIN1", "TRAIL", "IL-20RA", "CXCL9",
    "CST5", "IL-2RB", "IL-1 alpha", "OSM", "IL2", "CXCL1", "TSLP", "CCL4", "CD6", "SCF",
    "IL18", "SLAMF1", "TGF-alpha", "MCP-4", "CCL11", "TNFSF14", "FGF-23", "IL-10RA", "FGF-5",
    "MMP-1", "LIF-R", "FGF-21", "CCL19", "IL-15RA", "IL-10RB", "IL-22 RA1", "IL-18R1", "PD-L1",
    "Beta-NGF", "CXCL5", "TRANCE", "HGF", "IL-12B", "IL-24", "IL13", "ARTN", "MMP-10", "IL10",
    "TNF", "CCL23", "CD5", "CCL3", "Flt3L", "CXCL6", "CXCL10", "4E-BP1", "IL-20", "SIRT2",
    "CCL28", "DNER", "EN-RAGE", "CD40", "IL33", "IFN-gamma", "FGF-19", "IL4", "LIF", "NRTN",
    "MCP-2", "CASP-8", "CCL25", "CX3CL1", "TNFRSF9", "NT-3", "TWEAK", "CCL20", "ST1A1",
    "STAMBP", "IL5", "ADA", "TNFB", "CSF-1",
]
# mostly below the detection limit, as on the real panel
LOW = {"IL-1 alpha", "IL2", "IL-2RB", "IL4", "IL5", "IL13", "IL-20", "IL33", "TSLP", "NRTN"}
NEAR_LOD = {"IL-20RA", "IL-22 RA1", "IL-24", "ARTN", "Beta-NGF"}    # around the LOD
HIGHER_IN_CASES = {"IL6": 1.0, "CXCL10": 0.8, "HGF": 0.7, "OSM": 0.9, "IL-17C": 1.4}
PLATES = ("Plate 1", "Plate 2")
CONTROLS = ([("NEGATIVE_CONTROL", f"NEG CTRL {i}") for i in (1, 2, 3)]
            + [("PLATE_CONTROL", f"PLATE CTRL {i}") for i in (1, 2, 3)]
            + [("SAMPLE_CONTROL", f"SAMPLE CTRL {i}") for i in (1, 2)])


def _people():
    ids = [f"P{k:02d}" for k in range(1, 61)]
    group = ["case" if k % 2 == 0 else "control" for k in range(1, 61)]
    return ids, group


def example_samples() -> pd.DataFrame:
    """The sample metadata that go with the example: one row per person.

    Columns: ``sample_id``, ``group`` (``"case"`` / ``"control"``), ``age``, ``sex``. Synthetic values.
    """
    rng = np.random.default_rng(7)
    ids, group = _people()
    return pd.DataFrame({
        "sample_id": ids,
        "group": group,
        "age": rng.integers(35, 76, len(ids)),
        "sex": rng.choice(["female", "male"], len(ids)),
    })


def example_data() -> pd.DataFrame:
    """Both example plates in one table, in long format, as ``read_npx`` returns them.

    The same values as ``read_npx(example_file(1))`` and ``read_npx(example_file(2))`` stacked.
    Columns: ``SampleID``, ``SampleType``, ``Index``, ``OlinkID``, ``Assay``, ``Panel``,
    ``PlateID``, ``QC_Warning``, ``LOD``, ``MissingFreq``, ``NPX``.
    """
    rng = np.random.default_rng(2028)
    ids, group = _people()
    is_case = dict(zip(ids, (g == "case" for g in group)))
    n = len(PROTEINS)
    level, lod1 = {}, {}
    for p in PROTEINS:
        if p in LOW:
            level[p], lod1[p] = rng.uniform(0.8, 1.6), 2.4 + rng.uniform(-0.2, 0.2)
        elif p in NEAR_LOD or p == "IL-17C":
            level[p] = rng.uniform(2.0, 3.0)
            lod1[p] = level[p] + rng.uniform(0.0, 0.4)
        else:
            level[p] = rng.uniform(3.0, 9.0)
            lod1[p] = level[p] - 2.6 + rng.uniform(0, 0.6)
    spread = {p: rng.uniform(0.5, 1.1) for p in PROTEINS}
    shift = dict(zip(PROTEINS, rng.normal(0.5, 0.45, n)))           # plate 2 minus plate 1
    person = {(s, p): rng.normal(0, spread[p]) for s in ids for p in PROTEINS}

    plan = ([(s, "Plate 1", "SAMPLE") for s in ids[:40]]
            + [(c, "Plate 1", t) for t, c in CONTROLS]
            + [(s, "Plate 2", "SAMPLE") for s in ids[28:]]
            + [(c, "Plate 2", t) for t, c in CONTROLS])
    status = {("P05", "Plate 1"): "WARN", ("P20", "Plate 1"): "FAIL", ("P41", "Plate 2"): "WARN"}
    broken = {("P07", "Plate 1"), ("P33", "Plate 2")}               # far too low on every assay

    rows, index = [], {"Plate 1": 0, "Plate 2": 0}
    for s, plate, kind in plan:
        index[plate] += 1
        on2 = plate == "Plate 2"
        for k, p in enumerate(PROTEINS):
            base = level[p] + (shift[p] if on2 else 0.0)
            if kind == "SAMPLE":
                npx = (base + person[(s, p)] + (HIGHER_IN_CASES.get(p, 0.0) if is_case[s] else 0.0)
                       + rng.normal(0, 0.25))
                if (s, plate) in broken:
                    npx -= 5.0
            elif kind == "NEGATIVE_CONTROL":                         # background: below the LOD
                npx = lod1[p] + (shift[p] if on2 else 0.0) - 1.5 + rng.normal(0, 0.3)
            else:                                                    # plate and sample controls
                npx = base + rng.normal(0, 0.15)
            rows.append({
                "SampleID": s, "SampleType": kind, "Index": index[plate],
                "OlinkID": f"EX{k + 1:04d}", "Assay": p, "Panel": "Inflammation (synthetic)",
                "PlateID": plate, "QC_Warning": status.get((s, plate), "PASS"),
                "LOD": round(lod1[p] + (shift[p] if on2 else 0.0), 4), "NPX": round(npx, 4),
            })
    df = pd.DataFrame(rows)
    # MissingFreq: the share of samples below the LOD, per assay and plate (controls left out)
    is_sample = df["SampleType"] == "SAMPLE"
    below = (df["NPX"] < df["LOD"])[is_sample].groupby([df["Assay"], df["PlateID"]]).mean()
    freq = pd.MultiIndex.from_frame(df[["Assay", "PlateID"]]).map(below.round(4).to_dict())
    df.insert(df.columns.get_loc("NPX"), "MissingFreq", freq)
    return df


def example_file(plate: int = 1) -> str:
    """Path to a small, synthetic NPX file that comes with the package: plate 1 or plate 2.

    Use it like your own file: ``read_npx(example_file())`` reads plate 1;
    ``read_npx(example_file(2))`` reads plate 2 (for bridging).
    """
    if plate not in (1, 2):
        raise ValueError("The example has two plates: example_file(1) or example_file(2).")
    return str(files("signatomics.olink") / "data" / f"example_plate{plate}.csv")
