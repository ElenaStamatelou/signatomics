"""Made-up example data, so every step can be tried without an NPX file of your own.

The values are invented (no person, no study) and the same every time. They look like a small
Olink Target 96 study: 60 people (30 cases, 30 controls), 24 assays, measured on two plates,
one NPX file per plate: P01-P40 on plate 1, P29-P60 on plate 2. People P29-P40 are on both
plates: the bridge samples. Plate 2 reads
higher than plate 1, by a different amount per assay.

Built in on purpose, so the checks have something to find:

- Plate 1: P05 with QC status WARN, P20 with FAIL, P07 far too low on every assay (an outlier);
- Plate 2: P41 with QC status WARN, P33 far too low on every assay (a broken bridge sample);
- three assays mostly below the LOD (IL13, IL4, IL5);
- four assays higher in cases (IL6, CXCL10, HGF, OSM).
"""

from __future__ import annotations

from importlib.resources import files

import numpy as np
import pandas as pd

PROTEINS = ["IL6", "IL8", "IL10", "IL-18", "TNF", "MCP-1", "CCL3", "CCL4", "CXCL9", "CXCL10",
            "CSF-1", "HGF", "VEGFA", "TGF-alpha", "OSM", "FGF-21", "IL-17A", "IL13", "IL4", "IL5",
            "IL2", "IFN-gamma", "CDCP1", "ADA"]
LOW = {"IL13", "IL4", "IL5"}                     # mostly below the detection limit
HIGHER_IN_CASES = {"IL6": 1.0, "CXCL10": 0.8, "HGF": 0.7, "OSM": 0.9}
PLATES = ("Plate 1", "Plate 2")


def _people():
    ids = [f"P{k:02d}" for k in range(1, 61)]
    group = ["case" if k % 2 == 0 else "control" for k in range(1, 61)]
    return ids, group


def example_samples() -> pd.DataFrame:
    """The sample metadata that go with the example: one row per person.

    Columns: ``sample_id``, ``group`` (``"case"`` / ``"control"``), ``age``, ``sex``. Made-up values.
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
    Columns: ``SampleID``, ``Index``, ``OlinkID``, ``Assay``, ``Panel``, ``PlateID``,
    ``QC_Warning``, ``LOD``, ``MissingFreq``, ``NPX``.
    """
    rng = np.random.default_rng(2026)
    ids, group = _people()
    is_case = dict(zip(ids, (g == "case" for g in group)))
    n = len(PROTEINS)
    level = {p: (rng.uniform(0.8, 1.8) if p in LOW else rng.uniform(3.0, 9.0)) for p in PROTEINS}
    spread = {p: rng.uniform(0.5, 1.1) for p in PROTEINS}
    shift = dict(zip(PROTEINS, rng.normal(0.5, 0.45, n)))           # plate 2 minus plate 1
    lod1 = {p: (2.0 + rng.uniform(-0.2, 0.2) if p in LOW else level[p] - 2.6 + rng.uniform(0, 0.6))
            for p in PROTEINS}
    person = {(s, p): rng.normal(0, spread[p]) for s in ids for p in PROTEINS}

    plan = [(s, "Plate 1") for s in ids[:40]] + [(s, "Plate 2") for s in ids[28:]]
    status = {("P05", "Plate 1"): "WARN", ("P20", "Plate 1"): "FAIL", ("P41", "Plate 2"): "WARN"}
    broken = {("P07", "Plate 1"), ("P33", "Plate 2")}               # far too low on every assay

    rows, index = [], {"Plate 1": 0, "Plate 2": 0}
    for s, plate in plan:
        index[plate] += 1
        on2 = plate == "Plate 2"
        for k, p in enumerate(PROTEINS):
            npx = (level[p] + person[(s, p)] + (HIGHER_IN_CASES.get(p, 0.0) if is_case[s] else 0.0)
                   + (shift[p] if on2 else 0.0) + rng.normal(0, 0.25))
            if (s, plate) in broken:
                npx -= 5.0
            rows.append({
                "SampleID": s, "Index": index[plate], "OlinkID": f"EX{k + 1:04d}", "Assay": p,
                "Panel": "Inflammation (example)", "PlateID": plate,
                "QC_Warning": status.get((s, plate), "PASS"),
                "LOD": round(lod1[p] + (shift[p] if on2 else 0.0), 4), "NPX": round(npx, 4),
            })
    df = pd.DataFrame(rows)
    below = (df["NPX"] < df["LOD"]).groupby([df["Assay"], df["PlateID"]]).transform("mean")
    df.insert(df.columns.get_loc("NPX"), "MissingFreq", below.round(4))
    return df


def example_file(plate: int = 1) -> str:
    """Path to a small, generated NPX file (made-up values) that comes with the package: plate 1 or plate 2.

    Use it like your own file: ``read_npx(example_file())`` reads plate 1;
    ``read_npx(example_file(2))`` reads plate 2 (for bridging).
    """
    if plate not in (1, 2):
        raise ValueError("The example has two plates: example_file(1) or example_file(2).")
    return str(files("signatomics.olink") / "data" / f"example_plate{plate}.csv")
