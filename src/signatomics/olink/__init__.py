"""Olink extension: read and clean Olink data (Target 96/48, Explore, Explore HT, Reveal).

One function per step; a step that changes data has a function that looks (changes nothing) and
one that acts (returns a new table). Tables found while looking are inside the look result.

1. Read            ``read_npx``         read an export (csv/txt/tsv, Excel, parquet; long or older
                                        wide layout)
2. Match samples   ``add_sample_info``  attach a sample list (clinical data), report what did not
                                        match
3. Check quality   ``qc_report``        look: warnings, empty, repeats, bridging samples, detection
                                        limit (per group), odd samples, group mixed up with plate
                   ``apply_qc``         act: remove low-quality data, every removal reported
4. Bridge plates   ``bridge_check``     look: the plate difference on the bridging samples, and the
                                        shift per protein and plate (``.factors``)
                   ``bridge_plates``    act: put every plate on the scale of the reference plate
5. Analysis table  ``to_wide``        one row per sample, one column per assay, sample metadata
                                        first: the input of ``signatomics.stats``

Charts (Plotly): step 3 ``plot_below_lod``, ``plot_npx_strip``, ``plot_group_by_plate``,
``plot_sample_qc``; step 4 ``plot_bridge_shift``, ``plot_bridge_bland_altman``, ``plot_bridge_scatter``,
``plot_bridge_pca``.

In one call: ``prepare`` runs steps 1-4 with the defaults and prints one report.
To try it all without a file of your own: ``example_file()`` (a synthetic NPX file),
``example_data()`` (the same data as a table) and ``example_samples()`` (its sample metadata).

What the look functions return: ``QCReport`` and ``BridgeCheck``.
"""

from .bridge import BridgeCheck, bridge_check, bridge_plates
from .wide import to_wide
from .example import example_data, example_file, example_samples
from .plots import (plot_below_lod, plot_bridge_bland_altman, plot_bridge_pca, plot_bridge_scatter,
                    plot_bridge_shift,
                    plot_group_by_plate, plot_npx_strip, plot_sample_qc)
from .prepare import prepare
from .qc import QCReport, apply_qc, qc_report
from .reader import read_npx
from .samples import add_sample_info

__all__ = [
    # 1. read, 2. match samples
    "read_npx", "add_sample_info",
    # 3. check quality
    "qc_report", "apply_qc",
    "plot_below_lod", "plot_npx_strip", "plot_group_by_plate", "plot_sample_qc",
    # 4. bridge plates
    "bridge_check", "bridge_plates",
    "plot_bridge_shift", "plot_bridge_bland_altman", "plot_bridge_scatter", "plot_bridge_pca",
    # 5. analysis table
    "to_wide",
    # in one call, and data to try it on
    "prepare", "example_file", "example_data", "example_samples",
    # results of the look functions
    "QCReport", "BridgeCheck",
]
