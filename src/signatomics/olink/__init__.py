"""Olink extension: read and clean Olink data (Target 96/48, Explore, Explore HT, Reveal).

Available:
- ``read_npx``         read an export (csv/txt/tsv, Excel, parquet; long or older wide layout)
- ``tidy_sample_ids``  make sample IDs comparable between files ("308 BIS" == "308_BIS")
- ``add_sample_info``  attach a sample list (clinical data) and report what did not match
- ``qc_report``        look at data quality (warnings, empty, repeats, bridging samples,
                       detection limit, group mixed up with plate) without changing anything
- ``apply_qc``         remove low-quality data — every choice visible, every removal reported
- ``bridge_check``     look at the difference between plates, on the bridging samples
- ``bridge_factors``   the shift per protein and plate (+ ``bridge_differences`` per sample)
- ``bridge_plates``    put every plate on the scale of the reference plate
- charts (Plotly):     ``plot_below_lod``, ``plot_npx_strip``, ``plot_group_by_plate``,
                       ``plot_sample_qc`` (+ ``sample_qc_table``), ``plot_bridge_shift``,
                       ``plot_bridge_pca``, ``plot_bridge_scatter``
"""

from .bridge import (BridgeCheck, bridge_check, bridge_differences, bridge_factors,
                     bridge_plates)
from .plots import (plot_below_lod, plot_bridge_pca, plot_bridge_scatter, plot_bridge_shift,
                    plot_group_by_plate, plot_npx_strip, plot_sample_qc, sample_qc_table)
from .qc import QCReport, apply_qc, below_lod_by_group, qc_report
from .reader import read_npx, tidy_sample_ids
from .samples import add_sample_info

__all__ = ["read_npx", "tidy_sample_ids", "add_sample_info", "qc_report", "apply_qc", "QCReport",
           "plot_below_lod", "plot_npx_strip", "plot_group_by_plate", "plot_sample_qc",
           "sample_qc_table", "below_lod_by_group",
           "bridge_check", "BridgeCheck", "bridge_differences", "bridge_factors", "bridge_plates",
           "plot_bridge_shift", "plot_bridge_pca", "plot_bridge_scatter"]
