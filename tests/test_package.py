"""The simplest possible test: can Python find and load the package?"""

import signatomics
import signatomics.olink


def test_package_has_a_version():
    assert signatomics.__version__ == "0.0.1"


def test_olink_public_functions_one_per_step():
    """The public list is deliberately short: helpers live as tables inside the reports."""
    assert set(signatomics.olink.__all__) == {
        "read_npx", "add_sample_info",
        "qc_report", "apply_qc",
        "plot_below_lod", "plot_npx_strip", "plot_group_by_plate", "plot_sample_qc",
        "bridge_check", "bridge_plates",
        "plot_bridge_shift", "plot_bridge_bland_altman", "plot_bridge_scatter", "plot_bridge_pca",
        "prepare", "example_file", "example_data", "example_samples",
        "QCReport", "BridgeCheck",
    }
    for name in signatomics.olink.__all__:
        assert hasattr(signatomics.olink, name)
