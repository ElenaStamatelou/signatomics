"""Biomarker signatures: a small set of proteins that predicts an outcome, with an honest estimate.

Make the analysis table with ``signatomics.olink.to_wide`` and pick the outcome column.

- ``find_signature``         nested cross-validation for a binary, multiclass or continuous outcome;
                             the proteins are chosen inside every training part
- ``protein_contributions``  what each signature protein adds: alone, and the panel without it
- ``learning_curve``         performance as the training part grows
- ``selection_curve``        performance for 1, 2, 3 ... proteins

Charts (Plotly): ``plot_signature_performance``, ``plot_selection_frequency``, ``plot_curve``
(learning and selection curve), ``plot_contributions``.
"""

from .checks import CheckResult, learning_curve, protein_contributions, selection_curve
from .find import SignatureResult, find_signature
from .plots import plot_contributions, plot_curve, plot_selection_frequency, plot_signature_performance

__all__ = ["find_signature", "protein_contributions", "learning_curve", "selection_curve",
           "plot_signature_performance", "plot_selection_frequency", "plot_curve", "plot_contributions",
           "SignatureResult", "CheckResult"]
