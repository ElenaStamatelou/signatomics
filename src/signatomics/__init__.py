"""signatomics — multi-omics biomarker discovery.

From raw omics files to a biomarker signature for a chosen outcome:

1. read and clean each data type with its own rules (extensions, e.g. ``signatomics.olink``)
2. bring everything into one standard table and combine data types
3. find a small, honest biomarker signature (the core)
"""

__version__ = "0.0.1"
