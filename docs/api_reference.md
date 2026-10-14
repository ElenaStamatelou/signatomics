# Function reference — moved

The function reference lives in the documentation website, one page per step:
<https://signatomics.elenastamatelou.com/reference/>. The source pages are in this repository:

| Page | Functions |
|---|---|
| [website/docs/reference/index.md](../website/docs/reference/index.md) | overview, rules, the standard table |
| [reading.md](../website/docs/reference/reading.md) | `read_npx`, `tidy_sample_ids` |
| [samples.md](../website/docs/reference/samples.md) | `add_sample_info` |
| [quality.md](../website/docs/reference/quality.md) | `qc_report`, `apply_qc`, `below_lod_by_group`, `sample_qc_table` |
| [bridging.md](../website/docs/reference/bridging.md) | `bridge_check`, `bridge_factors`, `bridge_plates`, `bridge_differences` |
| [charts.md](../website/docs/reference/charts.md) | all `plot_*` functions |

Rule: when a function changes, its docstring and its reference page change in the same step.
Preview the website locally with `mkdocs serve -f website/mkdocs.yml`.
