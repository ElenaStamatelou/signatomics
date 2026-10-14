# Function reference — moved

The function reference lives in the documentation website, one page per step:
<https://signatomics.elenastamatelou.com/reference/>. The source pages are in this repository:

| Page | Functions |
|---|---|
| [website/docs/reference/index.md](../website/docs/reference/index.md) | overview by step, rules, the standard table |
| [reading.md](../website/docs/reference/reading.md) | `read_npx` |
| [samples.md](../website/docs/reference/samples.md) | `add_sample_info` |
| [quality.md](../website/docs/reference/quality.md) | `qc_report`, `apply_qc` and the quality charts |
| [bridging.md](../website/docs/reference/bridging.md) | `bridge_check`, `bridge_plates` and the bridging charts |

Rule: when a function changes, its docstring and its reference page change in the same step.
Preview the website locally with `mkdocs serve -f website/mkdocs.yml`.
