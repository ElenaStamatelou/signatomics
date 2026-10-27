# signatomics

**Python tools for Olink proteomics data: read, quality control, bridging and sample metadata, with
every step reported.**

signatomics is a family of Python tools for finding biomarker signatures in omics data. The
first module, `signatomics.olink`, works with **Olink (PEA) proteomics** data. More modules will
build on it later.

📖 **Documentation:** <https://signatomics.elenastamatelou.com>

> **Status:** in development (version 0.0.1). The Olink module reads data, controls its quality,
> bridges plates and adds your sample metadata; `signatomics.stats` compares groups and correlates.
> Functions may still change before the first release.

## What it does today

```python
import signatomics.olink as so

data  = so.read_npx(so.example_file())    # a synthetic example NPX file; your own: so.read_npx("my_npx_file.csv")
so.qc_report(data)                        # look at quality: your data stays as it is
clean = so.prepare(so.example_file())     # read, quality control, bridging and sample metadata in one call
```

Each step also exists on its own, with every choice visible: `read_npx`,
`qc_report` → `apply_qc`, `bridge_check` → `bridge_plates`, `add_sample_info`.

- **Reading:** NPX files in long format and older wide (NPX Manager) files; .csv, .txt, .xlsx and .parquet;
  European decimal commas, "No Data" cells, TRUE/FALSE written in other languages.
- **Quality control:** a quality check and charts that leave your data as it is, then your choices
  applied, with every removal listed.
- **Bridging:** per-protein plate correction, with the detection limits shifted along, and an
  independent PCA check.
- **Sample metadata:** adds your sample metadata (demographics, clinical data, ...) and reports every ID
  without a match.
- **Statistics:** the analysis table (`to_matrix`: one row per sample), then `compare_groups` and
  `correlate` on every assay, with effect sizes, 95% intervals and an FDR; `plot_volcano`.
- **Charts:** interactive Plotly charts in one style, colour-blind safe; sample IDs hidden by default.

## Install

```bash
pip install "signatomics[olink] @ git+https://github.com/ElenaStamatelou/signatomics"
```

Requires Python 3.10 or newer. See [Start in 5 minutes](https://signatomics.elenastamatelou.com/start/).

## The family

| Module | Data | Status |
|---|---|---|
| `signatomics.olink` | Olink (PEA) proteomics: Target 96/48, Explore, Explore HT, Reveal | in development |
| `signatomics.stats` | compare groups, correlate: every assay, effect sizes, FDR | in development |
| core | one standard table, combining data types, biomarker signature with nested cross-validation | planned |
| more modules | e.g. ELISA, microbiome | ideas |

## Repository layout

| Folder | What |
|---|---|
| `src/signatomics/` | the package |
| `tests/` | automatic tests (`pytest`) |
| `analyses/` | worked analyses of public data (Jupyter notebooks) |
| `website/` | the documentation website (MkDocs Material) |
| `docs/` | function reference pointer, notes on the example datasets |

Public datasets stay outside this repository; the only data inside are the small, synthetic example
files of the package. See [docs/datasets.md](docs/datasets.md).

## Author

[Elena Stamatelou](https://elenastamatelou.com)

*signatomics is an independent project. It is not affiliated with or endorsed by Olink or Thermo
Fisher Scientific. Olink is a trademark of Thermo Fisher Scientific.*
