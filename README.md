# signatomics

**Python tools for Olink proteomics data: read, check, bridge and analyse, with every step
reported.**

signatomics is a family of Python tools for finding biomarker signatures in omics data. The
first module, `signatomics.olink`, works with **Olink (PEA) proteomics** data. More modules will
build on it later.

📖 **Documentation:** <https://signatomics.elenastamatelou.com>

> **Status:** in development (version 0.0.1). The Olink module reads, matches, quality-checks
> and bridges data. Functions may still change before the first release.

## What it does today

```python
import pandas as pd
import signatomics.olink as so

npx     = pd.concat([so.read_npx("plate1.csv"), so.read_npx("plate2.csv")], ignore_index=True)
data    = so.add_sample_info(npx, "samples.xlsx", on="Sample ID")
so.qc_report(data, group="diabetes")          # look
clean   = so.apply_qc(data)                   # decide
bridged = so.bridge_plates(clean)             # line up the plates
so.plot_bridge_pca(clean, bridged)            # check that it worked
```

- **Reading:** current long exports and older wide (NPX Manager) files; csv, Excel and parquet;
  European decimal commas, "No Data" cells, TRUE/FALSE written in other languages.
- **Sample information:** attaches your clinical data and reports every ID that did not match.
- **Quality:** a report that changes nothing, then a quality step where every removal is listed.
- **Bridging:** per-protein plate correction, with the detection limits shifted along, and an
  independent PCA check.
- **Charts:** interactive Plotly charts, colour-blind safe; sample IDs hidden by default.

## Install

```bash
pip install "signatomics[olink] @ git+https://github.com/ElenaStamatelou/signatomics"
```

Requires Python 3.10 or newer. See [Getting started](https://signatomics.elenastamatelou.com/getting-started/).

## The family

| Module | Data | Status |
|---|---|---|
| `signatomics.olink` | Olink (PEA) proteomics: Target 96/48, Explore, Explore HT, Reveal | in development |
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

Example data is public and is never stored in this repository. See [docs/datasets.md](docs/datasets.md).

## Author

[Elena Stamatelou](https://elenastamatelou.com)

*signatomics is an independent project. It is not affiliated with or endorsed by Olink or Thermo
Fisher Scientific. Olink is a trademark of Thermo Fisher Scientific.*
