# Datasets

All from **PRIDE** (EBI public archive, affinity proteomics section), licence **CC0 (public
domain)**: free to use for anything. Still: cite the dataset and its paper in every tutorial,
and never try to identify participants. Downloaded files go in `data/<ID>/` (never saved into
the history). Download links: `https://ftp.pride.ebi.ac.uk/pride/data/archive/<year>/<month>/<ID>/`.

## 1. PAD000040 — in use (downloaded 2 Oct 2026)

**Inflammatory proteomic heterogeneity beyond glycemia status in severe obesity.**
Milito MM, ..., Giraudi PJ. *Int J Mol Sci* 2026;27(9):4152. doi:10.3390/ijms27094152.
PRIDE: https://www.ebi.ac.uk/pride/archive/projects/PAD000040

- **Olink Target 96 Inflammation**, 92 proteins, **142 people**, 2 plates from 2 studies
  (NASH study, plate 21/03/2024; PRAESIIDIUM study, plate 04/03/2025).
- **Outcome labels included:** diabetes status — 69 no diabetes, 25 prediabetes, 48 type 2
  diabetes; also sex and liver histology (fibrosis, inflammation, ballooning, steatosis, MASH).
- Files:
  | File | What |
  |---|---|
  | `dataset_all.npx.csv` | Final NPX values, both plates combined (has an adjustment factor → already bridged) |
  | `inflam_NASH_dt1.raw.csv`, `inflam_PRA_dt2.raw.csv` | The two plates separately, before combining |
  | `2025.05.08_sample_manifest.xlsx` | 142 samples with clinical labels |
  | `samplemanifest_bis.xlsx` | 193 rows incl. bridging samples and study labels |
  | `PAD000040_community_annotated.sdrf.tsv` | Standard sample description (incl. technical replicates) |
- **What we found (notebook `analyses/01_PAD000040_first_look.ipynb`):** comma-separated;
  Italian spreadsheet values ("VERO"/"FALSO" = TRUE/FALSE); `MissingFreq` mixes decimals and
  percentages; protein file `SampleID` matches sample list `FIF ID` — all 142 after tidying IDs
  (`308_BIS` vs `308 BIS`); 177 measurements: 34 bridging samples on both plates, sample 309
  rerun on the same plate; 17 measurements with a QC warning, 2 completely empty (196, 292);
  14 of 92 proteins below LOD in >50% of measurements. **Prediabetes was never measured on
  plate 1 alone** → group partly mixed up with plate (confounding by batch).
- **Good for:** reading Olink csv exports, QC report, detection-limit handling, checking the
  bridging of the two plates (raw vs combined files), a first biomarker signature for diabetes.

## 2. For later (not downloaded)

| ID | What | Format / size | Use |
|---|---|---|---|
| **PAD000062** | Down syndrome, tofacitinib treatment | Explore HT: parquet 27 MB + csv 327 MB + QC file | Newest format, ~5,400 proteins |
| **PAD000009** | NIST reference plasma + spike-in protein standards | Reveal: csv ~13 MB + analysis report | Reveal format; known added amounts = built-in truth |
| PAD000002 | Technical evaluation of plasma proteomics technologies | Explore HT: parquet 20 MB | Format testing |

Other public sources: MGH COVID-19 study (Explore 1536, registration
needed), Figshare "Raw Olink data" (Explore 1536 + Expansion, CC BY 4.0, small).
