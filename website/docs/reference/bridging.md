---
title: "Bridging Olink plates in Python"
description: "Bridge Olink plates and projects in Python: per-protein plate shifts from bridging samples, LOD shifted along, independent checks and Plotly charts."
---

# Bridging Olink plates

When a study is measured on several plates, or in several projects, the plates read slightly
differently. **Bridging samples**, the same tube measured on both plates, show by how much.
{ .lede }

For each protein, the **median** difference over the bridging samples is the plate shift. It is
added to every sample on the other plate, which puts that plate on the scale of the **reference
plate**. The reference plate itself is not changed.

NPX is on a log2 scale, so bridging is an *addition*: a difference of 1 NPX is a factor 2. Only
NPX data can be bridged. Skip this step when all samples are on one plate.

Like the quality step, bridging has two halves: **look** ([`bridge_check`](#bridge_check) and the
[charts](#charts-for-this-step)), then **act** ([`bridge_plates`](#bridge_plates)).

**Recommended order:** `apply_qc` → `bridge_check` → `bridge_plates` → `bridge_check` again,
plus `plot_bridge_pca`.

## Shared arguments

All bridging functions take these arguments:

| Argument | Type | Default | Meaning |
|---|---|---|---|
| `batch` | str | `"PlateID"` | the column that tells the plates apart. Use e.g. `"Project"` to bridge whole projects whose plate names may overlap |
| `reference` | plate name | `None` | the plate that stays as it is; default the plate with the most samples |
| `bridges` | list of IDs | `None` | the samples to use as bridges; default every sample on the reference *and* the other plate |
| `exclude_warnings` | bool | `True` | do not use bridging measurements with a quality warning |
| `exclude_below_lod` | bool | `False` | use only pairs where both values are above LOD |

Bridging samples are matched by tidy ID, so `"B 01"` on one plate matches `"B_01"` on the
other. A sample measured twice on the same plate is averaged first. Control samples and control
assays are never used.

---

## `bridge_check`

Look at the difference between plates using the bridging samples. It changes nothing. Run it
before bridging (how big is the difference?) and after (is it gone?).

```python
bridge_check(npx, *, batch="PlateID", reference=None, group=None, bridges=None,
             exclude_warnings=True, exclude_below_lod=False, min_bridges=8, tolerance=0.1,
             show_ids=True, verbose=True) -> BridgeCheck
```

**Inputs.** The [shared arguments](#shared-arguments), plus:

| Argument | Type | Default | Meaning |
|---|---|---|---|
| `npx` | DataFrame | — | all plates in one table (`pd.concat` separate files) |
| `group` | str | `None` | a sample-information column: checks whether the bridging samples cover every group, and whether the shift is the same in each group |
| `min_bridges` | int | `8` | fewer bridging samples than this gives a warning; proteins with fewer usable pairs are flagged in `.factors` |
| `tolerance` | float | `0.1` | a protein whose shift is larger than this (NPX) counts as "still shifted" |
| `show_ids` | bool | `True` | list sample IDs |
| `verbose` | bool | `True` | print the report |

**Output.** A `BridgeCheck` object, which prints as a report:

| Attribute | What it holds |
|---|---|
| `.lines` | the report lines |
| `.differences` | one row per bridging sample, protein and non-reference plate: `sample` (the tidy ID), the protein columns, `plate`, `reference`, `NPX_reference`, `NPX_plate` and `difference` (= plate − reference). The reference used is in `.attrs["reference"]` |
| `.proteins` | per plate and protein: `shift` (median difference), `spread` (IQR), `bridges` (number of pairs) |
| `.samples` | per plate and bridging sample: `shift` (median over proteins), to find odd bridging samples |
| `.factors` | per plate and protein: the `adjustment` that `bridge_plates` adds (see below) |
| `.by_group` | with `group`: per protein, the shift in each group and their `difference` |

The **factors** table holds, for each protein and each plate other than the reference:

```text
adjustment = − median( NPX on the plate − NPX on the reference )     over the bridging samples
```

| Column | What it holds |
|---|---|
| `plate`, `reference` | the plate that will be adjusted, and the reference |
| `OlinkID`, `Assay` | the protein |
| `adjustment` | NPX to **add** to this protein on this plate (empty if no bridging pair) |
| `bridges` | number of bridging pairs used |
| `spread` | IQR of the bridge differences: noise that bridging cannot remove |
| `few_bridges` | True when `bridges < min_bridges` |

You can inspect it, save it (`.to_csv`) and apply it later with
`bridge_plates(..., factors=...)`.

For each plate, the report says:
- how many bridging samples there are;
- the typical shift and its range over proteins;
- how many proteins are shifted by more than `tolerance`;
- the noise (the typical spread);
- **which bridging samples shift differently from the others**: more than 3 robust SDs and
  0.25 NPX away, which can mean a mix-up or a degraded tube;
- with `group`, the shift in each group.

!!! note "Checking after bridging"
    On a table that is already bridged, the shifts are about 0 *by construction*, because the
    same samples were used to compute the adjustment. The report says so. The independent check
    is [`plot_bridge_pca`](#plot_bridge_pca).

**Errors**

- `ValueError`: there is no NPX column, or only one plate.
- `KeyError`: `batch` or `group` is not a column.
- `ValueError`: `reference` is not one of the plates.

The table `.differences` is empty, and the report warns, when no sample is on two plates.

**Example** (PAD000040, the two raw plate files)

```python
chk = so.bridge_check(clean, group="diabetes", show_ids=False)
chk.factors.sort_values("adjustment").head()
```

```text
Bridging check
  reference: Inflammation F.I.F 21/03/2024
  Inflammation PRAESIIDIUM Giraudi 04_03_2025: 29 bridging samples; reads +0.52 NPX vs the
    reference (median over proteins; per protein -0.96 to +2.30; 1.44 times)
    proteins shifted by more than 0.1 NPX: 83 of 92 - bridging needed
    noise: the bridge differences of one protein typically spread 0.44 NPX (IQR)
    bridging samples that shift differently from the others: 2 (IDs hidden)
  'diabetes' among the bridging samples: No diabetes: 9 of 62; T2DM: 9 of 46; Prediabetes: 10 of 23
    overall shift per group: No diabetes +0.56; Prediabetes +0.57; T2DM +0.53 (largest difference 0.04 NPX)
```

---

## `bridge_plates`

Put every plate on the scale of the reference plate.

```python
bridge_plates(npx, *, factors=None, batch="PlateID", reference=None, bridges=None,
              exclude_warnings=True, exclude_below_lod=False, min_bridges=8, verbose=True)
    -> pandas.DataFrame
```

**Inputs.** The [shared arguments](#shared-arguments), plus:

| Argument | Type | Default | Meaning |
|---|---|---|---|
| `npx` | DataFrame | — | all plates in one table, best after `apply_qc` |
| `factors` | DataFrame | `None` | a factors table to apply, e.g. `bridge_check(...).factors` checked or saved earlier; computed here when `None` |
| `min_bridges` | int | `8` | proteins with fewer pairs are named in the report |
| `verbose` | bool | `True` | print the report |

**Output.** The same table, with:

- `NPX` **and** `LOD`, `MaxLOD` shifted by the adjustment, so a value below LOD stays below LOD;
- `NPX_unbridged` and `LOD_unbridged`: the original values;
- `Bridge_adjustment`: the NPX added to the row. It is 0 on the reference plate, and empty if the
  protein had no bridging pair; in that case the value is left unchanged and the report warns.

The report is stored in `df.attrs["bridge_report"]`.

**Errors**

- `ValueError`: the table is already bridged (it has a `Bridge_adjustment` column). Bridge the
  original table instead.
- `ValueError`: there is no `NPX` column, or only one plate.
- `KeyError`: `batch` is not a column.
- `ValueError`: `reference` is not one of the plates.

**Example**

```python
bridged = so.bridge_plates(clean)
```

```text
Bridging report
  reference (not changed): Inflammation F.I.F 21/03/2024
  Inflammation PRAESIIDIUM Giraudi 04_03_2025: 92 proteins adjusted, median -0.52 NPX
      (from -2.30 to +0.96); bridging samples per protein: 29-29
  NPX, LOD and MaxLOD shifted together; originals in NPX_unbridged / LOD_unbridged
  next: bridge_check(...) on the result - the shifts should now be about 0
```

---

## Charts for this step

Every chart returns a Plotly figure; see
[the rules for charts](index.md#rules-that-hold-for-every-function).

### `plot_bridge_shift`

For each protein, how much higher the other plate reads than the reference, for the same samples:
- grey dots: single bridging samples;
- coloured marker: the median, which is the shift `bridge_plates` removes;
- bar: the middle 50%.

Before bridging the markers sit away from 0; after bridging they sit on 0.

```python
plot_bridge_shift(npx, *, batch="PlateID", reference=None, bridges=None, exclude_warnings=True,
                  exclude_below_lod=False, show_ids=False) -> Figure
```

Arguments as in the [shared arguments](#shared-arguments). At most 3 plates besides the reference.

**Errors.** `ValueError` when there are no bridging samples, or more than 3 other plates.

### `plot_bridge_scatter`

Bridging samples, with NPX on the reference plate (x) against the other plate (y), one panel per
protein. On the dashed diagonal both plates agree; dotted lines mark each plate's LOD.

```python
plot_bridge_scatter(npx, proteins, *, batch="PlateID", reference=None, exclude_warnings=True,
                    show_ids=False, cols=3) -> Figure
```

| Argument | Default | Meaning |
|---|---|---|
| `proteins` | — | 1 to 12 protein names (`Assay`) |
| `cols` | `3` | panels per row |

**Errors.** `ValueError` for 0 or more than 12 proteins, or a protein without bridging pairs.

### `plot_bridge_pca`

A PCA of the bridging samples, coloured by plate, before and (optionally) after bridging. Each
bridging sample appears once per plate, and a thin line joins its measurements. After bridging,
the plates should overlap and the lines shrink.

**This is the independent check of bridging.** The shift was fitted per protein; the PCA looks at
all proteins together. The subtitle gives the median distance between a sample's measurements,
before and after.

```python
plot_bridge_pca(before, after=None, *, batch="PlateID", samples="bridges",
                exclude_warnings=True, show_ids=False) -> Figure
```

| Argument | Default | Meaning |
|---|---|---|
| `before` | — | table before bridging |
| `after` | `None` | the same table after `bridge_plates`; adds a second panel |
| `samples` | `"bridges"` | `"bridges"` (only samples on more than one plate) or `"all"` |

Proteins with more than half of their values missing are left out. Other missing values are
filled with the protein's median, and proteins are scaled before the PCA.

**Errors.** `ValueError` with fewer than 3 measurements, or an unknown `samples` value.

!!! example "On PAD000040"
    Before bridging, a sample's two measurements lie **6.5** apart (median). After bridging they
    lie **1.8** apart.
