---
title: "Bridging Olink plates in Python"
description: "Bridge Olink plates and projects in Python: per-protein plate shifts from bridging samples, LOD shifted along, and independent checks."
---

# Bridging Olink plates

When a study is measured on several plates, or in several projects, the plates read slightly
differently. **Bridging samples**, the same tube measured on both plates, show by how much.
{ .lede }

For each protein, the **median** difference over the bridging samples is the plate shift. It is
added to every sample on the other plate, which puts that plate on the scale of the **reference
plate**. The reference plate itself is not changed.

NPX is on a log2 scale, so bridging is an *addition*: a difference of 1 NPX is a factor 2. Only
NPX data can be bridged.

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
             exclude_warnings=True, exclude_below_lod=False, tolerance=0.1,
             show_ids=True, verbose=True) -> BridgeCheck
```

**Inputs.** The [shared arguments](#shared-arguments), plus:

| Argument | Type | Default | Meaning |
|---|---|---|---|
| `npx` | DataFrame | — | all plates in one table (`pd.concat` separate files) |
| `group` | str | `None` | a sample-information column: checks whether the bridging samples cover every group, and whether the shift is the same in each group |
| `tolerance` | float | `0.1` | a protein whose shift is larger than this (NPX) counts as "still shifted" |
| `show_ids` | bool | `True` | list sample IDs |
| `verbose` | bool | `True` | print the report |

**Output.** A `BridgeCheck` object, which prints as a report:

| Attribute | What it holds |
|---|---|
| `.lines` | the report lines |
| `.differences` | one row per bridging sample, protein and plate (see [`bridge_differences`](#bridge_differences)) |
| `.proteins` | per plate and protein: `shift` (median difference), `spread` (IQR), `bridges` (number of pairs) |
| `.samples` | per plate and bridging sample: `shift` (median over proteins), to find odd bridging samples |
| `.by_group` | with `group`: per protein, the shift in each group and their `difference` |

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
    is [`plot_bridge_pca`](charts.md#plot_bridge_pca).

**Errors.** The same as `bridge_differences`, plus `KeyError` when `group` is not a column.

**Example** (PAD000040, the two raw plate files)

```python
so.bridge_check(clean, group="diabetes", show_ids=False)
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

## `bridge_factors`

Compute the shift per protein and plate. It changes nothing: you get a table you can inspect,
save (`.to_csv`) and later apply with `bridge_plates(..., factors=...)`.

```python
bridge_factors(npx, *, batch="PlateID", reference=None, bridges=None, exclude_warnings=True,
               exclude_below_lod=False, min_bridges=8) -> pandas.DataFrame
```

For each protein and each plate other than the reference:

```text
adjustment = − median( NPX on the plate − NPX on the reference )     over the bridging samples
```

**Inputs.** The [shared arguments](#shared-arguments), plus `min_bridges` (default 8). Proteins
with fewer usable pairs are flagged.

**Output.** One row per plate and protein:

| Column | What it holds |
|---|---|
| `plate`, `reference` | the plate that will be adjusted, and the reference |
| `OlinkID`, `Assay` | the protein |
| `adjustment` | NPX to **add** to this protein on this plate (empty if no bridging pair) |
| `bridges` | number of bridging pairs used |
| `spread` | IQR of the bridge differences: noise that bridging cannot remove |
| `few_bridges` | True when `bridges < min_bridges` |

```python
factors = so.bridge_factors(clean)
factors.sort_values("adjustment").head()
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
| `factors` | DataFrame | `None` | a table from `bridge_factors` to apply; computed here when `None` |
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

## `bridge_differences`

The building block of the other bridging functions. For each bridging sample and protein, it
gives how much higher the other plate reads than the reference.

```python
bridge_differences(npx, *, batch="PlateID", reference=None, bridges=None,
                   exclude_warnings=True, exclude_below_lod=False) -> pandas.DataFrame
```

**Output.** One row per bridging sample, protein and non-reference plate, with these columns:
- `sample` (the tidy ID) and the protein columns;
- `plate`, `reference`, `NPX_reference` and `NPX_plate`;
- `difference` (= plate − reference).

The reference that was used is stored in `df.attrs["reference"]`. The table is empty when no
sample is on two plates.

**Errors**

- `ValueError`: there is no NPX column, or only one plate.
- `KeyError`: `batch` is not a column.
- `ValueError`: `reference` is not one of the plates.
