# Charts

Interactive Plotly charts for every quality and bridging step. Colour-blind safe, with sample IDs
hidden unless you ask for them.
{ .lede }

All charts return a `plotly.graph_objects.Figure`:
- it is interactive in a notebook;
- `fig.write_image("chart.png")` saves a static file;
- `fig.write_html("chart.html")` saves a file you can share.

Colours come from a colour-blind-checked palette in a fixed order. Identity is never shown by
colour alone: legends, marker shapes, and hollow versus filled dots carry it too. Sample IDs
appear in the hover text only with `show_ids=True`. Control samples and control assays are always
left out.

## `plot_below_lod`

Horizontal bars showing the % of measurements below the detection limit, per protein.

```python
plot_below_lod(npx, *, group=None, max_below_lod=None, only_with_below=True) -> Figure
```

| Argument | Default | Meaning |
|---|---|---|
| `npx` | — | table with a `LOD` column |
| `group` | `None` | one bar per group (at most 3), proteins sorted by the difference between groups |
| `max_below_lod` | `None` | draw a threshold line; proteins above it in red (what `apply_qc(max_below_lod=...)` would remove) |
| `only_with_below` | `True` | leave out proteins that are never below LOD |

**Errors.** `ValueError` without a `LOD` column, or with more than 3 groups.

---

## `plot_npx_strip`

Every measurement as a dot, per protein, with groups side by side and the LOD as a dark dash.
Values below LOD are drawn **hollow**.

```python
plot_npx_strip(npx, *, group=None, proteins=None, per_row=30, order="name", show_ids=False, seed=0)
    -> Figure
```

| Argument | Default | Meaning |
|---|---|---|
| `group` | `None` | colour by this column (at most 3 groups) |
| `proteins` | `None` | list of protein names to show; default all |
| `per_row` | `30` | proteins per row; long lists are split over rows |
| `order` | `"name"` | `"name"` (alphabetical) or `"below_lod"` (most below LOD first) |
| `show_ids` | `False` | sample IDs in the hover text |
| `seed` | `0` | the random side-to-side spread of the dots (same seed = same picture) |

**Errors.** `ValueError` when no protein is left to plot, or with more than 3 groups.

---

## `plot_group_by_plate`

Stacked bars showing how many people of each group were measured on each plate. If a group is
missing from a plate, plate differences can look like group differences.

```python
plot_group_by_plate(npx, group) -> Figure
```

**Input.** `group`: a sample-information column. Uses
[`qc_report(...).group_by_plate`](quality.md#qc_report).

---

## `plot_sample_qc`

Each measurement as a dot: its median (x) against its spread (IQR, y) over all proteins. This is
Olink's standard sample check. Dashed lines mark `sd` standard deviations, and outliers are drawn
as large red rings. Colour and marker shape show the plate; there is one panel per Olink panel
(up to 4).

```python
plot_sample_qc(npx, *, sd=3.0, show_ids=False) -> Figure
```

The table behind it: [`sample_qc_table`](quality.md#sample_qc_table).

---

## `plot_bridge_shift`

For each protein, how much higher the other plate reads than the reference, for the same samples:
- grey dots: single bridging samples;
- coloured marker: the median, which is the shift `bridge_plates` removes;
- bar: the middle 50%.

Before bridging the markers sit away from 0; after bridging they sit on 0.

```python
plot_bridge_shift(npx, *, batch="PlateID", reference=None, bridges=None, exclude_warnings=True,
                  exclude_below_lod=False, show_ids=False) -> Figure
```

Arguments as in [bridging](bridging.md#shared-arguments). At most 3 plates besides the reference.

**Errors.** `ValueError` when there are no bridging samples, or more than 3 other plates.

---

## `plot_bridge_pca`

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

---

## `plot_bridge_scatter`

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
