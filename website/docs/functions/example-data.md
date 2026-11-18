---
title: "Example data: synthetic Olink NPX data to try signatomics"
description: "example_file, example_data and example_samples: a synthetic two-plate Olink Target 96 study with sample metadata, to try every function without data of your own."
---

# Example data

Synthetic Olink data that comes with the package, so every function can be tried without data of
your own.
{ .lede }

## `example_file`

The path to a small, synthetic NPX file that comes with the package. Read it like your own file.

```python
example_file(plate=1) -> str
```

| Argument | Type | Default | Meaning |
|---|---|---|---|
| `plate` | int | `1` | `1` or `2`: the example has two plates, one file each |

```python
data = so.read_npx(so.example_file())        # plate 1: 40 samples + 8 controls, 92 proteins
plate2 = so.read_npx(so.example_file(2))     # plate 2: 32 samples + 8 controls; 12 also on plate 1
```

---

## `example_data`

Both example plates in one table, in the [standard table](index.md#the-standard-table), as
`read_npx` returns them. The values are synthetic and the same every time.

```python
example_data() -> pandas.DataFrame
```

60 people (30 cases, 30 controls), 92 proteins named after the Target 96 Inflammation panel, two
plates: P01–P40 on Plate 1, P29–P60 on Plate 2, so P29–P40 are bridge samples. Each plate also
holds the 8 controls of a Target 96 plate (3 negative, 3 plate and 2 sample controls), marked in
`SampleType`. Plate 2 reads higher, by a different amount per protein. Built in on purpose, so the
checks have something to find:

- Plate 1: QC status WARN for P05 and FAIL for P20; P07 far too low on every protein (an outlier);
- Plate 2: QC status WARN for P41; P33 far too low on every protein (a broken bridge sample);
- proteins mostly below the LOD (such as IL4, IL5 and IL13), as on the real panel;
- IL-17C detectable mostly in cases;
- four proteins higher in cases (IL6, CXCL10, HGF, OSM).

---

## `example_samples`

The sample metadata that go with `example_data` and `example_file`: one row per person, with `sample_id`, `group`
(`"case"` / `"control"`), `age` and `sex`. Synthetic values.

```python
example_samples(plate=None) -> pandas.DataFrame
```

| Argument | Type | Default | Meaning |
|---|---|---|---|
| `plate` | int | `None` | `None`: all 60 people; `1` or `2`: only the people on that plate (P01 to P40, or P29 to P60), to go with `example_file(1)` or `example_file(2)` |

```python
data = so.add_sample_info(so.example_data(), so.example_samples(), on="sample_id")
```
