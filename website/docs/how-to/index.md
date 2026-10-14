---
title: "How-to guides for Olink data in Python"
description: "Short recipes for Olink data in Python: read several plates, hide sample IDs, remove outliers, save and reuse bridging factors."
---

# How-to guides

Short recipes for one task each. You know what you want to do; here is how.
{ .lede }

## Ready to use

**Read several plates into one table**

```python
npx = pd.concat([so.read_npx(f) for f in ["plate1.csv", "plate2.csv"]], ignore_index=True)
```

**Read a file with unusual column names**

```python
npx = so.read_npx("export.csv", column_map={"Patient": "SampleID", "Protein": "Assay"})
```

**Make reports you can share without sample IDs**

```python
so.qc_report(data, show_ids=False)
clean = so.apply_qc(data, show_ids=False)
```

**Remove an outlier you found yourself, on one plate only**

```python
clean = so.apply_qc(data, drop_samples={"249": "Inflammation PRAESIIDIUM Giraudi 04_03_2025"})
```

**Check how much the quality step matters (sensitivity check)**

```python
strict  = so.apply_qc(data)
lenient = so.apply_qc(data, drop_warnings=False)
```

**Save the bridging factors and apply them later**

```python
factors = so.bridge_check(clean).factors
factors.to_csv("bridge_factors.csv", index=False)
bridged = so.bridge_plates(clean, factors=pd.read_csv("bridge_factors.csv"))
```

**Bridge whole projects instead of plates**

```python
bridged = so.bridge_plates(both_projects, batch="Project")
```

## Planned

<ul class="soon">
  <li>Read an older wide export (NPX Manager)</li>
  <li>Read Dutch or other European csv files (semicolons, decimal commas)</li>
  <li>Read Explore HT and Reveal parquet files</li>
  <li>Save a chart for a paper or a slide</li>
</ul>
