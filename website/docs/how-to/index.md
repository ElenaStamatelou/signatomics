---
title: "Recipes for Olink data in Python"
description: "Short recipes for Olink data in Python: read several plates, hide sample IDs, remove outliers, save and reuse bridging factors."
---

# Recipes

One task each.
{ .lede }

### Read several plates into one table

```python
files = ["plate1.csv", "plate2.csv"]
npx = pd.concat([so.read_npx(f) for f in files], ignore_index=True)
```

Semicolons, decimal commas and percentages are detected; nothing to set.

### Read a file with unusual column names

```python
npx = so.read_npx("export.csv",
                  column_map={"Patient": "SampleID", "Protein": "Assay"})
```

### Share reports without sample IDs

```python
so.qc_report(data, show_ids=False)
clean = so.apply_qc(data, show_ids=False)
```

### Remove an outlier on one plate

```python
plate2 = "Inflammation PRAESIIDIUM Giraudi 04_03_2025"
clean = so.apply_qc(data, drop_samples={"249": plate2})
```

A list (`drop_samples=["249"]`) removes the sample on every plate. Either way, the removal is
reported.

### Check how much the quality step matters

```python
strict  = so.apply_qc(data)
lenient = so.apply_qc(data, drop_warnings=False)
```

### Save the bridging factors and apply them later

```python
factors = so.bridge_check(clean).factors
factors.to_csv("bridge_factors.csv", index=False)
bridged = so.bridge_plates(clean, factors=pd.read_csv("bridge_factors.csv"))
```

### Bridge whole projects instead of plates

```python
bridged = so.bridge_plates(both_projects, batch="Project")
```

### Save a chart

```python
fig = so.plot_sample_qc(clean)
fig.write_html("sample_qc.html")      # interactive, opens in any browser
```
