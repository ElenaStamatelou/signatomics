"""Quality charts for Olink data (Plotly).

Built on the standard long table from ``read_npx``: detection limits (LOD, also per plate) are
already inside, so no separate LOD tables have to be merged by hand. Every function returns a
Plotly figure — interactive in a notebook, ``fig.write_image("x.png")`` for a static file.

Colours: a colour-blind-checked palette in a fixed order (never cycled); identity is never
carried by colour alone (legends, hollow vs filled markers, labels). Sample IDs appear in hover
text only with ``show_ids=True``.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from .qc import (_below_lod_by_group, _control_masks, _protein_col, _sample_qc_table, _value_col,
                 qc_report)

# Colour-blind-checked categorical palette (fixed order). In scatter/strip charts only the first
# three are used together: those three stay distinguishable for every pair.
SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#4a3aa7", "#e34948"]
CRITICAL = "#d03b3b"
INK, INK_2, MUTED = "#0b0b0b", "#52514e", "#898781"
GRID, AXIS, SURFACE = "#e1e0d9", "#c3c2b7", "#fcfcfb"
FONT = "system-ui, -apple-system, 'Segoe UI', sans-serif"


def _style(fig: go.Figure, title: str, subtitle: str | None = None) -> go.Figure:
    text = title + (f"<br><sup style='color:{INK_2}'>{subtitle}</sup>" if subtitle else "")
    fig.update_layout(
        title=dict(text=text, x=0.01, xanchor="left", font=dict(size=17, color=INK)),
        font=dict(family=FONT, size=12, color=INK_2),
        plot_bgcolor=SURFACE, paper_bgcolor=SURFACE,
        hoverlabel=dict(bgcolor="white", font=dict(family=FONT, color=INK)),
        legend=dict(orientation="h", yanchor="bottom", y=1.0, xanchor="left", x=0,
                    font=dict(color=INK_2), title=dict(font=dict(color=INK_2))),
    )
    fig.update_xaxes(gridcolor=GRID, linecolor=AXIS, zeroline=False, tickfont=dict(color=MUTED))
    fig.update_yaxes(gridcolor=GRID, linecolor=AXIS, zeroline=False, tickfont=dict(color=MUTED))
    return fig


def _samples_only(npx: pd.DataFrame) -> pd.DataFrame:
    ctrl_samples, ctrl_assays = _control_masks(npx)
    return npx[~(ctrl_samples | ctrl_assays)]


def _group_colors(values, max_groups: int = 3) -> dict:
    values = [v for v in pd.unique(pd.Series(values).dropna())]
    if len(values) > max_groups:
        raise ValueError(
            f"{len(values)} groups found; colours stay distinguishable for at most {max_groups} "
            f"groups in this chart. Combine groups or plot a subset.")
    return {v: SERIES[i] for i, v in enumerate(sorted(values, key=str))}


# ------------------------------------------------------------------------------------------
# 1. Share below the detection limit, per protein
# ------------------------------------------------------------------------------------------

def plot_below_lod(npx: pd.DataFrame, *, group: str | None = None, samples=None, on=None,
                   max_below_lod: float | None = None, only_with_below: bool = True) -> go.Figure:
    """Horizontal bars: % of measurements below the detection limit, per protein.

    With ``group``: one bar per group, sorted by the largest difference between groups - a
    protein detectable in one group but not in another stands out at the top.
    With ``max_below_lod``: a dashed threshold line, and proteins above it drawn in red (the
    ones ``apply_qc(max_below_lod=...)`` would remove). By default nothing is marked as removed.
    """
    if group is not None:
        from .qc import _with_groups
        return _plot_below_lod_by_group(_with_groups(npx, samples, on), group, only_with_below)
    rep = qc_report(npx, verbose=False)
    if rep.proteins is None:
        raise ValueError("No LOD column: the share below the detection limit cannot be computed.")
    p = rep.proteins.copy()
    label = "Assay" if "Assay" in p.columns else _protein_col(npx)
    p["pct"] = p["below_lod"] * 100
    never = int((p["pct"] == 0).sum())
    if only_with_below:
        p = p[p["pct"] > 0]
    p = p.sort_values("pct")

    fig = go.Figure()
    if max_below_lod is None:
        fig.add_bar(x=p["pct"], y=p[label].astype(str), orientation="h", name="proteins",
                    showlegend=False, marker=dict(color=SERIES[0], line=dict(color=SURFACE, width=1)),
                    hovertemplate="<b>%{y}</b><br>%{x:.1f}% below LOD<extra></extra>")
    else:
        removed = p["pct"] > max_below_lod * 100
        for flag, name, colour in [(False, f"kept (≤ {max_below_lod:.0%} below LOD)", SERIES[0]),
                                   (True, f"removed by apply_qc (> {max_below_lod:.0%})", CRITICAL)]:
            sub = p[removed == flag]
            fig.add_bar(x=sub["pct"], y=sub[label].astype(str), orientation="h", name=name,
                        marker=dict(color=colour, line=dict(color=SURFACE, width=1)),
                        hovertemplate="<b>%{y}</b><br>%{x:.1f}% below LOD<extra></extra>")
        fig.add_vline(x=max_below_lod * 100, line=dict(dash="dash", color=INK_2, width=1.5))
    subtitle = (f"{len(p)} proteins shown" + (f"; the other {never} are never below LOD"
                                               if only_with_below and never else ""))
    fig.update_layout(barmode="overlay", height=max(320, 22 * len(p) + 160), width=780,
                      margin=dict(l=120, r=30, t=110, b=50), bargap=0.25,
                      yaxis=dict(categoryorder="array", categoryarray=p[label].astype(str).tolist()))
    fig.update_xaxes(title="% of measurements below LOD", range=[0, 102])
    return _style(fig, "Share of measurements below the detection limit", subtitle)


def _plot_below_lod_by_group(npx: pd.DataFrame, group: str, only_with_below: bool) -> go.Figure:
    t = _below_lod_by_group(npx, group)
    groups = [c for c in t.columns if c != "difference"]
    colors = _group_colors(groups)
    if only_with_below:
        t = t[t[groups].max(axis=1) > 0]
    t = t.sort_values("difference")                      # largest difference at the top
    fig = go.Figure()
    for g in sorted(groups, key=str):
        fig.add_bar(x=t[g] * 100, y=t.index.astype(str), orientation="h", name=str(g),
                    marker=dict(color=colors[g], line=dict(color=SURFACE, width=1)),
                    customdata=t["difference"] * 100,
                    hovertemplate=f"<b>%{{y}}</b><br>{g}: %{{x:.0f}}% below LOD"
                                  "<br>difference between groups: %{customdata:.0f} points<extra></extra>")
    fig.update_layout(barmode="group", height=max(360, 30 * len(t) + 170), width=800,
                      margin=dict(l=120, r=30, t=120, b=50), bargap=0.25, bargroupgap=0.05,
                      legend=dict(title=dict(text=f"{group}: ")))
    fig.update_xaxes(title="% of measurements below LOD", range=[0, 102])
    return _style(fig, f"Below the detection limit, per protein and '{group}' group",
                  "sorted by the difference between groups (largest at the top) - "
                  "a protein detectable in only one group can be a signal")


# ------------------------------------------------------------------------------------------
# 2. NPX per protein, by group, with the detection limit
# ------------------------------------------------------------------------------------------

def plot_npx_strip(npx: pd.DataFrame, *, group: str | None = None, samples=None, on=None,
                   proteins=None, per_row: int = 30, order: str = "name", show_ids: bool = False,
                   seed: int = 0) -> go.Figure:
    """Every measurement as a dot, per protein, groups side by side; the LOD as a red line.

    Values below the detection limit are drawn as **hollow** dots, so it is visible how much of
    a protein's signal is unreliable. Long protein lists are split into rows of ``per_row``.

    Parameters
    ----------
    group : sample-information column for colours (at most 3 groups), e.g. ``"diabetes"``.
    proteins : list of proteins (Assay names) to show; default all.
    order : ``"name"`` (alphabetical) or ``"below_lod"`` (most below LOD first).
    show_ids : include sample IDs in the hover text.
    samples, on : your sample metadata and its ID column, when ``group`` is not yet a column of
        ``npx`` (only looked up, as in ``qc_report``).
    """
    from .qc import _with_groups
    df = _samples_only(_with_groups(npx, samples, on))
    value = _value_col(df)
    label = "Assay" if "Assay" in df.columns else _protein_col(df)
    if proteins is not None:
        df = df[df[label].isin(proteins)]
    has_lod = "LOD" in df.columns
    df = df.assign(_below=(df[value] < df["LOD"]) if has_lod else False)
    if order == "below_lod" and has_lod:
        names = df.groupby(label)["_below"].mean().sort_values(ascending=False).index.tolist()
    else:
        names = sorted(df[label].dropna().unique(), key=str)
    if not names:
        raise ValueError("No proteins to plot.")

    groups = [None] if group is None else sorted(df[group].dropna().unique(), key=str)
    colors = {None: SERIES[0]} if group is None else _group_colors(df[group])
    width = 0.8 / len(groups)
    n_rows = int(np.ceil(len(names) / per_row))
    size = int(np.ceil(len(names) / n_rows))            # spread proteins evenly over the rows
    rows = [names[i:i + size] for i in range(0, len(names), size)]
    rng = np.random.default_rng(seed)

    fig = make_subplots(rows=len(rows), cols=1, vertical_spacing=0.16 / max(len(rows), 1))
    for r, chunk in enumerate(rows, start=1):
        pos = {p: i for i, p in enumerate(chunk)}
        sub = df[df[label].isin(chunk)]
        for gi, g in enumerate(groups):
            gsub = sub if g is None else sub[sub[group] == g]
            offset = (gi - (len(groups) - 1) / 2) * width
            x = gsub[label].map(pos) + offset + rng.uniform(-width * 0.35, width * 0.35, len(gsub))
            for below in (False, True):
                part = gsub[gsub["_below"] == below]
                if part.empty:
                    continue
                xs = x[part.index]
                name = ("all samples" if g is None else str(g)) + (" (below LOD)" if below else "")
                hover = (f"<b>%{{customdata[0]}}</b><br>{value}: %{{y:.2f}}"
                         + ("" if g is None else f"<br>{group}: {g}")
                         + ("<br>below LOD" if below else "")
                         + ("<br>sample: %{customdata[1]}" if show_ids else "") + "<extra></extra>")
                custom = np.column_stack([part[label].astype(str),
                                          part["SampleID"].astype(str) if show_ids else [""] * len(part)])
                fig.add_scatter(
                    x=xs, y=part[value], mode="markers", name=name, legendgroup=name,
                    showlegend=(r == 1), customdata=custom, hovertemplate=hover,
                    marker=dict(size=6, color=colors[g] if not below else "rgba(0,0,0,0)",
                                line=dict(color=colors[g], width=1.2), opacity=0.8),
                    row=r, col=1)
        if has_lod:
            lod = sub.groupby(label)["LOD"].agg(["median", "min", "max"]).reindex(chunk)
            fig.add_scatter(
                x=[pos[p] for p in chunk], y=lod["median"], mode="markers",
                name="detection limit (LOD)", legendgroup="lod", showlegend=(r == 1),
                marker=dict(symbol="line-ew", size=26, line=dict(color=CRITICAL, width=2.5)),
                customdata=np.column_stack([chunk, lod["min"], lod["max"]]),
                hovertemplate="<b>%{customdata[0]}</b><br>LOD %{y:.2f}"
                              "<br>(per plate %{customdata[1]:.2f} – %{customdata[2]:.2f})<extra></extra>",
                row=r, col=1)
        fig.update_xaxes(tickvals=list(range(len(chunk))), ticktext=chunk, tickangle=-60,
                         range=[-0.6, len(chunk) - 0.4], row=r, col=1, showgrid=False)
        fig.update_yaxes(title=value, row=r, col=1)

    subtitle = "red line = detection limit (LOD); hollow dots = below the LOD"
    fig.update_layout(height=430 * len(rows) + 90, width=max(700, 34 * min(per_row, len(names)) + 160),
                      margin=dict(l=70, r=30, t=120, b=90))
    return _style(fig, f"{value} per protein" + (f" by {group}" if group else ""), subtitle)


# ------------------------------------------------------------------------------------------
# 3. Groups per plate
# ------------------------------------------------------------------------------------------

def plot_group_by_plate(npx: pd.DataFrame, group: str, *, samples=None, on=None) -> go.Figure:
    """Stacked bars: how many people of each group were measured on each plate.

    A group that is missing from a plate (or strongly over-represented) means plate
    differences can look like group differences.

    samples, on : your sample metadata and its ID column, when ``group`` is not yet a column of
        ``npx`` (only looked up, as in ``qc_report``).
    """
    table = qc_report(npx, group=group, samples=samples, on=on, verbose=False).group_by_plate
    plates = list(table.columns)
    if len(plates) > len(SERIES):
        raise ValueError(f"{len(plates)} plates: too many to colour distinctly.")
    fig = go.Figure()
    for i, plate in enumerate(plates):
        counts = table[plate]
        fig.add_bar(x=[str(g) for g in table.index], y=counts, name=str(plate),
                    marker=dict(color=SERIES[i], line=dict(color=SURFACE, width=2)),
                    text=[str(c) if c else "" for c in counts], textposition="inside",
                    insidetextfont=dict(color="white"),
                    hovertemplate=f"<b>%{{x}}</b><br>{plate}: %{{y}}<extra></extra>")
    fig.update_layout(barmode="stack", height=480, width=760, margin=dict(l=70, r=30, t=150, b=60),
                      legend=dict(title=dict(text="measured on: ")))
    fig.update_yaxes(title="number of people")
    fig.update_xaxes(title=group)
    return _style(fig, f"Where each '{group}' group was measured",
                  "a group missing from a plate → plate differences can look like group differences")


# ------------------------------------------------------------------------------------------
# 4. Per-sample check: median vs spread (Olink's standard sample QC plot)
# ------------------------------------------------------------------------------------------

def plot_sample_qc(npx: pd.DataFrame, *, sd: float = 3.0, show_ids: bool = False) -> go.Figure:
    """Each sample measurement as a dot: median (x) against spread (IQR, y) over all proteins.

    Dashed lines mark ``sd`` standard deviations from the mean; outliers are drawn as large red
    rings. Colour = plate (secondary encoding: marker shape), facets per panel if there are
    several. The numbers behind the chart are in ``qc_report(...).samples``.
    """
    t = _sample_qc_table(npx, sd=sd)
    panels = sorted(t["Panel"].dropna().unique(), key=str) if "Panel" in t.columns else [None]
    panels = panels[:4]
    plates = sorted(t["PlateID"].dropna().unique(), key=str) if "PlateID" in t.columns else [None]
    if len(plates) > len(SERIES):
        raise ValueError(f"{len(plates)} plates: too many to colour distinctly.")
    symbols = ["circle", "square", "diamond", "triangle-up", "cross", "x", "star", "hexagon"]
    several = len(panels) > 1
    fig = make_subplots(rows=1, cols=len(panels), horizontal_spacing=0.08,
                        subplot_titles=[str(p) for p in panels] if several else None)
    for c, panel in enumerate(panels, start=1):
        tp = t if panel is None else t[t["Panel"] == panel]
        for i, plate in enumerate(plates):
            part = tp if plate is None else tp[tp["PlateID"] == plate]
            if part.empty:
                continue
            custom = part["SampleID"].astype(str) if show_ids else pd.Series([""] * len(part))
            fig.add_scatter(
                x=part["median"], y=part["iqr"], mode="markers",
                name=str(plate) if plate is not None else "measurements",
                legendgroup=str(plate), showlegend=(c == 1),
                marker=dict(size=8, color=SERIES[i], symbol=symbols[i],
                            line=dict(color=SURFACE, width=1), opacity=0.85),
                customdata=custom,
                hovertemplate=("median %{x:.2f}<br>IQR %{y:.2f}"
                               + ("<br>sample: %{customdata}" if show_ids else "") + "<extra></extra>"),
                row=1, col=c)
        out = tp[tp["outlier"]]
        if len(out):
            fig.add_scatter(x=out["median"], y=out["iqr"], mode="markers",
                            name=f"outlier (> {sd:g} SD)", legendgroup="outlier", showlegend=(c == 1),
                            marker=dict(size=16, color="rgba(0,0,0,0)", line=dict(color=CRITICAL, width=2)),
                            hoverinfo="skip", row=1, col=c)
        first = tp.iloc[0]
        for col, axis in (("median", "x"), ("iqr", "y")):
            for edge in ("low", "high"):
                v = first[f"{col}_{edge}"]
                kw = dict(line=dict(dash="dash", color=MUTED, width=1), row=1, col=c)
                fig.add_vline(x=v, **kw) if axis == "x" else fig.add_hline(y=v, **kw)
        fig.update_xaxes(title="median over proteins", row=1, col=c)
        fig.update_yaxes(title="IQR over proteins" if c == 1 else None, row=1, col=c)
    n_out = int(t["outlier"].sum())
    panel_note = f"{panels[0]} — " if (not several and panels[0] is not None) else ""
    fig = _style(fig, "Sample check: median vs spread of each measurement",
                 f"{panel_note}{n_out} outlier(s) beyond {sd:g} SD (dashed lines) — "
                 f"check them for technical problems")
    fig.update_layout(height=580, width=max(680, 480 * len(panels)),
                      margin=dict(l=70, r=30, t=90 if not several else 120, b=150),
                      legend=dict(title=dict(text="plate: "), orientation="h",
                                  yanchor="top", y=-0.2, xanchor="left", x=0))
    return fig


# ------------------------------------------------------------------------------------------
# 5. Bridging: plate shift per protein, PCA before/after, bridging-sample scatter
# ------------------------------------------------------------------------------------------

def _plate_style(plates) -> dict:
    symbols = ["circle", "square", "diamond", "triangle-up", "cross", "x", "star", "hexagon"]
    if len(plates) > len(SERIES):
        raise ValueError(f"{len(plates)} plates: too many to colour distinctly.")
    return {p: (SERIES[i], symbols[i]) for i, p in enumerate(plates)}


def plot_bridge_shift(npx: pd.DataFrame, *, batch: str = "PlateID", reference=None,
                      bridges=None, exclude_warnings: bool = True,
                      exclude_below_lod: bool = False, show_ids: bool = False) -> go.Figure:
    """Per protein: how much higher the other plate reads than the reference, for the same samples.

    Each small grey dot is one bridging sample (plate minus reference, in NPX); the coloured
    marker is the median (= the shift that ``bridge_plates`` removes) and the bar the middle 50%
    (IQR). Before bridging the markers sit away from 0; after bridging they sit on 0, and what
    remains is the spread - measurement noise that bridging cannot remove.
    Arguments as in ``bridge_check``.
    """
    from .bridge import _bridge_differences

    d = _bridge_differences(npx, batch=batch, reference=reference, bridges=bridges,
                           exclude_warnings=exclude_warnings, exclude_below_lod=exclude_below_lod)
    if d.empty:
        raise ValueError("No bridging samples found (no sample on the reference and another plate).")
    ref = d.attrs["reference"]
    label = "Assay" if "Assay" in d.columns else _protein_col(d)
    plates = list(pd.unique(d["plate"]))
    if len(plates) > 3:
        raise ValueError(f"{len(plates)} plates besides the reference; plot at most 3 at a time "
                         f"(pass a subset of the data).")
    stats = (d.groupby(["plate", label])["difference"]
             .agg(median="median", q1=lambda s: s.quantile(0.25), q3=lambda s: s.quantile(0.75),
                  n="size").reset_index())
    order = stats.groupby(label)["median"].median().sort_values().index.astype(str).tolist()
    style = _plate_style([ref] + plates)
    offsets = {p: (i - (len(plates) - 1) / 2) * 0.28 for i, p in enumerate(plates)}
    pos = {p: i for i, p in enumerate(order)}
    rng = np.random.default_rng(0)

    fig = go.Figure()
    fig.add_vline(x=0, line=dict(color=INK_2, width=1.5))
    for p in plates:
        dp = d[d["plate"] == p]
        y = dp[label].astype(str).map(pos) + offsets[p] + rng.uniform(-0.12, 0.12, len(dp))
        custom = np.column_stack([dp[label].astype(str),
                                  dp["sample"].astype(str) if show_ids else [""] * len(dp)])
        fig.add_scatter(
            x=dp["difference"], y=y, mode="markers", name="one bridging sample",
            legendgroup="dots", showlegend=(p == plates[0]),
            marker=dict(size=4, color=MUTED, opacity=0.45), customdata=custom,
            hovertemplate="<b>%{customdata[0]}</b><br>difference %{x:+.2f} NPX"
                          + ("<br>sample: %{customdata[1]}" if show_ids else "") + "<extra></extra>")
    for p in plates:
        sp = stats[stats["plate"] == p]
        y = sp[label].astype(str).map(pos) + offsets[p]
        colour, symbol = style[p]
        fig.add_scatter(
            x=sp["median"], y=y, mode="markers", name=f"{p}: median (bar = middle 50%)",
            marker=dict(size=9, color=colour, symbol=symbol, line=dict(color=SURFACE, width=2)),
            error_x=dict(type="data", symmetric=False, array=sp["q3"] - sp["median"],
                         arrayminus=sp["median"] - sp["q1"], color=colour, thickness=2, width=0),
            customdata=np.column_stack([sp[label].astype(str), sp["q1"], sp["q3"], sp["n"]]),
            hovertemplate="<b>%{customdata[0]}</b><br>shift (median) %{x:+.2f} NPX"
                          "<br>middle 50%: %{customdata[1]:+.2f} to %{customdata[2]:+.2f}"
                          "<br>bridging samples: %{customdata[3]}<extra></extra>")
    med = stats.groupby("plate")["median"].median()
    subtitle = ("same samples, other plate minus reference ("
                + "; ".join(f"median {round(med[p], 2) + 0.0:+.2f} NPX" for p in plates)
                + ") - 0 = lined up")
    fig.update_layout(height=max(380, 17 * len(order) + 190), width=820,
                      margin=dict(l=120, r=30, t=130, b=60))
    fig.update_yaxes(tickvals=list(range(len(order))), ticktext=order,
                     range=[-0.8, len(order) - 0.2], showgrid=False)
    fig.update_xaxes(title=f"difference in NPX (vs reference {ref})")
    return _style(fig, "Plate shift per protein, measured on the bridging samples", subtitle)


def _pca_scores(x: pd.DataFrame) -> tuple[np.ndarray, np.ndarray]:
    """First two principal components of a samples x proteins table (median-imputed, scaled)."""
    x = x.loc[:, x.notna().mean() >= 0.5]
    x = x.fillna(x.median())
    sd = x.std(ddof=0).replace(0, np.nan)
    z = ((x - x.mean()) / sd).dropna(axis=1).to_numpy()
    u, s, _ = np.linalg.svd(z, full_matrices=False)
    explained = s ** 2 / (s ** 2).sum()
    return u[:, :2] * s[:2], explained[:2]


def plot_bridge_pca(before: pd.DataFrame, after: pd.DataFrame | None = None, *,
                    batch: str = "PlateID", samples: str = "bridges",
                    exclude_warnings: bool = True, show_ids: bool = False) -> go.Figure:
    """PCA of the bridging samples before (and after) bridging, coloured by plate.

    Each bridging sample appears once per plate; a thin line connects its measurements. Before
    bridging the plates form two clouds and the lines are long; after bridging the clouds should
    overlap and the lines shrink to measurement noise. This is the independent check of bridging:
    the shift was fitted per protein, the PCA looks at all proteins together.

    Parameters
    ----------
    before : table before bridging (several plates in one table).
    after : the same table after ``bridge_plates``; if given, a second panel is drawn.
    samples : ``"bridges"`` (default: only samples on more than one plate) or ``"all"``.
    """
    from .bridge import _sample_key
    from .qc import _is_warning

    if samples not in ("bridges", "all"):
        raise ValueError("samples must be 'bridges' or 'all'")
    tables = [("before bridging", before)] + ([("after bridging", after)] if after is not None else [])
    fig = make_subplots(rows=1, cols=len(tables), horizontal_spacing=0.1,
                        subplot_titles=[t for t, _ in tables])
    notes = []
    for c, (title, data) in enumerate(tables, start=1):
        df = _samples_only(data)
        if exclude_warnings and "QC_Warning" in df.columns:
            df = df[~_is_warning(df["QC_Warning"])]
        df = df.assign(_sample=_sample_key(df)).dropna(subset=["_sample", batch])
        protein = _protein_col(df)
        wide = df.pivot_table(index=["_sample", batch], columns=protein, values="NPX", aggfunc="mean")
        n_plates = wide.reset_index().groupby("_sample")[batch].transform("nunique").to_numpy()
        if samples == "bridges":
            wide = wide[n_plates > 1]
        if len(wide) < 3:
            raise ValueError("Fewer than 3 measurements to plot; are there bridging samples?")
        scores, explained = _pca_scores(wide)
        pts = wide.index.to_frame(index=False).assign(pc1=scores[:, 0], pc2=scores[:, 1])
        plates = sorted(pts[batch].unique(), key=str)
        style = _plate_style(plates)
        on_several = pts.groupby("_sample")[batch].transform("nunique") > 1
        xs, ys = [], []                      # lines between the measurements of the same sample
        for _, g in pts[on_several].groupby("_sample"):
            for i in range(1, len(g)):
                xs += [g["pc1"].iloc[0], g["pc1"].iloc[i], None]
                ys += [g["pc2"].iloc[0], g["pc2"].iloc[i], None]
        if xs:
            fig.add_scatter(x=xs, y=ys, mode="lines", name="same sample", legendgroup="pairs",
                            showlegend=(c == 1), line=dict(color=AXIS, width=1.2),
                            hoverinfo="skip", row=1, col=c)
        for p in plates:
            pp = pts[pts[batch] == p]
            colour, symbol = style[p]
            fig.add_scatter(
                x=pp["pc1"], y=pp["pc2"], mode="markers", name=str(p), legendgroup=str(p),
                showlegend=(c == 1),
                marker=dict(size=9, color=colour, symbol=symbol, line=dict(color=SURFACE, width=1.5)),
                customdata=pp["_sample"].astype(str) if show_ids else None,
                hovertemplate=(f"{p}<br>PC1 %{{x:.1f}}, PC2 %{{y:.1f}}"
                               + ("<br>sample: %{customdata}" if show_ids else "") + "<extra></extra>"),
                row=1, col=c)
        fig.update_xaxes(title=f"PC1 ({explained[0]:.0%} of variation)", row=1, col=c)
        fig.update_yaxes(title=f"PC2 ({explained[1]:.0%})", row=1, col=c)
        paired = pts[on_several].groupby("_sample")[["pc1", "pc2"]].agg(lambda s: s.max() - s.min())
        if len(paired):
            notes.append(f"{title}: a sample's two measurements lie "
                         f"{np.hypot(paired['pc1'], paired['pc2']).median():.1f} apart (median)")
    fig = _style(fig, "Bridging samples in a PCA, coloured by plate",
                 "; ".join(notes) + " - after bridging, the plates should overlap")
    fig.update_layout(height=540, width=520 * len(tables) + 60, margin=dict(l=70, r=30, t=130, b=110),
                      legend=dict(title=dict(text="plate: "), orientation="h", yanchor="top", y=-0.16,
                                  xanchor="left", x=0))
    return fig


def plot_bridge_scatter(npx: pd.DataFrame, proteins, *, batch: str = "PlateID", reference=None,
                        exclude_warnings: bool = True, show_ids: bool = False,
                        cols: int = 3) -> go.Figure:
    """Bridging samples: NPX on the reference plate (x) against the other plate (y), per protein.

    Points on the dashed diagonal (y = x) mean both plates give the same value. Before bridging
    the cloud lies above or below the diagonal by the plate shift; after bridging it lies on it.
    Dotted lines mark the detection limit (LOD) on each plate. 1 to 12 proteins (Assay names).
    """
    from .bridge import _bridge_differences

    proteins = [proteins] if isinstance(proteins, str) else list(proteins)
    if not 1 <= len(proteins) <= 12:
        raise ValueError("Give 1 to 12 proteins.")
    d = _bridge_differences(npx, batch=batch, reference=reference, exclude_warnings=exclude_warnings)
    ref = d.attrs["reference"]
    label = "Assay" if "Assay" in d.columns else _protein_col(d)
    missing = [p for p in proteins if p not in set(d[label])]
    if missing:
        raise ValueError(f"Not found among the bridging-sample proteins: {', '.join(missing)}")
    plates = list(pd.unique(d["plate"]))
    style = _plate_style([ref] + plates)
    lod = npx.groupby([batch, label])["LOD"].median() if "LOD" in npx.columns else None
    ncol = min(cols, len(proteins))
    rows = int(np.ceil(len(proteins) / ncol))
    fig = make_subplots(rows=rows, cols=ncol, subplot_titles=proteins,
                        horizontal_spacing=0.08, vertical_spacing=0.25 / rows)
    for k, prot in enumerate(proteins):
        r, c = k // ncol + 1, k % ncol + 1
        dp = d[d[label] == prot]
        values = dp[["NPX_reference", "NPX_plate"]].to_numpy()
        lo, hi = float(np.nanmin(values)) - 0.3, float(np.nanmax(values)) + 0.3
        fig.add_scatter(x=[lo, hi], y=[lo, hi], mode="lines", name="same value (y = x)",
                        legendgroup="id", showlegend=(k == 0),
                        line=dict(color=INK_2, width=1, dash="dash"), hoverinfo="skip", row=r, col=c)
        for p in plates:
            pp = dp[dp["plate"] == p]
            colour, symbol = style[p]
            fig.add_scatter(
                x=pp["NPX_reference"], y=pp["NPX_plate"], mode="markers", name=f"{p} vs reference",
                legendgroup=str(p), showlegend=(k == 0),
                marker=dict(size=8, color=colour, symbol=symbol, line=dict(color=SURFACE, width=1.5)),
                customdata=pp["sample"].astype(str) if show_ids else None,
                hovertemplate=("reference %{x:.2f}<br>other plate %{y:.2f}"
                               + ("<br>sample: %{customdata}" if show_ids else "") + "<extra></extra>"),
                row=r, col=c)
            if lod is not None and (p, prot) in lod.index:
                fig.add_hline(y=lod[(p, prot)], line=dict(color=colour, width=1, dash="dot"),
                              row=r, col=c)
        if lod is not None and (ref, prot) in lod.index:
            fig.add_vline(x=lod[(ref, prot)], line=dict(color=style[ref][0], width=1, dash="dot"),
                          row=r, col=c)
        fig.update_xaxes(range=[lo, hi], title="reference" if r == rows else None, row=r, col=c)
        fig.update_yaxes(range=[lo, hi], title="other plate" if c == 1 else None, row=r, col=c)
    fig = _style(fig, "Bridging samples: the same tube on two plates",
                 f"x = reference ({ref}); dotted = LOD per plate; on the dashed line = same value")
    fig.update_layout(height=330 * rows + 210, width=330 * ncol + 120,
                      margin=dict(l=70, r=30, t=120, b=120),
                      legend=dict(orientation="h", yanchor="top", y=-0.1 / rows - 0.05,
                                  xanchor="left", x=0))
    return fig
