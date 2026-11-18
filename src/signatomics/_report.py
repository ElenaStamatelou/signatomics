"""One printed style for every report: a title, then short ``label: value`` lines.

Rules of the style:

- the title says what was done, e.g. ``Group comparison: case vs control (Welch's t-test)``;
- one fact per line, ``label: value``, with the values lined up;
- a list of results goes underneath its line, indented, as a small aligned table;
- ``WARNING`` lines stay as they are.
"""

from __future__ import annotations

_MAX_LABEL = 16          # longer "labels" are sentences: they are not lined up


def _label(line: str) -> str | None:
    if line.startswith((" ", "WARNING")) or ": " not in line:
        return None
    label = line.split(": ", 1)[0]
    return label if len(label) <= _MAX_LABEL else None


def format_report(title: str, lines: list[str]) -> str:
    """The report as printed text."""
    labels = [lab for lab in map(_label, lines) if lab]
    width = max((len(lab) for lab in labels), default=0) + 1
    subs = [lab for lab in (_label(line.strip()) for line in lines if line.startswith(" ")) if lab]
    sub_width = max((len(lab) for lab in subs), default=0) + 1
    out = [title]
    for line in lines:
        if line.startswith(" "):
            lab = _label(line.strip())
            out.append("    " + (f"{lab + ':':<{sub_width}}  {line.strip().split(': ', 1)[1]}" if lab
                                 else line.strip()))
            continue
        lab = _label(line)
        out.append(f"  {lab + ':':<{width}}  {line.split(': ', 1)[1]}" if lab else "  " + line)
    return "\n".join(out)


def table_lines(rows: list[list[str]]) -> list[str]:
    """Rows of cells as sub-lines of a report, every column lined up (text left, numbers right)."""
    if not rows:
        return []
    widths = [max(len(r[i]) for r in rows) for i in range(len(rows[0]))]

    def cell(text: str, i: int) -> str:
        numeric = text[:1] in "+-0123456789"
        return text.rjust(widths[i]) if numeric and i else text.ljust(widths[i])

    return ["  " + "   ".join(cell(c, i) for i, c in enumerate(r)).rstrip() for r in rows]
