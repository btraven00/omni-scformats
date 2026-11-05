"""Figures and tables for the SUMMARY gather, in the style of scx-bench-omni's
make_figures.py (same palette, layout, formatters). One benchmark here has 55
datasets, not one, so per tool and direction: bar = median over the cases the
tool converted of each case's median-of-3; ticks = each case. The fidelity
score sits under each tool's name.

    write(rows, paper, out_dir)
      rows: {(test_id, tool): perf.json + check fields}; paper: paper_fidelity.py rows
    -> tables.md, tools.{png,svg}, phases.{png,svg}

`type_test` (the paper's synthetic type object) is left out of the case counts
and only feeds the paper's type-retention table, as in the paper.
"""
import statistics as st
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker

TOOLS = ["scx", "crosscell", "anndataR", "zellkonverter", "convert2anndata", "easySCF"]
DIRS = {"rds_to_h5ad": "Seurat .rds → h5ad", "h5ad_to_rds": "h5ad → Seurat-readable"}
SURFACE, INK, INK2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e4e3df"
DIR_COLOR = {"rds_to_h5ad": "#2a78d6", "h5ad_to_rds": "#eb6834"}
PHASES = [  # (display name, convert.py phase names, colour)
    ("start-up (R + packages / binary)", ["load"], "#2a78d6"),
    ("read", ["read", "read_fallback"], "#eb6834"),
    ("convert", ["convert", "convert_fallback"], "#eda100"),
    ("matrix stream (scx: read, convert, write in one pass)", ["stream"], "#e87ba4"),
    ("write (saveRDS / h5ad; scx: finalize, flush, hashes)", ["write"], "#008300"),
]
PHASE_OF = {raw: name for name, raws, _ in PHASES for raw in raws}


def fmt_s(s):
    return f"{s:.1f} s" if s < 100 else f"{s:.0f} s"


def fmt_mb(mb):
    return f"{mb:.0f} MB" if mb < 1024 else f"{mb / 1024:.1f} GB"


def collect(rows):
    """(dir, tool) -> dict(n, ok, faithful, wall[], rss[], phases{name: (secs[], mb[])})"""
    agg = defaultdict(lambda: {"n": 0, "ok": 0, "faithful": 0, "wall": [], "rss": [],
                               "phases": defaultdict(lambda: ([], []))})
    for (case, _), r in rows.items():
        if case == "type_test":
            continue
        a = agg[(r["direction"], r["tool"])]
        a["n"] += 1
        a["faithful"] += r.get("check") == "valid"
        if r["status"] != "success":
            continue
        a["ok"] += 1
        a["wall"].append(r["conversion_time_seconds"])
        a["rss"].append(r["peak_memory_mb"])
        per = defaultdict(lambda: [0.0, 0.0])
        for raw, (secs, mb) in r.get("phases", {}).items():
            p = per[PHASE_OF.get(raw, raw)]
            p[0] += secs
            p[1] = max(p[1], mb)
        for name, (secs, mb) in per.items():
            a["phases"][name][0].append(secs)
            if mb > 0:  # 0 = no denet sample fell inside a phase shorter than 50 ms
                a["phases"][name][1].append(mb)
    return agg


# The paper's type test: expected class per column, in its notebook's names.
TYPE_COLS = {"int_col": "int", "dbl_col": "float", "bool_col": "bool", "str_col": "str", "cat_col": "category"}
OURS = {"int": "int", "float": "float", "bool": "bool", "string": "str", "categorical": "category"}


def type_retention(rows, tool):
    """The paper's type check on our type_test verdicts: only the class counts."""
    r = rows.get(("type_test", tool))
    if not r or r.get("check") in (None, "not_run", "error") or not r.get("properties"):
        return None
    props = {p["property"]: p for p in r["properties"]}
    got = {}
    for col, want in TYPE_COLS.items():
        p = props.get(f"obs/{col}")
        if p is None or "dropped" in p["detail"]:
            got[col] = "MISSING"
        elif p["status"] == "fail" and p["detail"].startswith("type "):
            got[col] = OURS.get(p["detail"].split("-> ")[-1].split()[0], "other")
        else:
            got[col] = want
    return got


def paper_tables(rows, paper, agg):
    md = ["## The paper's own checks, as published", "",
          "Pass rates as in its README: converted = every run exited 0 (no output check).", "",
          "| Tool | rds_to_h5ad | h5ad_to_rds |", "|---|--:|--:|"]
    for t in TOOLS:
        c = [f"{agg[(d, t)]['ok']}/{agg[(d, t)]['n']}" if (d, t) in agg else "N/A" for d in DIRS]
        if c != ["N/A", "N/A"]:
            md.append(f"| {t} | " + " | ".join(c) + " |")
    md += ["", "Matrix fidelity (compute_fidelity.ipynb): each tool's h5ad X against crosscell's h5ad X,",
           "matched by position, means over cases; crosscell's row is its own h5ad → rds → h5ad",
           "round trip. Cells and genes are not matched by name and the reference is a tool, not",
           "the input, so this cannot see renamed cells or a wrong reference.", "",
           "| tool | cases | Pearson r | Spearman ρ | MSE | exact match | not computed |",
           "|---|--:|--:|--:|--:|--:|---|"]
    for t in TOOLS:
        rs = [r for r in paper if r["tool"] == t and r["test_id"] != "type_test"]
        if not rs:
            continue
        ok = [r for r in rs if r["status"] == "success"]
        other = sorted({r["status"] for r in rs if r["status"] != "success"})
        m = lambda k, f: f.format(st.mean(r[k] for r in ok)) if ok else ""
        md.append(f"| {t} | {len(ok)}/{len(rs)} | {m('pearson_r', '{:.6f}')} | {m('spearman_rho', '{:.6f}')} | "
                  f"{m('mse', '{:.2e}')} | {m('exact_match_pct', '{:.1f}%')} | {', '.join(other)} |")
    md += ["", "Type retention (compute_type_fidelity.ipynb): the 50-cell object with one column per",
           "type, converted to h5ad; the class read back per column (values not checked).", "",
           "| tool | retention | " + " | ".join(TYPE_COLS) + " |", "|---|--:|" + "---|" * len(TYPE_COLS)]
    for t in TOOLS:
        if ("type_test", t) not in rows:
            continue
        got = type_retention(rows, t)
        if got is None:
            md.append(f"| {t} | 0% | " + " | ".join(["READ_ERROR"] * len(TYPE_COLS)) + " |")
            continue
        n = sum(got[c] == w for c, w in TYPE_COLS.items())
        md.append(f"| {t} | {100 * n / len(TYPE_COLS):.0f}% | " + " | ".join(
            got[c] if got[c] == w else f"{w}→{got[c]} ✖" for c, w in TYPE_COLS.items()) + " |")
    return md + [""]


def tables(agg):
    md = ["# CrossCell benchmark, reproduced", "",
          "Per tool and direction: cases converted / cases, cases passing every fidelity property,",
          "then wall time and peak RSS (GNU time ru_maxrss) as the median over converted cases of",
          "each case's median of 3 runs; range = min–max over those cases.", ""]
    for d, label in DIRS.items():
        md += [f"## {label}", "", "| tool | converted | faithful | wall | range | peak RSS |",
               "|---|--:|--:|--:|--:|--:|"]
        for t in TOOLS:
            a = agg.get((d, t))
            if not a:
                continue
            w = (f"{fmt_s(st.median(a['wall']))} | {min(a['wall']):.1f}–{max(a['wall']):.1f} | "
                 f"{fmt_mb(st.median(a['rss']))}") if a["wall"] else "| |"
            md.append(f"| {t} | {a['ok']}/{a['n']} | {a['faithful']}/{a['n']} | {w} |")
        md.append("")
    md += ["## By phase", "",
           "Median over converted cases of each phase's seconds, and of its denet-sampled peak.",
           "R tools: obkit anchors in each step. scx: its own log timestamps. crosscell: its",
           "printed step durations laid end to end (long steps are rounded to whole seconds).", "",
           "| tool | direction | " + " | ".join(n for n, _, _ in PHASES) + " |",
           "|---|---|" + "--:|" * len(PHASES)]
    for t in TOOLS:
        for d in DIRS:
            a = agg.get((d, t))
            if not a or not a["phases"]:
                continue
            cells = [(f"{fmt_s(st.median(s))} / {fmt_mb(st.median(m))}" if m else fmt_s(st.median(s))) if s else ""
                     for s, m in (a["phases"].get(n, ([], [])) for n, _, _ in PHASES)]
            md.append(f"| {t} | {d} | " + " | ".join(cells) + " |")
    return "\n".join(md) + "\n"


def style(ax):
    ax.set_facecolor(SURFACE)
    for s in ["top", "right", "left"]:
        ax.spines[s].set_visible(False)
    ax.spines["bottom"].set_color(INK2)
    ax.tick_params(colors=INK2, labelsize=9, length=0)
    ax.grid(axis="x", color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)


def save(fig, path):
    for ext in ["png", "svg"]:
        fig.savefig(path.with_suffix(f".{ext}"), dpi=200, facecolor=SURFACE)
    plt.close(fig)


def fig_tools(agg, path):
    tools = [t for t in TOOLS if any((d, t) in agg for d in DIRS)]
    fig, axes = plt.subplots(1, 2, figsize=(11, 0.62 * len(tools) + 1.6), sharey=True)
    fig.patch.set_facecolor(SURFACE)
    h = 0.36
    for ax, key, title, fmt in [(axes[0], "wall", "Wall time (log scale)", fmt_s),
                                (axes[1], "rss", "Peak RSS (log scale)", fmt_mb)]:
        style(ax)
        for row, t in enumerate(tools):
            for j, d in enumerate(DIRS):
                y = row + (j - 0.5) * (h + 0.04)
                a = agg.get((d, t))
                if not a or not a[key]:
                    note = "not supported" if not a else f"failed {a['n']}/{a['n']}"
                    ax.text(0.01, y, note, va="center", fontsize=8, color=INK2,
                            transform=ax.get_yaxis_transform())
                    continue
                vals = a[key]
                med = st.median(vals)
                ax.barh(y, med, height=h, color=DIR_COLOR[d], linewidth=0)
                ax.scatter(vals, [y] * len(vals), marker="|", s=60, color=INK, linewidths=1, zorder=3)
                ax.text(max(vals) * 1.12, y, f"{fmt(med)}  ({a['ok']}/{a['n']})", va="center",
                        fontsize=8, color=INK2)
        ax.set_xscale("log")
        if key == "wall":
            ax.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:g} s"))
        else:
            ax.xaxis.set_major_locator(matplotlib.ticker.FixedLocator([10, 100, 1024, 10240]))
            ax.xaxis.set_major_formatter(matplotlib.ticker.FixedFormatter(["10 MB", "100 MB", "1 GB", "10 GB"]))
        ax.set_title(title, loc="left", fontsize=11, color=INK, pad=8)
        ax.set_xlim(right=ax.get_xlim()[1] * 6)

    def label(t):
        f = [f"{agg[(d, t)]['faithful']}/{agg[(d, t)]['n']}" for d in DIRS if (d, t) in agg]
        return f"{t}\nfaithful {' · '.join(f)}"
    axes[0].set_yticks(range(len(tools)))
    axes[0].set_yticklabels([label(t) for t in tools], fontsize=9, color=INK)
    axes[0].invert_yaxis()
    handles = [plt.Rectangle((0, 0), 1, 1, color=c) for c in DIR_COLOR.values()]
    fig.legend(handles, DIRS.values(), loc="upper right", ncol=2, frameon=False, fontsize=9)
    fig.suptitle("CrossCell test cases: median over converted cases (ticks: each case)",
                 x=0.01, ha="left", fontsize=12, color=INK)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    save(fig, path)


def fig_phases(agg, path):
    rows = [(t, d) for t in TOOLS for d in DIRS if (d, t) in agg and agg[(d, t)]["phases"]]
    if not rows:
        return
    fig, axes = plt.subplots(1, 2, figsize=(12, 0.42 * len(rows) + 1.8), sharey=True,
                             gridspec_kw={"width_ratios": [1.6, 1]})
    fig.patch.set_facecolor(SURFACE)
    ax, ax2 = axes
    for a in axes:
        style(a)
    used = []
    # widest bar, so a segment is labelled only where its number fits
    xmax = max(sum(st.median(v[0]) for v in agg[(d, t)]["phases"].values() if v[0]) for t, d in rows)
    for y, (t, d) in enumerate(rows):
        left, ph = 0.0, agg[(d, t)]["phases"]
        for k, (pname, _, color) in enumerate(PHASES):
            secs, mb = ph.get(pname, ([], []))
            if not secs:
                continue
            w = st.median(secs)
            ax.barh(y, w, left=left, height=0.62, color=color, edgecolor=SURFACE, linewidth=2)
            if w > 0.06 * xmax:
                ax.text(left + w / 2, y, f"{w:.1f}", ha="center", va="center", fontsize=8, color="#ffffff")
            if mb:
                dy = (k - (len(PHASES) - 1) / 2) * 0.14
                ax2.scatter(st.median(mb) / 1024, y + dy, s=36, color=color, edgecolors=SURFACE,
                            linewidths=2, zorder=3)
            left += w
            used.append(pname)
        ax.text(left * 1.02 + 0.1, y, fmt_s(left), va="center", fontsize=8, color=INK2)
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([f"{t}  ·  {DIRS[d]}" for t, d in rows], fontsize=9, color=INK)
    ax.invert_yaxis()
    ax.set_xlim(right=ax.get_xlim()[1] * 1.12)
    ax.set_title("Wall time by phase (s)", loc="left", fontsize=11, color=INK, pad=8)
    ax2.set_title("Peak RSS within each phase (GB)", loc="left", fontsize=11, color=INK, pad=8)
    ax2.set_xlim(left=0)
    handles = [plt.Rectangle((0, 0), 1, 1, color=c) for n, _, c in PHASES if n in used]
    fig.suptitle("By phase, median over converted cases (bar ends: sum of the phase medians)",
                 x=0.01, ha="left", fontsize=12, color=INK)
    fig.legend(handles, [n for n, _, _ in PHASES if n in used], loc="upper left",
               bbox_to_anchor=(0.01, 0.95), ncol=min(3, len(handles)), frameon=False, fontsize=9)
    fig.tight_layout(rect=(0, 0, 1, 0.9))
    save(fig, path)


def write(rows, paper, out):
    out = Path(out)
    agg = collect(rows)
    (out / "tables.md").write_text(tables(agg) + "\n" + "\n".join(paper_tables(rows, paper, agg)))
    fig_tools(agg, out / "tools")
    fig_phases(agg, out / "phases")
