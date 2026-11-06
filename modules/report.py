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

TOOLS = ["scx", "scConvert", "crosscell", "anndataR", "zellkonverter", "convert2anndata", "easySCF"]
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
UNIT = "case"  # "dataset" while drawing the by-dataset versions (write())


def fmt_s(s):
    return f"{s:.1f} s" if s < 100 else f"{s:.0f} s"


def fmt_mb(mb):
    return f"{mb:.0f} MB" if mb < 1024 else f"{mb / 1024:.1f} GB"


def collect(rows):
    """(dir, tool) -> dict(n, ok, faithful, wall[], rss[], phases{name: (secs[], mb[])})"""
    agg = defaultdict(lambda: {"n": 0, "ok": 0, "faithful": 0, "wall": [], "rss": [], "work": [],
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
        # wall minus start-up: R + package loading for the R tools (obkit `load`);
        # the binaries have none to subtract
        a["work"].append(r["conversion_time_seconds"] - r.get("phases", {}).get("load", [0])[0])
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
    if UNIT == "dataset":
        md = ["# By dataset", "",
              "As above, but one row per dataset: the 4 variants of an rds dataset (v4/v5 × raw/processed)",
              "count once, with the median over its converted variants. Faithful = every variant faithful.",
              "Converted = at least one variant converted. n = 11 rds and 13 h5ad datasets.", ""]
        md += tables_body(agg)
        return "\n".join(md) + "\n"
    md = ["# CrossCell benchmark, reproduced", "",
          "Per tool and direction: cases converted / cases, cases passing every fidelity property,",
          "then wall time and peak RSS (GNU time ru_maxrss) as the median over converted cases of",
          "each case's median of 3 runs; range = min–max over those cases. \"excl. start-up\" subtracts",
          "R + package loading (4–10 s per run), which a user converting inside an R session",
          "doesn't pay. The 42 rds cases are 11 datasets in up to 4 variants (v4/v5 × raw/processed),",
          "so per-dataset n is 11, not 42.", ""]
    return "\n".join(md + tables_body(agg)) + "\n"


def tables_body(agg):
    md = []
    for d, label in DIRS.items():
        md += [f"## {label}", "", "| tool | converted | faithful | wall | excl. start-up | range | peak RSS |",
               "|---|--:|--:|--:|--:|--:|--:|"]
        for t in TOOLS:
            a = agg.get((d, t))
            if not a:
                continue
            w = (f"{fmt_s(st.median(a['wall']))} | {fmt_s(st.median(a['work']))} | "
                 f"{min(a['wall']):.1f}–{max(a['wall']):.1f} | {fmt_mb(st.median(a['rss']))}") if a["wall"] else "| | |"
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
    return md


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
                note = f"{a['ok']}/{a['n']} converted, {a['faithful']} faithful" if key == "wall" else f"{a['ok']}/{a['n']}"
                ax.text(max(vals) * 1.12, y, f"{fmt(med)}  ({note})", va="center", fontsize=8, color=INK2)
        ax.set_xscale("log")
        if key == "wall":
            ax.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:g} s"))
        else:
            ax.xaxis.set_major_locator(matplotlib.ticker.FixedLocator([10, 100, 1024, 10240]))
            ax.xaxis.set_major_formatter(matplotlib.ticker.FixedFormatter(["10 MB", "100 MB", "1 GB", "10 GB"]))
        ax.set_title(title, loc="left", fontsize=11, color=INK, pad=8)
        ax.set_xlim(right=ax.get_xlim()[1] * 6)

    axes[0].set_yticks(range(len(tools)))
    axes[0].set_yticklabels(tools, fontsize=10, color=INK)
    axes[0].invert_yaxis()
    handles = [plt.Rectangle((0, 0), 1, 1, color=c) for c in DIR_COLOR.values()]
    fig.legend(handles, DIRS.values(), loc="upper right", ncol=2, frameon=False, fontsize=9)
    fig.suptitle(f"CrossCell test cases: median over converted {UNIT}s (ticks: each {UNIT})",
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
                v = st.median(mb) / 1024
                ax2.hlines(y + dy, 0, v, color=color, linewidth=1.4, alpha=0.7, zorder=2)  # stem to the value
                ax2.scatter(v, y + dy, s=36, color=color, edgecolors=SURFACE, linewidths=2, zorder=3)
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
    fig.suptitle(f"By phase, median over converted {UNIT}s (bar ends: sum of the phase medians)",
                 x=0.01, ha="left", fontsize=12, color=INK)
    fig.legend(handles, [n for n, _, _ in PHASES if n in used], loc="upper left",
               bbox_to_anchor=(0.01, 0.95), ncol=min(3, len(handles)), frameon=False, fontsize=9)
    fig.tight_layout(rect=(0, 0, 1, 0.9))
    save(fig, path)


# ---- the paper's fig3 (notebooks/fig3_performance.ipynb), from our runs, plus scx ----
TOOL_COLOR = {"scx": "#2a78d6", "scConvert": "#7b5ea7", "crosscell": "#eb6834", "anndataR": "#1baf7a",
              "zellkonverter": "#eda100", "convert2anndata": "#e87ba4", "easySCF": "#008300"}
TOOL_MARKER = {"scx": "o", "scConvert": "P", "crosscell": "o", "zellkonverter": "s", "anndataR": "^",
               "convert2anndata": "D", "easySCF": "v"}
# nnz and v4/v5 x raw/processed per case, from the paper's config. Its nnz match
# scx's own total_nnz on every case scx converted (43); cellxgene_pbmc_15k has
# none in the config, so it takes scx's count.
NNZ_MISSING = {"cellxgene_pbmc_15k": 14341801}


def case_meta():
    import json
    cfg = Path(__file__).resolve().parent.parent / "config" / "benchmark_testcases.json"
    meta = {}
    for cases in json.loads(cfg.read_text())["test_cases"].values():
        for c in cases:
            tid = c["test_id"].replace(".", "_")
            meta[tid] = {"nnz": c.get("nnz") or NNZ_MISSING.get(tid), "dataset": c["dataset"].replace(".", "_"),
                         "cat": f"{c['seurat_version']} {c['type']}" if c.get("seurat_version") else ""}
    for ds in {m["dataset"] for m in list(meta.values())}:  # by-dataset rows: median nnz of the variants
        nn = [m["nnz"] for m in meta.values() if m.get("dataset") == ds and m["nnz"]]
        meta.setdefault(ds, {"nnz": st.median(nn) if nn else None, "dataset": ds, "cat": "dataset"})
    return meta


def by_dataset(rows):
    """One row per (dataset, tool): converted if any variant converted, time / RSS / phases as the
    median over converted variants, faithful only if every variant is."""
    ds = {k: m["dataset"] for k, m in case_meta().items()}
    groups = defaultdict(list)
    for (case, tool), r in rows.items():
        groups[(ds.get(case, case), tool)].append(r)
    out = {}
    for key, rs in groups.items():
        ok = [r for r in rs if r["status"] == "success"]
        r = {**rs[0], "test_id": key[0], "status": "success" if ok else rs[0]["status"],
             "check": "valid" if all(r.get("check") == "valid" for r in rs) else "invalid"}
        if ok:
            r["conversion_time_seconds"] = st.median(x["conversion_time_seconds"] for x in ok)
            r["peak_memory_mb"] = st.median(x["peak_memory_mb"] for x in ok)
            names = {n for x in ok for n in x.get("phases", {})}
            r["phases"] = {n: [st.median(x["phases"][n][i] for x in ok if n in x.get("phases", {})) for i in (0, 1)]
                           for n in names}
        out[key] = r
    return out


def fig_paper3(rows, path):
    from scipy import stats
    meta = case_meta()
    ok = defaultdict(dict)  # (direction, tool) -> {case: perf}
    for (case, tool), r in rows.items():
        if r["status"] == "success" and meta.get(case, {}).get("nnz"):
            ok[(r["direction"], tool)][case] = r
    fig, axes = plt.subplots(2, 2, figsize=(14, 11.5))
    fig.patch.set_facecolor(SURFACE)
    for ax in axes.flat:
        style(ax)
        ax.grid(axis="y", color=GRID, linewidth=0.8)
    ax_a, ax_b, ax_c, ax_d = axes.flat
    r2h = "rds_to_h5ad"

    def fit(ax, xs, ys, color, label=None):
        if len(xs) >= 3:
            f = stats.linregress(xs, ys)
            x = [min(xs), max(xs)]
            ax.plot(x, [f.slope * v + f.intercept for v in x], "--", color=color, linewidth=1.3, alpha=0.8,
                    label=label and f"{label} fit (R²={f.rvalue ** 2:.3f})")

    # A: time vs nnz. crosscell as published (colour = Seurat version x type), plus scx.
    cat_color = {"V4 raw": "#025e8d", "V4 processed": "#f58220", "V5 raw": "#be1818", "V5 processed": "#c75301"}
    cc = ok[(r2h, "crosscell")]
    for case, r in cc.items():
        ax_a.scatter(meta[case]["nnz"] / 1e6, r["conversion_time_seconds"], s=46,
                     color=cat_color.get(meta[case]["cat"], INK2), edgecolors=SURFACE, linewidths=0.8, zorder=3)
    fit(ax_a, [meta[c]["nnz"] / 1e6 for c in cc], [r["conversion_time_seconds"] for r in cc.values()], INK, "crosscell")
    sx = ok[(r2h, "scx")]
    ax_a.scatter([meta[c]["nnz"] / 1e6 for c in sx], [r["conversion_time_seconds"] for r in sx.values()], s=46,
                 facecolors="none", edgecolors=TOOL_COLOR["scx"], linewidths=1.4, zorder=3, label="scx")
    fit(ax_a, [meta[c]["nnz"] / 1e6 for c in sx], [r["conversion_time_seconds"] for r in sx.values()],
        TOOL_COLOR["scx"], "scx")
    cats = [k for k in cat_color if any(meta[c]["cat"] == k for c in cc)]
    handles = [plt.Rectangle((0, 0), 1, 1, color=cat_color[k]) for k in cats]
    h2, l2 = ax_a.get_legend_handles_labels()
    ax_a.legend(handles + h2, [f"crosscell {k}" for k in cats] + l2, loc="upper left", frameon=False, fontsize=8)
    ax_a.set_xlabel("Non-zero elements (millions)", color=INK2)
    ax_a.set_ylabel("Conversion time (s)", color=INK2)
    ax_a.set_title("A   Time vs nnz (rds → h5ad)", loc="left", fontsize=11, color=INK)

    # B: time ratio tool / crosscell, per shared case.
    tools_b = ["scx", "zellkonverter", "anndataR", "convert2anndata", "easySCF"]
    data, labels = [], []
    for t in tools_b:
        tr = ok[(r2h, t)]
        ratios = [tr[c]["conversion_time_seconds"] / cc[c]["conversion_time_seconds"] for c in tr if c in cc]
        if ratios:
            data.append(ratios)
            labels.append(t)
    bp = ax_b.boxplot(data, widths=0.5, patch_artist=True, showfliers=False)
    for i, (patch, t) in enumerate(zip(bp["boxes"], labels)):
        patch.set_facecolor(TOOL_COLOR[t]); patch.set_alpha(0.35); patch.set_edgecolor(INK2)
        xs = [i + 1 + 0.08 * ((j * 7919) % 11 - 5) / 5 for j in range(len(data[i]))]  # fixed jitter
        ax_b.scatter(xs, data[i], s=16, color=TOOL_COLOR[t], alpha=0.8, edgecolors=SURFACE, linewidths=0.4, zorder=3)
        ax_b.text(i + 1, max(data[i]) * 1.08, f"{st.median(data[i]):.1f}×\n(n={len(data[i])})",
                  ha="center", fontsize=8, color=INK)
    for m in bp["medians"]:
        m.set_color(INK)
    ax_b.axhline(1, color=INK2, linestyle="--", linewidth=1, alpha=0.6)
    ax_b.set_yscale("log")
    ax_b.set_ylim(top=max(max(d) for d in data) * 3)  # room for the median labels
    ax_b.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:g}×"))
    ax_b.set_xticks(range(1, len(labels) + 1))
    ax_b.set_xticklabels(labels, fontsize=9, color=INK)
    ax_b.set_ylabel("Time ratio (tool / crosscell), log scale", color=INK2)
    ax_b.set_title("B   Time ratio vs crosscell (rds → h5ad)", loc="left", fontsize=11, color=INK)

    # C: peak memory vs nnz, every tool, linear fits, MB per million nnz.
    eff = []
    for t in ["scx", "crosscell", "zellkonverter", "anndataR", "convert2anndata", "easySCF"]:
        tr = ok[(r2h, t)]
        if not tr:
            continue
        xs = [meta[c]["nnz"] / 1e6 for c in tr]
        ys = [r["peak_memory_mb"] for r in tr.values()]
        ax_c.scatter(xs, ys, s=34, color=TOOL_COLOR[t], marker=TOOL_MARKER[t], alpha=0.75,
                     edgecolors=SURFACE, linewidths=0.5, label=t, zorder=3)
        fit(ax_c, xs, ys, TOOL_COLOR[t])
        eff.append(f"{t:>15}: {st.median(y / x for x, y in zip(xs, ys)):.0f}")
    ax_c.text(0.98, 0.02, "MB / M nnz (median)\n" + "\n".join(eff), transform=ax_c.transAxes, ha="right",
              va="bottom", fontsize=8, family="monospace", color=INK2)
    ax_c.legend(loc="upper left", frameon=False, fontsize=8)
    ax_c.set_xlabel("Non-zero elements (millions)", color=INK2)
    ax_c.set_ylabel("Peak memory (MB)", color=INK2)
    ax_c.set_title("C   Peak memory vs nnz (rds → h5ad)", loc="left", fontsize=11, color=INK)

    # D: h5ad -> Seurat, paired times against zellkonverter (log-log).
    zk = ok[("h5ad_to_rds", "zellkonverter")]
    pts = []
    for t in ["crosscell", "scx", "scConvert"]:
        tr = ok[("h5ad_to_rds", t)]
        pr = [(tr[c]["conversion_time_seconds"], zk[c]["conversion_time_seconds"]) for c in tr if c in zk]
        if not pr:
            continue
        pts += [v for p in pr for v in p]
        faster = sum(a < b for a, b in pr)
        ax_d.scatter([a for a, _ in pr], [b for _, b in pr], s=50, color=TOOL_COLOR[t], edgecolors=SURFACE,
                     linewidths=1, zorder=3,
                     label=f"{t}: faster in {faster}/{len(pr)}, median {st.median(b / a for a, b in pr):.1f}×")
    if pts:
        lo, hi = min(pts) * 0.5, max(pts) * 2
        ax_d.plot([lo, hi], [lo, hi], "--", color=INK2, linewidth=1, alpha=0.6)
        ax_d.set_xlim(lo, hi); ax_d.set_ylim(lo, hi)
    ax_d.set_xscale("log"); ax_d.set_yscale("log")
    ax_d.grid(color=GRID, linewidth=0.8, which="major")
    ax_d.legend(loc="upper left", frameon=False, fontsize=8)
    ax_d.set_xlabel("Tool time (s)", color=INK2)
    ax_d.set_ylabel("zellkonverter time (s)", color=INK2)
    ax_d.set_title("D   h5ad → Seurat-readable: paired with zellkonverter", loc="left", fontsize=11, color=INK)

    fig.suptitle(f"The paper's fig3 from this run, scx added (one point per {UNIT}; only converted {UNIT}s)",
                 x=0.01, ha="left", fontsize=12, color=INK)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    save(fig, path)


def write(rows, paper, out):
    out = Path(out)
    agg = collect(rows)
    (out / "tables.md").write_text(tables(agg) + "\n" + "\n".join(paper_tables(rows, paper, agg)))
    write_figures(rows, out)


def write_figures(rows, out):
    """Figures, per case and by dataset, plus tables_by_dataset.md: everything results.tsv can
    rebuild (report/report.qmd calls this; the paper's tables need paper_fidelity JSONs)."""
    global UNIT
    agg = collect(rows)
    fig_tools(agg, out / "tools")
    fig_phases(agg, out / "phases")
    fig_paper3(rows, out / "fig3")
    # The same, one row per dataset (pseudoreplication: 42 rds cases are 11 datasets).
    UNIT, ds_rows = "dataset", by_dataset(rows)
    try:
        agg = collect(ds_rows)
        (out / "tables_by_dataset.md").write_text(tables(agg))
        fig_tools(agg, out / "tools_by_dataset")
        fig_phases(agg, out / "phases_by_dataset")
        fig_paper3(ds_rows, out / "fig3_by_dataset")
    finally:
        UNIT = "case"
