#!/usr/bin/env python3
"""Convert one test case with one tool, `runs` times, measured.

Port of crosscell's benchmark/scripts/run_benchmark.py (git a910b12). The tool
commands are its cmd_* with three changes, none to what a tool computes:
  - each R step is wrapped in obkit phase anchors (load/read/convert/write),
    so time-to-object can be separated from .rds/.h5ad write (gzip saveRDS is
    most of an h5ad -> rds run; see scx docs/rds-write.md);
  - zellkonverter's reticulate fallback uses this env's python, not /opt/venv;
  - library() calls are a timed `load` phase rather than untimed noise.
One arm is new: scx, the binary alone, both directions. h5ad -> Seurat
writes scx's default .h5seurat (BPCells-packed counts + meta.data), a file R
reads with BPCells; no R step, no .rds (scx docs/rds-write.md).

The TO_SEURAT output is a symlink to what the tool wrote (<tool>.rds or
<tool>.h5seurat), so every tool keeps its natural format and its command.

Measurement: GNU time (exact ru_maxrss, the paper's number) around denet
(50 ms samples, children included), as ../omni-scx/modules/measure.sh. Same
rule as the original: the first failing run stops the case. Median of runs.

A failure exits 0 with an empty output: the pass rate is the result, so a
failed conversion is a row in SUMMARY, not a broken DAG.
"""
import json, os, re, shutil, statistics, subprocess, sys, time
from datetime import datetime
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from ob import parse

# Phase anchors (obkit wire format 0.1). phase() returns its value, so
# `x <- phase("read", readRDS(f))` times exactly the read.
PRE = """obkit::logger_init(Sys.getenv("OBKIT_DIR"))
phase <- function(e, x) { obkit::logger_emit(e, "start"); force(x); obkit::logger_emit(e, "end"); invisible(x) }
"""


def R(body):
    return ["Rscript", "-e", PRE + body.strip()]


def cmd_crosscell(i, o, d):
    fmt = "anndata" if d == "rds_to_h5ad" else "seurat"
    return [os.environ.get("CROSSCELL_BIN", "crosscell"), "convert", "-i", i, "-o", o, "-f", fmt]


def cmd_scx(i, o, d):
    # Same flags on every case, as every other tool. Format follows the extension.
    return [os.environ.get("SCX_BIN", "scx"), "convert", i, o]


def cmd_scConvert(i, o, d):
    # scConvert's C CLI (github.com/mianaz/scConvert, ecab18f): h5ad -> h5seurat
    # with no R, as in scx-bench-omni. Its .rds routes need its R package: not here.
    if d != "h5ad_to_rds":
        return None
    return [os.environ.get("SCCONVERT_BIN", "scconvert"), i, o, "--overwrite"]


NATIVE = {"scx": ".h5seurat", "scConvert": ".h5seurat"}  # h5ad -> Seurat: what each tool writes; default .rds


def cmd_zellkonverter(i, o, d):
    if d == "rds_to_h5ad":
        return R(f"""
phase("load", {{ library(Seurat); library(zellkonverter); library(SingleCellExperiment) }})
obj <- phase("read", readRDS("{i}"))
sce <- phase("convert", {{ obj <- UpdateSeuratObject(obj); as.SingleCellExperiment(obj) }})
phase("write", writeH5AD(sce, "{o}"))
""")
    return R(f"""
phase("load", {{ library(zellkonverter); library(Seurat); library(SingleCellExperiment) }})

fix_dimnames <- function(sce) {{
  if (is.null(rownames(sce)) || length(rownames(sce)) == 0) {{
    if ("_index" %in% names(rowData(sce))) {{
      rownames(sce) <- rowData(sce)[["_index"]]
    }} else if ("gene_ids" %in% names(rowData(sce))) {{
      rownames(sce) <- rowData(sce)[["gene_ids"]]
    }} else {{
      rownames(sce) <- paste0("Gene_", seq_len(nrow(sce)))
    }}
  }}
  if (is.null(colnames(sce)) || length(colnames(sce)) == 0) {{
    if ("_index" %in% names(colData(sce))) {{
      colnames(sce) <- colData(sce)[["_index"]]
    }} else if ("obs_names" %in% names(colData(sce))) {{
      colnames(sce) <- colData(sce)[["obs_names"]]
    }} else {{
      colnames(sce) <- paste0("Cell_", seq_len(ncol(sce)))
    }}
  }}
  sce
}}

convert_ok <- FALSE
tryCatch({{
  sce <- phase("read", readH5AD("{i}", reader = "R"))
  obj <- phase("convert", {{ sce <- fix_dimnames(sce); as.Seurat(sce, counts = "X", data = NULL) }})
  phase("write", saveRDS(obj, "{o}"))
  convert_ok <- TRUE
}}, error = function(e) {{
  cat("R reader failed:", conditionMessage(e), "\\n")
}})

if (!convert_ok) {{
  cat("Falling back to reticulate Python reader...\\n")
  Sys.setenv(RETICULATE_PYTHON = "{sys.executable}")
  library(reticulate)
  use_python("{sys.executable}", required = TRUE)
  ad <- import("anndata")
  adata <- phase("read_fallback", ad$read_h5ad("{i}"))
  obj <- phase("convert_fallback", {{ sce <- fix_dimnames(AnnData2SCE(adata)); as.Seurat(sce, counts = "X", data = NULL) }})
  phase("write", saveRDS(obj, "{o}"))
}}
""")


def cmd_anndataR(i, o, d):
    if d == "rds_to_h5ad":
        return R(f"""
phase("load", {{ library(anndataR); library(Seurat) }})
seurat_obj <- phase("read", readRDS("{i}"))
ad <- phase("convert", {{ seurat_obj <- UpdateSeuratObject(seurat_obj); as_AnnData(seurat_obj) }})
phase("write", write_h5ad(ad, "{o}"))
""")
    return R(f"""
phase("load", {{ library(anndataR); library(Seurat) }})
ad <- phase("read", read_h5ad("{i}"))
seurat_obj <- phase("convert", ad$as_Seurat())
phase("write", saveRDS(seurat_obj, "{o}"))
""")


def cmd_convert2anndata(i, o, d):
    if d != "rds_to_h5ad":
        return None
    return R(f"""
phase("load", {{ library(convert2anndata); library(Seurat) }})
obj <- phase("read", readRDS("{i}"))
ad <- phase("convert", {{ obj <- UpdateSeuratObject(obj); sce <- convert_seurat_to_sce(obj); convert_to_anndata(sce) }})
phase("write", anndata::write_h5ad(ad, "{o}"))
""")


def cmd_easySCF(i, o, d):
    # easySCF writes its own H5 layout, not standard H5AD (original's note).
    if d != "rds_to_h5ad":
        return None
    return R(f"""
phase("load", {{ library(easySCFr); library(Seurat) }})
obj <- phase("read", readRDS("{i}"))
obj <- phase("convert", UpdateSeuratObject(obj))
phase("write", saveH5(obj, "{o}"))
""")


def ts(s):
    return datetime.fromisoformat(s.replace("Z", "+00:00")).timestamp()


def obkit_windows(obkit_dir):
    """[(event, t0, t1)] from obkit anchors (R steps). Pairs by (pid, event) in
    file order; an unmatched start (a step that raised) is dropped."""
    f = Path(obkit_dir) / "obkit-events.jsonl"
    open_, out = {}, []
    for line in f.read_text().splitlines() if f.exists() else []:
        try:
            e = json.loads(line)
        except json.JSONDecodeError:
            continue
        k = (e.get("pid"), e["event"])
        if e["phase"] == "start":
            open_[k] = ts(e["ts"])
        elif k in open_:
            out.append((e["event"], open_.pop(k), ts(e["ts"])))
    return out


def scx_windows(log, start, end):
    """scx's own INFO timestamps. "starting conversion" is logged after the input
    is opened (for .rds, the whole file parsed), so process start -> metadata/
    obsm read is `read`; then the matrix stream (read, convert, write in one
    pass); then finalize, output hash and sidecar up to process exit = `write`."""
    at = {}
    for line in log.splitlines():
        m = re.match(r"(?:\x1b\[[0-9;]*m)*(\d{4}-\d\d-\d\dT[\d:.]+Z)", line)
        for key, marker in (("start", "starting conversion"), ("meta", "metadata loaded"),
                            ("aux", "read obsm, varm"), ("matrix", "matrix streamed"), ("end", "conversion complete")):
            if m and marker in line:
                at[key] = ts(m.group(1))
    if not {"start", "matrix", "end"} <= at.keys():
        return []
    read_end = at.get("aux", at.get("meta", at["start"]))
    return [("read", start, read_end), ("stream", read_end, at["matrix"]), ("write", at["matrix"], end)]


def crosscell_windows(log, start):
    """crosscell's 4 steps with their printed durations, laid end to end from
    process start (its long steps print whole seconds, e.g. "[6s]")."""
    name = {"1": "read", "2": "convert", "3": "convert", "4": "write"}
    out, t = [], start
    for step, dur in re.findall(r"\[(\d)/4\] \S+ .*?\[([^\]]+)\]\s*$", log, re.M):
        secs = sum(float(v) * {"ms": 1e-3, "s": 1, "m": 60, "h": 3600}[u] for v, u in re.findall(r"([\d.]+)(ms|s|m|h)\b", dur))
        out.append((name[step], t, t + secs))
        t += secs
    return out


def phases(windows, denet_file):
    """{event: [seconds, peak_rss_mb]}: summed per event, peak = denet samples
    inside the event's windows."""
    samples = []
    for line in Path(denet_file).read_text().splitlines() if Path(denet_file).exists() else []:
        try:
            a = json.loads(line).get("aggregated")
        except json.JSONDecodeError:
            continue
        if a:
            samples.append((a["ts_ms"] / 1000, a["mem_rss_kb"] / 1024))
    out = {}
    for event, t0, t1 in windows:
        peak = max((m for t, m in samples if t0 <= t <= t1), default=0.0)  # 0: no sample inside
        s, p = out.get(event, [0.0, 0.0])
        out[event] = [round(s + t1 - t0, 3), round(max(p, peak), 1)]
    return out


def measure(cmd, timeout, tag, out_dir):
    """One run under GNU time + denet. -> (ok, seconds, maxrss_mb, err, phases)."""
    raw = out_dir / ".runs"
    raw.mkdir(exist_ok=True)
    denet_file, time_file, status = raw / f"{tag}.denet.ndjson", raw / f"{tag}.time.txt", raw / f"{tag}.status"
    obkit_dir = raw / f"{tag}.obkit"
    shutil.rmtree(obkit_dir, ignore_errors=True)  # obkit appends
    obkit_dir.mkdir()
    status.unlink(missing_ok=True)
    log = raw / f"{tag}.log"  # denet -q swallows the child's stdout/stderr, so keep them here
    denet = os.environ.get("DENET") or shutil.which("denet")
    if not denet:
        sys.exit("denet not found: set DENET (static binary, side-loaded like ../omni-scx)")
    # denet exits 0 whatever its child did, so the child's status is captured
    # inside and read back (../omni-scx/modules/measure.sh).
    full = ["/usr/bin/time", "-v", "-o", str(time_file),
            denet, "-j", "-o", str(denet_file), "-q", "-i", "50", "-m", "50", "run", "--",
            "bash", "-c", 'f="$1"; l="$2"; shift 2; "$@" > "$l" 2>&1; echo $? > "$f"', "_", str(status), str(log), *cmd]
    start = time.time()
    try:
        r = subprocess.run(full, capture_output=True, text=True, timeout=timeout,
                           # one python for R tools too, as the paper's Docker image
                           env={**os.environ, "OBKIT_DIR": str(obkit_dir), "RETICULATE_PYTHON": sys.executable})
    except subprocess.TimeoutExpired:
        return False, timeout, 0.0, "TIMEOUT", {}
    elapsed = time.time() - start
    rc = int(status.read_text()) if status.exists() else 127
    m = re.search(r"Maximum resident set size \(kbytes\): (\d+)", time_file.read_text() if time_file.exists() else "")
    peak = int(m.group(1)) / 1024.0 if m else 0.0
    err = "" if rc == 0 else (log.read_text(errors="replace") if log.exists() else r.stderr).strip()
    windows = obkit_windows(obkit_dir)
    if not windows and log.exists():  # the binaries: phases from their own log
        text = log.read_text(errors="replace")
        windows = scx_windows(text, start, start + elapsed) if tag.startswith("scx.") else \
            crosscell_windows(text, start) if tag.startswith("crosscell.") else []
    return rc == 0, elapsed, peak, err, phases(windows, denet_file)


def main():
    one, _, out = parse()
    tool, runs, timeout = one["tool"], int(one["runs"]), int(one["timeout"])
    if "rds" in one:
        direction, inp, out_path = "rds_to_h5ad", one["rds"], out["out_h5ad"]
    else:
        direction, inp, out_path = "h5ad_to_rds", one["h5ad"], out["out_seurat"]
    out_path = os.path.abspath(out_path)
    declared = out_path
    if direction == "h5ad_to_rds":
        out_path = str(Path(out_path).with_suffix(NATIVE.get(tool, ".rds")))
    test_id = Path(inp).name.rsplit(".", 1)[0]
    rec = {"test_id": test_id, "tool": tool, "direction": direction}

    cmd = globals()[f"cmd_{tool}"](os.path.abspath(inp), out_path, direction)
    times, mems, ph = [], [], []
    if cmd is None:
        rec.update(status="skipped", reason=f"{tool} does not support {direction}")
    for run in range(runs if cmd else 0):
        Path(out_path).unlink(missing_ok=True)
        ok, t, mem, err, p = measure(cmd, timeout, f"{tool}.run{run + 1}", Path(out_path).parent)
        if not ok:
            rec.update(status="failed", error=err[-5000:], run=run + 1)
            break
        times.append(t); mems.append(mem); ph.append(p)
    else:
        if cmd:
            names = {k for p in ph for k in p}
            rec.update(status="success",
                       conversion_time_seconds=round(statistics.median(times), 3),
                       peak_memory_mb=round(statistics.median(mems), 1),
                       output_file_size_mb=round(os.path.getsize(out_path) / 1048576, 2),
                       all_times=[round(t, 3) for t in times], all_mems=[round(m, 1) for m in mems],
                       # median seconds / median denet peak MB, per phase
                       phases={k: [statistics.median(p.get(k, [0, 0])[0] for p in ph),
                                   statistics.median(p.get(k, [0, 0])[1] for p in ph)] for k in sorted(names)})

    if rec["status"] != "success":
        Path(out_path).unlink(missing_ok=True)
        Path(out_path).touch()  # declared output must exist; empty == no output
    if declared != out_path:
        Path(declared).unlink(missing_ok=True)
        Path(declared).symlink_to(Path(out_path).name)
    rec["output"] = Path(out_path).name
    Path(out["perf"]).write_text(json.dumps(rec, indent=1) + "\n")
    print(f"{test_id} {tool} {direction}: {rec['status']}")


if __name__ == "__main__":
    main()
