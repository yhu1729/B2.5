#!/usr/bin/env python
"""Temporal-order analysis of MMS slab runs.

Reads run_*/run.json and run_*/mms_time_error.dat (written by
b2mod_mms_analytic: one line per step), takes the final-time row of each
run, and for every scheme fits log-log slopes of error versus dt per field.
Prints a table and writes error_vs_dt.png and error_vs_sweeps.png.

  analyze_mms_time.py [--root sim/mms_slab] [--prefix run] [--fields rel_na_001,rel_te,...]
                      [--floor 1e-10] [--output DIR]
"""
import argparse
import glob
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT_ROOT = os.path.normpath(os.path.join(HERE, "..", "..", "..", "..", "..", "sim", "mms_slab"))
DEFAULT_FIELDS = ["rel_na_000", "rel_na_001", "rms_ua_000", "rms_ua_001", "rel_te", "rel_ti", "rms_po"]


def read_error_file(path):
    cols = None
    rows = []
    for line in open(path):
        if line.startswith("#"):
            if "columns:" in line:
                cols = line.split("columns:")[1].split()
            continue
        if line.strip():
            rows.append([float(x) for x in line.split()])
    if cols is None or not rows:
        raise ValueError("no data in %s" % path)
    data = np.array(rows)
    return cols, data


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", default=DEFAULT_ROOT)
    ap.add_argument("--prefix", default="run")
    ap.add_argument("--fields", default=",".join(DEFAULT_FIELDS))
    ap.add_argument("--floor", type=float, default=0.0,
                    help="ignore errors below this value in the slope fits")
    ap.add_argument("--output", default=None, help="directory for the plots (default: root)")
    ap.add_argument("--no-plot", action="store_true")
    a = ap.parse_args()

    fields = a.fields.split(",")
    runs = []
    for manifest in sorted(glob.glob(os.path.join(a.root, "%s_*" % a.prefix, "run.json"))):
        rdir = os.path.dirname(manifest)
        errf = os.path.join(rdir, "mms_time_error.dat")
        if not os.path.exists(errf):
            print("skip %s: no mms_time_error.dat" % os.path.relpath(rdir))
            continue
        m = json.load(open(manifest))
        cols, data = read_error_file(errf)
        last = data[-1]
        t = last[cols.index("t")]
        if abs(t - m["T"]) > 0.5 * m["dt"]:
            print("skip %s: final time %g differs from T = %g (run incomplete?)"
                  % (os.path.relpath(rdir), t, m["T"]))
            continue
        rec = {"scheme": m["scheme"], "dt": m["dt"], "dir": rdir,
               "sweeps": last[cols.index("n_solve_sweeps")],
               "evals": last[cols.index("n_eval_sweeps")],
               "sdc_res": last[cols.index("sdc_res")]}
        # total cost: node/step solves plus the evaluation-only sweeps (one
        # full b2mndt sweep each), which n_solve_sweeps does not include
        rec["cost"] = rec["sweeps"] + rec["evals"]
        for f in fields:
            rec[f] = last[cols.index(f)]
        runs.append(rec)
    if not runs:
        sys.exit("no complete runs found under %s" % a.root)

    schemes = sorted(set(r["scheme"] for r in runs))
    print("%-14s %-12s %6s %10s %10s %10s %s" % ("scheme", "field", "npts", "LS slope",
                                                  "min err", "max err", "pairwise slopes"))
    for sch in schemes:
        rs = sorted([r for r in runs if r["scheme"] == sch], key=lambda r: -r["dt"])
        dt = np.array([r["dt"] for r in rs])
        for f in fields:
            e = np.array([r[f] for r in rs])
            mask = e > a.floor
            if mask.sum() >= 2:
                slope = np.polyfit(np.log(dt[mask]), np.log(e[mask]), 1)[0]
            else:
                slope = float("nan")
            pw = []
            for i in range(len(rs) - 1):
                if e[i] > a.floor and e[i + 1] > a.floor and dt[i] != dt[i + 1]:
                    pw.append(np.log(e[i] / e[i + 1]) / np.log(dt[i] / dt[i + 1]))
            print("%-14s %-12s %6d %10.3f %10.3e %10.3e %s"
                  % (sch, f, mask.sum(), slope, e.min(), e.max(),
                     " ".join("%.2f" % p for p in pw)))
        print("%-14s sweeps per run (largest dt first): %s" %
              (sch, " ".join("%d" % r["sweeps"] for r in rs)))
        print("%-14s solve + eval sweeps per run:       %s" %
              (sch, " ".join("%d" % r["cost"] for r in rs)))

    if a.no_plot:
        return
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    out = a.output or a.root
    ncol = 3
    nrow = (len(fields) + ncol - 1) // ncol
    for xkey, fname, xlabel in [("dt", "error_vs_dt.png", "dt [s]"),
                                ("sweeps", "error_vs_sweeps.png", "solve sweeps"),
                                ("cost", "error_vs_cost.png", "solve + eval sweeps")]:
        fig, axes = plt.subplots(nrow, ncol, figsize=(4.5 * ncol, 3.6 * nrow), squeeze=False)
        for k, f in enumerate(fields):
            ax = axes[k // ncol][k % ncol]
            for sch in schemes:
                rs = sorted([r for r in runs if r["scheme"] == sch], key=lambda r: r["dt"])
                x = [r[xkey] for r in rs]
                y = [r[f] for r in rs]
                ax.loglog(x, y, "o-", label=sch)
            if xkey == "dt":
                dts = np.array(sorted(set(r["dt"] for r in runs)))
                ymax = max(r[f] for r in runs)
                for p, ls in [(1, ":"), (2, "--"), (3, "-."), (5, ":")]:
                    ax.loglog(dts, ymax * (dts / dts.max()) ** p, "k" + ls, lw=0.8,
                              label="slope %d" % p if k == 0 else None)
            ax.set_title(f)
            ax.set_xlabel(xlabel)
            ax.grid(True, which="both", alpha=0.3)
        axes[0][0].legend(fontsize=7)
        for k in range(len(fields), nrow * ncol):
            axes[k // ncol][k % ncol].axis("off")
        fig.tight_layout()
        path = os.path.join(out, fname)
        fig.savefig(path, dpi=120)
        print("wrote %s" % path)


if __name__ == "__main__":
    main()
