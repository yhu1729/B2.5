#!/usr/bin/env python
"""Set up (and optionally run) a dt-convergence study on the MMS slab case.

Creates sim/mms_slab/run_<scheme>_dt<value>/ from sim/mms_slab/templates/,
one directory per (scheme, dt), each with its own b2mn.dat, the four
b2.*.parameters files, clean.csh and a run.json manifest, then calls
sim/mms_slab/run_mms.csh on each directory in sequence.

Schemes: ie (implicit Euler), cn (Crank-Nicolson), bdf2, sdc:M:K[:tol]
(Radau IIA with M nodes, K sweeps, optional residual tolerance) and
sdclu:M:K[:tol] (same, with the LU-trick Q_delta sweep, b2mndt_sdc_qdelta = 1).

Example:
  mms_study.py --dts 1e-4,5e-5,2.5e-5,1.25e-5 --schemes ie,bdf2,sdc:2:3 --T 1e-3
"""
import argparse
import json
import os
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT_ROOT = os.path.normpath(os.path.join(HERE, "..", "..", "..", "..", "..", "sim", "mms_slab"))
PARAM_FILES = ["b2.boundary.parameters", "b2.numerics.parameters",
               "b2.transport.parameters", "b2.neutrals.parameters", "clean.csh"]


def dt_tag(dt):
    return ("%.3e" % dt).replace("-", "m").replace("+", "p").replace(".", "_")


def parse_scheme(s):
    parts = s.split(":")
    name = parts[0].lower()
    sw = {"CKN": "1", "BDF": "1", "SDC": "0", "SDC_NODES": "2", "SDC_SWEEPS": "0", "SDC_TOL": "0.0",
          "SDC_QDELTA": "0"}
    if name == "ie":
        pass
    elif name == "cn":
        sw["CKN"] = "2"
    elif name == "bdf2":
        sw["BDF"] = "2"
    elif name in ("sdc", "sdclu"):
        if len(parts) < 3:
            sys.exit("%s scheme needs %s:M:K[:tol]" % (name, name))
        sw["SDC"] = "1"
        sw["SDC_NODES"] = str(int(parts[1]))
        sw["SDC_SWEEPS"] = str(int(parts[2]))
        sw["SDC_QDELTA"] = "1" if name == "sdclu" else "0"
        if len(parts) > 3:
            sw["SDC_TOL"] = "%g" % float(parts[3])
        name = "%s%s_%s" % (name, parts[1], parts[2])
        if len(parts) > 3:
            name += "_tol%s" % parts[3].replace("-", "m").replace("+", "p").replace(".", "_")
    else:
        sys.exit("unknown scheme %s" % s)
    return name, sw


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dts", required=True, help="comma-separated time steps [s]")
    ap.add_argument("--schemes", required=True, help="comma-separated: ie,cn,bdf2,sdc:M:K[:tol]")
    ap.add_argument("--T", type=float, default=1.0e-3, help="final time [s]")
    ap.add_argument("--period", type=float, default=1.0e-3, help="MMS temporal period [s]")
    ap.add_argument("--nstg2", type=int, default=1, help="inner sweeps per (node) step")
    ap.add_argument("--nstg0", type=int, default=2,
                    help="outer iterations per (node) step: >1 re-linearises sources and rates at the current iterate")
    ap.add_argument("--rxf", type=float, default=1.0, help="b2mndt_rxf")
    ap.add_argument("--kinsol", type=int, default=1, help="b2mndt_use_kinsol")
    ap.add_argument("--kinsol-tol", type=float, default=1.0e-9, help="b2mndt_kinsol_fnorm_tol")
    ap.add_argument("--kinsol-sral", type=int, default=1,
                    help="b2mndt_kinsol_sral: re-linearise sources inside the KINSOL map (1 = nonlinear fixed point)")
    ap.add_argument("--pbig", type=float, default=1.0e6,
                    help="b2stbc_pbig_factor: multiplier of the BCPOT=7 penalty (code default 1 is too weak for MMS)")
    ap.add_argument("--root", default=DEFAULT_ROOT, help="sim/mms_slab directory")
    ap.add_argument("--prefix", default="run", help="run directory prefix")
    ap.add_argument("--no-run", action="store_true", help="only create the directories")
    ap.add_argument("--force", action="store_true", help="overwrite an existing b2mn.dat")
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()

    root = os.path.abspath(a.root)
    tdir = os.path.join(root, "templates")
    template = open(os.path.join(tdir, "b2mn.dat.template")).read()
    dts = [float(x) for x in a.dts.split(",")]
    schemes = [parse_scheme(s) for s in a.schemes.split(",")]
    launcher = os.path.join(root, "run_mms.csh")

    made = []
    for name, sw in schemes:
        for dt in dts:
            ntim = int(round(a.T / dt))
            if abs(ntim * dt - a.T) > 1e-9 * a.T:
                sys.exit("T = %g is not a multiple of dt = %g" % (a.T, dt))
            rdir = os.path.join(root, "%s_%s_dt%s" % (a.prefix, name, dt_tag(dt)))
            subst = dict(sw)
            subst.update({"SCHEME_NAME": name, "NTIM": str(ntim), "DTIM": "%.12e" % dt,
                          "NSTG2": str(a.nstg2), "NSTG0": str(a.nstg0), "RXF": "%g" % a.rxf,
                          "PERIOD": "%.12e" % a.period,
                          "USE_KINSOL": str(a.kinsol), "KINSOL_TOL": "%g" % a.kinsol_tol,
                          "PBIG": "%g" % a.pbig, "KINSOL_SRAL": str(a.kinsol_sral)})
            text = template
            for k, v in subst.items():
                text = text.replace("@%s@" % k, v)
            if "@" in text:
                sys.exit("unreplaced placeholder in template: %s" %
                         [l for l in text.splitlines() if "@" in l][0])
            manifest = {"scheme": name, "dt": dt, "ntim": ntim, "T": a.T, "period": a.period,
                        "nstg2": a.nstg2, "nstg0": a.nstg0, "rxf": a.rxf, "kinsol": a.kinsol,
                        "kinsol_tol": a.kinsol_tol, "pbig": a.pbig, "kinsol_sral": a.kinsol_sral, "switches": sw}
            print("%s: ntim = %d" % (os.path.relpath(rdir), ntim))
            if a.dry_run:
                continue
            os.makedirs(rdir, exist_ok=True)
            b2mn = os.path.join(rdir, "b2mn.dat")
            if os.path.exists(b2mn) and not a.force:
                print("  b2mn.dat exists, kept (use --force to overwrite)")
            else:
                open(b2mn, "w").write(text)
            for f in PARAM_FILES:
                dst = os.path.join(rdir, f)
                if not os.path.exists(dst) or a.force:
                    shutil.copy(os.path.join(tdir, f), dst)
            json.dump(manifest, open(os.path.join(rdir, "run.json"), "w"), indent=1)
            made.append(rdir)

    if a.dry_run or a.no_run:
        return
    for rdir in made:
        print("=== running %s" % os.path.relpath(rdir))
        rc = subprocess.call(["tcsh", launcher, rdir])
        if rc != 0:
            print("run_mms.csh failed with status %d in %s" % (rc, rdir))
            sys.exit(rc)


if __name__ == "__main__":
    main()
