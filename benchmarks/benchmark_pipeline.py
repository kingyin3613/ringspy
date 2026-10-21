# -*- coding: utf-8 -*-
"""
Wall time per stage of a RingsPy template script (default: tests/test_wood_cube.py).

    python benchmarks/benchmark_pipeline.py --iter-max 500 2000 --seed 0

Runs the template once per iter_max value, with np.random seeded, writing to a
temporary folder (or --out), and prints the time spent in each MeshGenTools
function the template calls, plus the matplotlib preview and the total wall time.
Not a test: it is meant for comparing versions of RingsPy on the same machine.
"""

import argparse
import json
import re
import sys
import tempfile
import time
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import scipy.spatial

import ringspy.MeshGenTools as rpgen

TESTS_DIR = Path(__file__).resolve().parent.parent / "tests"

# MeshGenTools functions called by the template scripts, in pipeline order
STAGES = ["CellPlacement_Binary", "CellPlacement_Binary_Lloyd", "CellPlacement_Honeycomb", "ReadSavedSites",
          "Clipping_Box", "RebuildVoronoi", "RebuildVoronoi_merge", "RidgeMidQuarterPts", "VertexandRidgeinfo",
          "GenerateBeamElement", "insert_precracks", "ConnectorMeshFile", "ModelInfo", "BezierExtraction",
          "BezierBeamFile", "VisualizationFiles", "StlModelFile", "AbaqusFile", "LogFile"]


def timed(name, func, times):
    def wrapper(*args, **kwargs):
        start = time.perf_counter()
        try:
            return func(*args, **kwargs)
        finally:
            times[name] = times.get(name, 0.0) + time.perf_counter() - start
    return wrapper


def run(template, iter_max, seed, out):
    """Run the template once, return {stage: seconds}, total seconds and the number of cells"""
    src = (TESTS_DIR / (template + ".py")).read_text()
    overrides = {"path": str(out), "iter_max": iter_max, "print_interval": 10**9}
    for key, value in overrides.items():
        src, n = re.subn(r"^(\s*)" + key + r" = .*$", lambda m: m.group(1) + key + " = " + repr(value),
                         src, count=1, flags=re.M)
        if n != 1:
            raise ValueError("{} not found in {}".format(key, template))

    times = {}
    originals = {}
    patches = [(rpgen, name) for name in STAGES] + \
              [(scipy.spatial, "Voronoi"), (scipy.spatial, "voronoi_plot_2d"), (plt, "savefig")]
    for module, name in patches:
        originals[(module, name)] = getattr(module, name)
        setattr(module, name, timed(name, getattr(module, name), times))
    try:
        plt.close("all")
        np.random.seed(seed)
        namespace = {"__file__": str(TESTS_DIR / (template + ".py")), "__name__": "benchmark"}
        exec(compile(src, template + ".py", "exec"), namespace)
        start = time.perf_counter()
        namespace[template]()
        total = time.perf_counter() - start
    finally:
        for (module, name), func in originals.items():
            setattr(module, name, func)
        plt.close("all")

    geo = next(p for p in Path(out).iterdir() if p.is_dir())
    log = next(geo.glob("*-report.log")).read_text()
    ncells = int(re.search(r"nParticles:\s+(\d+)", log).group(1))
    size_mb = sum(f.stat().st_size for f in geo.iterdir()) / 1024**2
    return times, total, ncells, size_mb


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--template", default="test_wood_cube", help="template script in tests/ (without .py)")
    parser.add_argument("--iter-max", type=int, nargs="+", default=[500])
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--out", help="output folder (default: a temporary folder, deleted afterwards)")
    parser.add_argument("--json", help="append one JSON line per run to this file")
    args = parser.parse_args()

    print("ringspy {} from {}".format(rpgen.__version__, Path(rpgen.__file__).parent))
    for iter_max in args.iter_max:
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(args.out, "iter_max_{}".format(iter_max)).resolve() if args.out else Path(tmp)
            out.mkdir(parents=True, exist_ok=True)
            times, total, ncells, size_mb = run(args.template, iter_max, args.seed, out)
        print("\n{} iter_max={} seed={}: {} cells, {:.0f} MB of output".format(
            args.template, iter_max, args.seed, ncells, size_mb))
        for name in sorted(times, key=times.get, reverse=True):
            print("  {:28s} {:9.3f} s".format(name, times[name]))
        print("  {:28s} {:9.3f} s".format("other", total - sum(times.values())))
        print("  {:28s} {:9.3f} s".format("TOTAL", total))
        rss_mb = peak_rss_mb()
        if rss_mb is not None:
            print("  peak memory so far: {:.0f} MB".format(rss_mb))
        sys.stdout.flush()
        if args.json:
            with open(args.json, "a") as f:
                f.write(json.dumps({"template": args.template, "iter_max": iter_max, "seed": args.seed,
                                    "ncells": ncells, "size_mb": size_mb, "total": total, "stages": times,
                                    "peak_rss_mb": rss_mb, "version": rpgen.__version__}) + "\n")


def peak_rss_mb():
    """Peak resident memory of this process (Linux and macOS only)"""
    try:
        import resource
    except ImportError:
        return None
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return rss / 1024**2 if sys.platform == "darwin" else rss / 1024


if __name__ == "__main__":
    main()
