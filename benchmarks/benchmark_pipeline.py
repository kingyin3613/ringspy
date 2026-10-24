# -*- coding: utf-8 -*-
"""
Wall time per stage of ringspy.generate() for an example configuration
(default: examples/wood_cube.toml).

    python benchmarks/benchmark_pipeline.py --iter-max 500 2000 --seed 0

Generates the mesh once per iter_max value, writing to a temporary folder (or --out),
and prints the time spent in each MeshGenTools function, the Voronoi tessellation and
the preview plot, plus the total wall time. Not a test: it is meant for comparing
versions of RingsPy on the same machine.
"""

import argparse
import json
import sys
import tempfile
import time
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.figure

import ringspy
import ringspy.MeshGenTools as rpgen
import ringspy.pipeline

EXAMPLES_DIR = Path(__file__).resolve().parent.parent / "examples"

# MeshGenTools functions called by ringspy.generate(), in pipeline order
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


def run(config, iter_max, seed, out):
    """Generate the mesh once, return {stage: seconds}, total seconds, the number of cells and the output size"""
    times = {}
    originals = {}
    patches = [(rpgen, name) for name in STAGES] + \
              [(ringspy.pipeline, "Voronoi"), (ringspy.pipeline, "voronoi_plot_2d"), (matplotlib.figure.Figure, "savefig")]
    for owner, name in patches:
        originals[(owner, name)] = getattr(owner, name)
        setattr(owner, name, timed(name, getattr(owner, name), times))
    try:
        start = time.perf_counter()
        result = ringspy.generate(config, iter_max=iter_max, seed=seed, path=str(out), print_interval=10**9)
        total = time.perf_counter() - start
    finally:
        for (owner, name), func in originals.items():
            setattr(owner, name, func)
    size_mb = sum(f.stat().st_size for f in result.files) / 1024**2
    return times, total, result.nParticles, size_mb


def peak_rss_mb():
    """Peak resident memory of this process (Linux and macOS only)"""
    try:
        import resource
    except ImportError:
        return None
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return rss / 1024**2 if sys.platform == "darwin" else rss / 1024


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--config", default=str(EXAMPLES_DIR / "wood_cube.toml"), help="TOML configuration file")
    parser.add_argument("--iter-max", type=int, nargs="+", default=[500])
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--out", help="output folder (default: a temporary folder, deleted afterwards)")
    parser.add_argument("--json", help="append one JSON line per run to this file")
    args = parser.parse_args()

    print("ringspy {} from {}".format(ringspy.__version__, Path(rpgen.__file__).parent))
    for iter_max in args.iter_max:
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(args.out, "iter_max_{}".format(iter_max)).resolve() if args.out else Path(tmp)
            times, total, ncells, size_mb = run(args.config, iter_max, args.seed, out)
        print("\n{} iter_max={} seed={}: {} cells, {:.0f} MB of output".format(
            Path(args.config).stem, iter_max, args.seed, ncells, size_mb))
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
                f.write(json.dumps({"config": args.config, "iter_max": iter_max, "seed": args.seed,
                                    "ncells": ncells, "size_mb": size_mb, "total": total, "stages": times,
                                    "peak_rss_mb": rss_mb, "version": ringspy.__version__}) + "\n")


if __name__ == "__main__":
    main()
