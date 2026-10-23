# -*- coding: utf-8 -*-
"""Shared helpers of the tests"""

from pathlib import Path

from matplotlib.path import Path as MplPath

ROOT = Path(__file__).resolve().parents[1]
EXAMPLES_DIR = ROOT / "examples"
EXAMPLES = sorted(EXAMPLES_DIR.glob("*.toml"))

# files written for every geometry (plus .stl and .inp if stlFlag and inpFlag are on)
OUTPUT_SUFFIXES = ["-vertex.mesh", "-ridge.mesh", "-mesh.txt", "IGA.txt", "_vertices.vtu", "_beams.vtu",
                   "_conns.vtu", "_conns_vol.vtu", ".png", "-report.log"]


def inside_polygon(points, polygon, tol=1e-9):
    path = MplPath(polygon)
    return path.contains_points(points, radius=tol) | path.contains_points(points, radius=-tol)


def without_timing(log_text):
    """Report log without the lines that change from run to run (date and run times)"""
    return [line for line in log_text.splitlines() if "local time" not in line and "Time:" not in line]


def assert_same_files(dir_a, dir_b):
    """The two output folders hold the same files with the same content (report logs: except run times)"""
    names_a = sorted(p.name for p in Path(dir_a).iterdir())
    names_b = sorted(p.name for p in Path(dir_b).iterdir())
    assert names_a == names_b
    for name in names_a:
        a, b = Path(dir_a, name), Path(dir_b, name)
        if name.endswith("-report.log"):
            assert without_timing(a.read_text()) == without_timing(b.read_text()), name
        else:
            assert a.read_bytes() == b.read_bytes(), name
