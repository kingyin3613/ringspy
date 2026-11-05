# -*- coding: utf-8 -*-
"""Shared helpers of the tests"""

import re
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


def without_ring_labels(name, data):
    """File content (bytes) without the growth ring labels that generate() adds since v1.4: the Ring/Dense
    cell data of the .vtu files, the ring element sets of the .inp file, and the transition and ring radii
    lines of the report log. Returns (content without labels, whether labels were found)."""
    if name.endswith(".vtu"):
        text = data.decode("ascii")
        # the files are written in text mode: \r\n line endings on Windows
        stripped = re.sub(r'<DataArray type="Int32" Name="(Ring|Dense)" format="ascii">\r?\n(?:-?\d+\r?\n)*</DataArray>\r?\n', "", text)
        stripped = re.sub(r'<CellData Tensors="" Vectors="" Scalars="Ring">\r?\n</CellData>\r?\n', "", stripped)
        return stripped.encode("ascii"), stripped != text
    if name.endswith(".inp"):
        text = data.decode("ascii")
        stripped = re.sub(r"\*\* Element sets of the growth rings.*?(?=\*\* Section: Section-1\r?\n)", "", text, flags=re.S)
        return stripped.encode("ascii"), stripped != text
    if name.endswith("-report.log"):
        lines = data.decode().splitlines()
        kept = [line for line in lines if not line.startswith(("transition:", "transition_width:", "ring radii:"))]
        return "\n".join(kept).encode(), len(kept) < len(lines)
    return data, False


def assert_same_files(dir_a, dir_b, ignore_ring_labels=False):
    """The two output folders hold the same files with the same content (report logs: except run times).
    ignore_ring_labels: the files of dir_b may carry the growth ring labels (and must, where they have them)"""
    names_a = sorted(p.name for p in Path(dir_a).iterdir())
    names_b = sorted(p.name for p in Path(dir_b).iterdir())
    assert names_a == names_b
    for name in names_a:
        a, b = Path(dir_a, name), Path(dir_b, name)
        data_a, data_b = a.read_bytes(), b.read_bytes()
        if ignore_ring_labels:
            _, labels_a = without_ring_labels(name, data_a)
            data_b, labels_b = without_ring_labels(name, data_b)
            assert not labels_a and labels_b == name.endswith((".vtu", ".inp", "-report.log")), name
        if name.endswith("-report.log"):
            assert without_timing(data_a.decode()) == without_timing(data_b.decode()), name
        elif name.endswith(".stl") and data_a != data_b:
            # the same triangles: an STL made from the ASCII .vtu file (9 significant digits) differs in the last
            # bits of its float32 vertices from one made from the exact coordinates
            import numpy as np
            record = np.dtype([("normal", "<f4", 3), ("vertices", "<f4", (3, 3)), ("attribute", "<u2")])
            ta, tb = np.frombuffer(data_a[84:], dtype=record), np.frombuffer(data_b[84:], dtype=record)
            assert len(ta) == len(tb), name
            np.testing.assert_allclose(ta["vertices"], tb["vertices"], rtol=1e-6, atol=1e-6, err_msg=name)
        else:
            assert data_a == data_b, name
