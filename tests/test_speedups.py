# -*- coding: utf-8 -*-
"""
Checks that the faster implementations give the same results as the RingsPy 1.1.0
implementations they replaced.
"""

import functools
import math
from pathlib import Path

import numpy as np
import pytest

import ringspy
import ringspy.MeshGenTools as rpgen
from helpers import EXAMPLES_DIR

pytestmark = pytest.mark.usefixtures("keep_global_random_state") # these tests seed np.random


def brute_force_binary_placement(generation_center, r_max, r_min, nrings, width_heart, width_sparse,
                                 width_dense, cellsize_sparse, cellsize_dense, iter_max, rng=np.random):
    """CellPlacement_Binary of RingsPy 1.1.0: every try is checked against all placed cells"""
    width = np.concatenate(([width_heart], np.tile([width_sparse, width_dense], nrings)))
    noise = rng.normal(1, 0.25, len(width))
    width = np.multiply(width, noise)
    radii = np.concatenate(([0], np.cumsum(width)))
    circles = []
    for ibin in range(0, nrings*2+1):
        for _ in range(iter_max):
            r = (radii[ibin+1]-radii[ibin])*math.sqrt(rng.random()) + radii[ibin]
            t = rng.random()*2*math.pi
            x = r*math.cos(t) + generation_center[0]
            y = r*math.sin(t) + generation_center[1]
            w = cellsize_dense if ibin % 2 == 0 else cellsize_sparse
            circC = [x, y, w]
            if rpgen.check_overlap(circles, circC) and \
                    rpgen.check_isinside_boundcircle(circC, generation_center, radii[ibin], radii[ibin+1]):
                circles.append(circC)
    sites = np.array(circles)
    outofbound = [i for i in range(sites.shape[0])
                  if (sites[i,0]-generation_center[0])**2 + (sites[i,1]-generation_center[1])**2 > (1.2*r_max)**2]
    return np.delete(sites, outofbound, 0), radii


@pytest.mark.parametrize("random_source", ["global", "generator"])
@pytest.mark.parametrize("seed", [0, 1, 2])
@pytest.mark.parametrize("generation_center,cellsize_sparse,cellsize_dense", [
    ((0.0, 0.0), 0.04, 0.02),
    ((0.3, -5.2), 0.03, 0.03),
])
def test_binary_placement_same_sites_as_brute_force(random_source, seed, generation_center, cellsize_sparse,
                                                     cellsize_dense):
    # the grid search only skips cells too far away to overlap: same accepted sites, same random numbers
    args = (generation_center, 1.0, 0.0, 2, 0.15, 0.35, 0.15, cellsize_sparse, cellsize_dense, 300)
    if random_source == "global":
        np.random.seed(seed)
        expected_sites, expected_radii = brute_force_binary_placement(*args)
        expected_next_random = np.random.random()
        np.random.seed(seed)
        sites, radii = rpgen.CellPlacement_Binary(*args, 10**9)
        next_random = np.random.random()
    else:
        rng = np.random.default_rng(seed)
        expected_sites, expected_radii = brute_force_binary_placement(*args, rng=rng)
        expected_next_random = rng.random()
        rng = np.random.default_rng(seed)
        sites, radii = rpgen.CellPlacement_Binary(*args, 10**9, rng=rng)
        next_random = rng.random()
    assert len(sites) > 100
    np.testing.assert_array_equal(sites, expected_sites)
    np.testing.assert_array_equal(radii, expected_radii)
    assert next_random == expected_next_random


@pytest.mark.parametrize("fmt", ["%.16g", "%d", "%11.8e"])
@pytest.mark.parametrize("comments", ["", "# "])
def test_savetxt_writes_same_text_as_numpy(tmp_path, fmt, comments):
    rng = np.random.default_rng(0)
    if fmt == "%d":
        X = rng.integers(-10**6, 10**6, size=(45001, 6))
    else:
        X = rng.normal(size=(45001, 6)) * 10.0**rng.integers(-8, 8, size=(45001, 6))
        X[0, 0:3] = [0.0, -0.0, 1.0]
    header = "Data\nNumber of rows\n" + str(len(X))
    np.savetxt(tmp_path / "numpy.txt", X, fmt=fmt, delimiter=" ", header=header, comments=comments)
    rpgen._savetxt(tmp_path / "ringspy.txt", X, fmt=fmt, delimiter=" ", header=header, comments=comments)
    assert (tmp_path / "ringspy.txt").read_bytes() == (tmp_path / "numpy.txt").read_bytes()


def read_vtu(path):
    """Points, cells and data arrays of a .vtu file, read with VTK"""
    from vtkmodules.vtkIOXML import vtkXMLUnstructuredGridReader
    from vtkmodules.util.numpy_support import vtk_to_numpy
    reader = vtkXMLUnstructuredGridReader()
    reader.SetFileName(str(path))
    reader.Update()
    grid = reader.GetOutput()
    arrays = {"points": vtk_to_numpy(grid.GetPoints().GetData()),
              "connectivity": vtk_to_numpy(grid.GetCells().GetConnectivityArray()),
              "offsets": vtk_to_numpy(grid.GetCells().GetOffsetsArray()),
              "types": vtk_to_numpy(grid.GetCellTypes() if hasattr(grid, "GetCellTypes") else grid.GetCellTypesArray())}
    for kind, data in [("point", grid.GetPointData()), ("cell", grid.GetCellData())]:
        for i in range(data.GetNumberOfArrays()):
            arrays[kind + ":" + data.GetArrayName(i)] = vtk_to_numpy(data.GetArray(i))
    return arrays


def test_binary_vtu_files_match_ascii(monkeypatch):
    pytest.importorskip("vtkmodules.vtkIOXML")
    honeycomb_hex = EXAMPLES_DIR / "honeycomb_hex.toml"
    ringspy.generate(honeycomb_hex, geoName="ascii", seed=0)
    monkeypatch.setattr(rpgen, "VisualizationFiles", functools.partial(rpgen.VisualizationFiles, vtk_format="binary"))
    ringspy.generate(honeycomb_hex, geoName="binary", seed=0)
    for suffix in ["_vertices", "_beams", "_conns", "_conns_vol"]:
        ascii_file = Path("meshes", "ascii", "ascii" + suffix + ".vtu")
        binary_file = Path("meshes", "binary", "binary" + suffix + ".vtu")
        assert binary_file.stat().st_size < ascii_file.stat().st_size
        expected, arrays = read_vtu(ascii_file), read_vtu(binary_file)
        assert arrays.keys() == expected.keys()
        for name in expected:
            if name in ["connectivity", "offsets", "types"]:
                np.testing.assert_array_equal(arrays[name], expected[name], err_msg=suffix + " " + name)
            else: # the ASCII files have 9 significant digits
                np.testing.assert_allclose(arrays[name], expected[name], rtol=1e-8, atol=1e-12, err_msg=suffix + " " + name)
    # the STL file is made from the binary connector volume file: same surface quads, each split
    # into two triangles (not always along the same diagonal, the coordinates differ slightly)
    quads = []
    for name in ["ascii", "binary"]:
        stl = Path("meshes", name, name + ".stl").read_bytes()
        triangles = np.frombuffer(stl[84:], dtype=np.dtype([("normal", "<f4", 3), ("vertices", "<f4", (3, 3)),
                                                            ("attribute", "<u2")]))["vertices"].astype(float)
        area = 0.5*np.linalg.norm(np.cross(triangles[:,1] - triangles[:,0], triangles[:,2] - triangles[:,0]), axis=1)
        centroid = triangles.mean(axis=1)
        quad_area = area[0::2] + area[1::2]
        quad_centroid = (area[0::2,None]*centroid[0::2] + area[1::2,None]*centroid[1::2])/quad_area[:,None]
        quads.append((quad_area, quad_centroid))
    assert len(quads[1][0]) == len(quads[0][0]) > 0
    np.testing.assert_allclose(quads[1][0], quads[0][0], rtol=1e-5)
    np.testing.assert_allclose(quads[1][1], quads[0][1], atol=1e-6)


def test_unknown_vtk_format_raises():
    with pytest.raises(ValueError):
        rpgen.VisualizationFiles(*([None]*22), vtk_format="xml")
