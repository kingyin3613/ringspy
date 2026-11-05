# -*- coding: utf-8 -*-
"""
The STL file without VTK: the outer surface of the connector volumes (MeshGenTools.HexahedraSurface,
WriteStl, ReadVtuCells, StlModelFile)
"""

import functools
from pathlib import Path

import numpy as np
import pytest

import ringspy
import ringspy.MeshGenTools as rpgen
from helpers import EXAMPLES_DIR

STL_RECORD = np.dtype([("normal", "<f4", 3), ("vertices", "<f4", (3, 3)), ("attribute", "<u2")])


def read_stl(file):
    data = Path(file).read_bytes()
    n = int.from_bytes(data[80:84], "little")
    records = np.frombuffer(data[84:], dtype=STL_RECORD)
    assert len(records) == n and len(data) == 84 + 50*n
    return records["normal"].astype(float), records["vertices"].astype(float)


def triangle_area(t):
    return 0.5*np.linalg.norm(np.cross(t[:, 1] - t[:, 0], t[:, 2] - t[:, 0]), axis=1)


def cube(origin, size=1.0):
    o = np.asarray(origin, dtype=float)
    return o + size*np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0], [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1]], dtype=float)


def test_hexahedra_surface_of_cubes():
    points = cube((0, 0, 0))
    triangles = rpgen.HexahedraSurface(points, [np.arange(8)])
    assert triangles.shape == (12, 3, 3) and triangle_area(triangles).sum() == pytest.approx(6.0)
    # outward normals: every triangle's normal points away from the centre
    normals = np.cross(triangles[:, 1] - triangles[:, 0], triangles[:, 2] - triangles[:, 0])
    assert np.all(np.einsum("ij,ij->i", normals, triangles.mean(axis=1) - 0.5) > 0)
    # the same cube with its vertices in the inside-out order is oriented outwards too
    triangles = rpgen.HexahedraSurface(points, [[0, 3, 2, 1, 4, 7, 6, 5]])
    normals = np.cross(triangles[:, 1] - triangles[:, 0], triangles[:, 2] - triangles[:, 0])
    assert np.all(np.einsum("ij,ij->i", normals, triangles.mean(axis=1) - 0.5) > 0)
    # two cubes sharing a face (the same point indices): the shared face is left out
    points = np.vstack((cube((0, 0, 0)), cube((1, 0, 0))[[1, 2, 5, 6]]))  # the second cube reuses 4 points of the first
    second = [1, 8, 9, 2, 5, 10, 11, 6]
    triangles = rpgen.HexahedraSurface(points, [np.arange(8), second])
    assert triangles.shape == (20, 3, 3) and triangle_area(triangles).sum() == pytest.approx(10.0)
    # two cubes touching with different points: both faces kept (as VTK does)
    points = np.vstack((cube((0, 0, 0)), cube((1, 0, 0))))
    triangles = rpgen.HexahedraSurface(points, [np.arange(8), np.arange(8, 16)])
    assert triangles.shape == (24, 3, 3)
    assert rpgen.HexahedraSurface(points, np.zeros((0, 8), dtype=int)).shape == (0, 3, 3)


def test_write_and_read_stl_and_vtu(tmp_path):
    points = cube((0, 0, 0))
    triangles = rpgen.HexahedraSurface(points, [np.arange(8)])
    rpgen.WriteStl(tmp_path / "cube.stl", triangles)
    normals, vertices = read_stl(tmp_path / "cube.stl")
    np.testing.assert_allclose(vertices, triangles)
    np.testing.assert_allclose(np.linalg.norm(normals, axis=1), 1.0)
    assert (tmp_path / "cube.stl").read_bytes()[:18] == b"RingsPy binary STL"
    # the .vtu files written by RingsPy read back, ASCII and binary
    cells = [np.arange(8).reshape(1, 8)]
    files = [("_vol", [points], cells, np.array([12]), [], [("Connector_width", np.array([0.1]))])]
    (tmp_path / "cube").mkdir()
    rpgen._write_vtu_files_ascii(tmp_path, "cube", files)  # writes tmp_path/cube/cube_vol.vtu
    rpgen._write_vtu_appended(tmp_path / "cube" / "cube_bin.vtu", points, cells, np.array([12]), None,
                              [("Connector_width", np.array([0.1]))])
    for name in ["cube/cube_vol.vtu", "cube/cube_bin.vtu"]:
        p, c = rpgen.ReadVtuCells(tmp_path / name)
        np.testing.assert_allclose(p, points, rtol=1e-8)
        np.testing.assert_array_equal(c, np.arange(8).reshape(1, 8))


def test_stl_of_a_mesh_without_vtk(tmp_path, monkeypatch):
    res = ringspy.generate(dict(geoName="stl", r_max=1.5, nrings=2, box_size=1.0, box_center=[0.75, 0], iter_max=200,
                                seed=1, inpFlag=False, stlFlag=True, path=str(tmp_path)))
    normals, triangles = read_stl(res.output_dir / "stl.stl")
    assert len(triangles) > 0
    # the signed volume enclosed by the surface is the volume of the connectors (the mass over the density)
    volume = np.sum(np.einsum("ij,ij->i", normals*triangle_area(triangles)[:, None], triangles.mean(axis=1)))/3
    assert volume == pytest.approx(res.mass/res.config.skeleton_density, rel=1e-3)
    # the same surface from the connector volume file alone, as old scripts call it: the ASCII file has 9
    # significant digits, the binary one is exact
    direct = (res.output_dir / "stl.stl").read_bytes()
    rpgen.StlModelFile("stl", path=str(tmp_path))
    _, again = read_stl(res.output_dir / "stl.stl")
    np.testing.assert_allclose(again, triangles, rtol=1e-6, atol=1e-7)
    monkeypatch.setattr(rpgen, "VisualizationFiles", functools.partial(rpgen.VisualizationFiles, vtk_format="binary"))
    res2 = ringspy.generate(dict(geoName="stlb", r_max=1.5, nrings=2, box_size=1.0, box_center=[0.75, 0], iter_max=200,
                                 seed=1, inpFlag=False, stlFlag=False, path=str(tmp_path)))
    rpgen.StlModelFile("stlb", path=str(tmp_path))
    assert (res2.output_dir / "stlb.stl").read_bytes() == direct


def test_stl_matches_vtk_surface(tmp_path):
    # with VTK installed: the same triangles (count and total area) as VTK's surface filter gives
    vtkIOXML = pytest.importorskip("vtkmodules.vtkIOXML")
    from vtkmodules.vtkFiltersGeometry import vtkDataSetSurfaceFilter
    from vtkmodules.vtkFiltersCore import vtkTriangleFilter
    from vtkmodules.vtkIOGeometry import vtkSTLWriter
    res = ringspy.generate(dict(geoName="cmp", r_max=1.5, nrings=2, box_size=1.0, box_center=[0.75, 0], iter_max=200,
                                seed=1, inpFlag=False, stlFlag=True, path=str(tmp_path)))
    reader = vtkIOXML.vtkXMLUnstructuredGridReader()
    reader.SetFileName(str(res.output_dir / "cmp_conns_vol.vtu"))
    surface = vtkDataSetSurfaceFilter()
    surface.SetInputConnection(reader.GetOutputPort())
    triangles = vtkTriangleFilter()
    triangles.SetInputConnection(surface.GetOutputPort())
    writer = vtkSTLWriter()
    writer.SetFileTypeToBinary()
    writer.SetFileName(str(tmp_path / "vtk.stl"))
    writer.SetInputConnection(triangles.GetOutputPort())
    writer.Write()
    _, ours = read_stl(res.output_dir / "cmp.stl")
    _, theirs = read_stl(tmp_path / "vtk.stl")
    assert len(ours) == len(theirs)
    assert triangle_area(ours).sum() == pytest.approx(triangle_area(theirs).sum(), rel=1e-5)


def test_fibers_and_examples_write_stl_without_vtk(tmp_path):
    res = ringspy.generate(dict(geoName="fib", r_max=1.0, nrings=2, box_size=0.5, box_center=[0.55, 0], z_max=1.0, seed=3,
                                iter_max=200, nsegments=4, fiber_length=0.4, inpFlag=False, stlFlag=True, path=str(tmp_path)))
    normals, triangles = read_stl(res.output_dir / "fib.stl")
    assert len(triangles) > 0 and np.all(np.isfinite(triangles))
