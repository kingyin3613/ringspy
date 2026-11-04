# -*- coding: utf-8 -*-
"""
Lattices with fibers of various lengths (ringspy.lattice), in the mesh generation (fiber_length)
"""

import numpy as np
import pytest

import ringspy
import ringspy.MeshGenTools as rpgen
from ringspy import lattice as L
from helpers import EXAMPLES_DIR, OUTPUT_SUFFIXES

SMALL = ringspy.Config(geoName="small", r_max=1.0, nrings=2, box_size=0.5, box_center=(0.55, 0.0), z_max=1.0, seed=3,
                       iter_max=300, nsegments=4, stlFlag=False, inpFlag=False)


def mesh2d(config):
    """The 2D mesh arrays of a Config, as _generate makes them"""
    import os
    from scipy.spatial import Voronoi
    from ringspy.pipeline import place_cells
    os.makedirs("mesh2d", exist_ok=True)
    c = config.resolved()
    sites, radii = place_cells(c, np.random.default_rng(c.seed))
    x_min, x_max, y_min, y_max, boundaries, boundary_points, _ = rpgen.Clipping_Box(c.box_shape, c.box_center, c.box_size, True)
    vor = Voronoi(sites[:, 0:2])
    out = rpgen.RebuildVoronoi_merge(vor, sites, boundaries, c.generation_center, x_min, x_max, y_min, y_max,
                                     c.box_center, c.box_shape, c.merge_tol, True)
    voronoi_vertices, _, _, finite, boundary, nvertex, nvertices_in, nfinite, nboundary, nb_pts, nb_feat, ridges, nridge = out
    all_pts_2D, all_ridges, npt, npt_normal, npt_vtk = rpgen.RidgeMidQuarterPts(
        voronoi_vertices, nvertex, nvertices_in, ridges, finite, boundary, nfinite, nboundary, nb_pts, nb_feat)
    all_vertices_2D, max_wings, flattened, all_ridges = rpgen.VertexandRidgeinfo(
        all_pts_2D, all_ridges, npt, npt_normal, npt_vtk, nridge, "mesh2d", radii, c.generation_center,
        c.cellwallthickness_sparse, c.cellwallthickness_dense, path=".")
    return c, all_pts_2D, all_ridges, nvertex, nridge, all_vertices_2D, max_wings, flattened, npt_vtk


def test_fiber_end_layers_statistics():
    ends = L.fiber_end_layers(2000, 40, 0.05, fiber_length=0.5, fiber_length_cov=0.3, rng=1)
    assert len(ends) == 2000 and all(np.all((e >= 1) & (e <= 39)) and np.all(np.diff(e) > 0) for e in ends)
    stats = L.fiber_statistics(ends, 40, 0.05)
    # about 2 m of grains cut into fibers of 0.5 mm: about 4 fibers per grain; the fibers cut by the faces are shorter
    assert 3.0 < stats["nfibers"]/2000 < 5.5
    assert 0.3 < stats["mean fiber length (cut by the faces)"] < 0.5
    assert stats["longest fiber"] <= 2.0
    # the ends are spread over the grid lines, not on common planes
    counts = np.bincount(np.concatenate(ends), minlength=41)[1:40]
    assert counts.min() > 0 and counts.max() < 3*counts.mean()
    with pytest.raises(ValueError):
        L.fiber_end_layers(5, 10, 0.1, fiber_length=0)
    same = L.fiber_end_layers(50, 20, 0.1, 0.5, rng=7)
    again = L.fiber_end_layers(50, 20, 0.1, 0.5, rng=7)
    assert all(np.array_equal(a, b) for a, b in zip(same, again))


def test_lattice_with_regular_ends_matches_regular_extrusion(tmp_path, monkeypatch):
    # fiber ends on every segment interface of every grain: the same beams and connectors as the regular extrusion
    monkeypatch.chdir(tmp_path)
    c, all_pts_2D, all_ridges, nvertex, nridge, all_vertices_2D, max_wings, flattened, npt_vtk = mesh2d(SMALL)
    nelem = 2*c.nsegments
    nhalf = 2*nelem + 1
    theta = np.linspace(c.theta_min, c.theta_max, nhalf)
    layer_points = rpgen._rotate_layers(all_pts_2D[:npt_vtk], theta, np.linspace(c.z_min, c.z_max, nhalf),
                                        c.generation_center)[:, :, 0:2]
    ends = [np.arange(2, nelem, 2)]*nvertex
    lat = L.build_lattice(all_ridges, nvertex, nridge, all_vertices_2D, max_wings, flattened, layer_points,
                          c.z_min, c.z_max, nelem, c.long_connector_ratio, ends)
    regular = ringspy.generate(SMALL, geoName="regular")
    header = (regular.output_dir / "regular-mesh.txt").read_text().splitlines()
    nbot, nreg, ntop, nlong = [int(header[i]) for i in (2, 4, 6, 8)]
    assert (len(lat.connector_t_bot), len(lat.connector_t_reg), len(lat.connector_t_top), len(lat.connector_l)) == \
        (nbot, nreg, ntop, nlong)
    assert lat.nbeamElem == regular.nbeamElem and lat.npatch == nvertex*c.nsegments
    assert len(lat.nodes) == nvertex*(2*c.nsegments*2 + c.nsegments)  # 5 nodes per segment of each grain
    mesh = np.loadtxt(regular.output_dir / "regular-mesh.txt", skiprows=10)
    # same heights and widths of the transverse connectors, and the same total mass
    np.testing.assert_allclose(np.sort(lat.Meshdata[:nbot + nreg + ntop, 12]), np.sort(mesh[:nbot + nreg + ntop, 12]))
    np.testing.assert_allclose(np.sort(lat.Meshdata[:, 11]), np.sort(mesh[:, 11]))
    mass = lambda m: np.sum(m[:, 11]*m[:, 12]*np.linalg.norm(m[:, 5:8] - m[:, 8:11], axis=1))  # noqa: E731
    assert mass(lat.Meshdata) == pytest.approx(mass(mesh), rel=1e-4)  # the regular extrusion shifts the layers above each joint by its length
    # the fiber-end nodes sit connector_l_length apart, symmetric about the grid line
    pairs = lat.connector_l - 1
    gap = lat.nodes[pairs[:, 1], 2] - lat.nodes[pairs[:, 0], 2]
    np.testing.assert_allclose(gap, lat.connector_l_length)
    np.testing.assert_allclose((lat.nodes[pairs[:, 1], 2] + lat.nodes[pairs[:, 0], 2])/2 % lat.h, 0, atol=1e-9)


def test_lattice_connectors_cover_every_wall_on_every_grid_line():
    c, all_pts_2D, all_ridges, nvertex, nridge, all_vertices_2D, max_wings, flattened, npt_vtk = mesh2d(SMALL)
    nelem = 6
    layer_points = np.repeat(all_pts_2D[np.newaxis, :npt_vtk, :], 2*nelem + 1, axis=0)
    ends = L.fiber_end_layers(nvertex, nelem, (c.z_max - c.z_min)/nelem, 0.3, 0.5, rng=2)
    lat = L.build_lattice(all_ridges, nvertex, nridge, all_vertices_2D, max_wings, flattened, layer_points,
                          c.z_min, c.z_max, nelem, 0.02, ends)
    # on every grid line, the connectors of every ridge add up to the element length (half at the faces)
    height = np.zeros((nelem + 1, nridge))
    for kind in range(3):
        rows, lines = lat.connector_t_ridge[kind], lat.connector_t_line[kind]
        heights = [len(lat.connector_t_bot), len(lat.connector_t_reg), len(lat.connector_t_top)]
        offset = sum(heights[:kind])
        np.add.at(height, (lines, rows), lat.Meshdata[offset:offset + heights[kind], 12])
    np.testing.assert_allclose(height[1:-1], lat.h)
    np.testing.assert_allclose(height[[0, -1]], lat.h/2)
    # each fiber end gives one longitudinal connector per wing of the grain, between its two nodes on the grid line
    nwings = all_vertices_2D[:, 2].astype(int)
    assert len(lat.connector_l) == sum(nwings[g]*len(e) for g, e in enumerate(ends))
    nodes = lat.connector_l - 1
    assert np.all(lat.node_grain[nodes[:, 0]] == lat.node_grain[nodes[:, 1]])
    assert np.all(lat.node_layer[nodes[:, 0]] == lat.node_layer[nodes[:, 1]])
    assert np.all(lat.node_layer[nodes[:, 0]] % 2 == 0)
    # every beam element joins three consecutive nodes of one grain; every grain has nelem elements
    beams = lat.beam_connectivity - 1
    assert np.all(lat.node_grain[beams[:, 0]] == lat.node_grain[beams[:, 2]])
    assert np.all(lat.node_layer[beams[:, 2]] - lat.node_layer[beams[:, 0]] == 2)
    assert lat.nbeamElem == nvertex*nelem
    assert lat.element_nodes.shape == (nvertex, nelem, 3)
    np.testing.assert_array_equal(np.sort(lat.element_nodes.reshape(-1, 3), axis=0), np.sort(beams, axis=0))
    with pytest.raises(ValueError):
        L.build_lattice(all_ridges, nvertex, nridge, all_vertices_2D, max_wings, flattened, layer_points,
                        c.z_min, c.z_max, nelem, 0.02, [np.array([nelem])]*nvertex)


@pytest.mark.parametrize("extra", [dict(), dict(precrackFlag=True, nsegments=1), dict(radial_growth_rule="binary_lloyd")])
def test_generate_with_fibers_writes_all_files(extra, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    config = SMALL.replace(geoName="fibers", fiber_length=0.3, fiber_length_cov=0.4, stlFlag=True, inpFlag=True, **extra)
    result = ringspy.generate(config)
    expected = ["fibers" + suffix for suffix in OUTPUT_SUFFIXES] + ["fibers.stl", "fibers.inp"]
    assert sorted(p.name for p in result.output_dir.iterdir()) == sorted(expected)
    log = (result.output_dir / "fibers-report.log").read_text()
    assert "FIBERS:" in log and "nfiber_ends:" in log and "fiber_length:                            0.3" in log
    assert result.nbeamElem > 0 and 0 < result.porosity < 1
    inp = (result.output_dir / "fibers.inp").read_text()
    assert "*Element, type=T3D2,elset=AllLongConns" in inp
    if extra.get("precrackFlag"):
        assert "AllPrecrackTConns" in inp and "nconnector_t_precrack" in log
    iga = (result.output_dir / "fibersIGA.txt").read_text().splitlines()
    sizes = [int(iga[i + 1]) for i, line in enumerate(iga) if line.startswith("Size of knot vectors")]
    if config.nsegments > 1:
        assert len(set(sizes)) > 1  # fibers of various lengths: knot vectors of various sizes
    assert max(sizes) == (int(iga[5]) - 1)//2 + 5  # a patch of n elements: n + 5 knots, 2n + 1 control points
    assert max(sizes) == int(iga[7]) + 5
    # the same seed gives the same lattice
    again = ringspy.generate(config, path="again")
    assert (again.output_dir / "fibers-mesh.txt").read_bytes() == (result.output_dir / "fibers-mesh.txt").read_bytes()
    assert (again.output_dir / "fibersIGA.txt").read_bytes() == (result.output_dir / "fibersIGA.txt").read_bytes()


def test_cells_of_vertices_and_taper():
    c, all_pts_2D, all_ridges, nvertex, nridge, all_vertices_2D, max_wings, flattened, npt_vtk = mesh2d(SMALL)
    from ringspy.pipeline import place_cells
    sites, _ = place_cells(c, np.random.default_rng(c.seed))
    _, _, _, _, boundaries, _, _ = rpgen.Clipping_Box(c.box_shape, c.box_center, c.box_size, True)
    vertex_cells = L.cells_of_vertices(all_pts_2D[:nvertex], all_ridges, sites, boundaries)
    assert len(vertex_cells) == nvertex
    corners_of_box = [v for v in range(nvertex) if np.isclose(abs(all_pts_2D[v, 0] - 0.55), 0.25) and np.isclose(abs(all_pts_2D[v, 1]), 0.25)]
    assert all(len(vertex_cells[v]) >= (1 if v in corners_of_box else 2) for v in range(nvertex))
    # the cells of an inside vertex are its nearest sites (a Voronoi vertex is equidistant to three sites)
    from scipy.spatial import cKDTree
    tree = cKDTree(sites[:, :2])
    kinds = np.zeros(nvertex, dtype=int)
    on_side = np.isclose(np.abs(all_pts_2D[:nvertex, 0] - 0.55), 0.25) | np.isclose(np.abs(all_pts_2D[:nvertex, 1]), 0.25)
    kinds[on_side] = 1
    inside = np.nonzero(~on_side)[0]
    nearest3 = tree.query(all_pts_2D[inside], k=3)[1]
    agree = [set(nearest3[i]) <= vertex_cells[v] for i, v in enumerate(inside) if len(vertex_cells[v]) == 3]
    assert np.mean(agree) > 0.95
    # taper: a cell ending at z = 0.5 shrinks to taper_ratio of its size there, and is untouched far away
    nl = 21
    z = np.linspace(0, 1, nl)
    V = np.repeat(all_pts_2D[np.newaxis, :nvertex], nl, axis=0)
    cell = next(cell for cell in sorted(set().union(*vertex_cells))
                if sum(1 for v in range(nvertex) if cell in vertex_cells[v] and not on_side[v]) >= 4
                and all(not on_side[v] for v in range(nvertex) if cell in vertex_cells[v]))
    verts = [v for v in range(nvertex) if cell in vertex_cells[v]]
    tapered = L.taper_vertices(V, z, vertex_cells, {cell: np.array([0.5])}, taper_length=0.2, taper_ratio=0.2,
                               kind=kinds, boundaries=boundaries)
    center = V[0, verts].mean(axis=0)
    size = lambda m: np.linalg.norm(tapered[m, verts] - center, axis=1).mean()  # noqa: E731
    assert size(10) == pytest.approx(0.2*size(0), rel=1e-6)   # at the end
    assert size(0) == pytest.approx(np.linalg.norm(V[0, verts] - center, axis=1).mean())
    np.testing.assert_allclose(tapered[[0, 5, 15, 20]], V[[0, 5, 15, 20]])  # outside the taper length
    assert size(9) < size(8) < size(7) < size(0)
    untouched = [v for v in range(nvertex) if v not in verts]
    np.testing.assert_allclose(tapered[:, untouched], V[:, untouched])
    # boundary points stay on their sides, corners of the box stay put
    side_cell = next(cell for cell in sorted(set().union(*vertex_cells)) if any(on_side[v] for v in range(nvertex) if cell in vertex_cells[v]))
    kinds2 = kinds.copy()
    corners = [v for v in range(nvertex) if np.isclose(abs(all_pts_2D[v, 0] - 0.55), 0.25) and np.isclose(abs(all_pts_2D[v, 1]), 0.25)]
    kinds2[corners] = 2
    tapered = L.taper_vertices(V, z, vertex_cells, {side_cell: np.array([0.5])}, 0.3, 0.1, kind=kinds2, boundaries=boundaries)
    rows = [v for v in range(nvertex) if side_cell in vertex_cells[v] and kinds2[v] == 1]
    moved = tapered[10, rows]
    assert np.all(np.isclose(np.abs(moved[:, 0] - 0.55), 0.25) | np.isclose(np.abs(moved[:, 1]), 0.25))
    np.testing.assert_allclose(tapered[:, corners], V[:, corners])
    with pytest.raises(ValueError):
        L.taper_vertices(V, z, vertex_cells, {cell: np.array([0.5])}, 0.2, taper_ratio=1.5)
    # cell ends and the fiber ends of the grains
    ends = L.cell_end_heights(sorted(set().union(*vertex_cells)), 0.0, 1.0, 0.4, 0.3, rng=3)
    grain_ends = L.grain_ends_from_cells(vertex_cells, ends, 0.0, 0.0625, 16)
    assert all(np.all((e >= 1) & (e <= 15)) for e in grain_ends)
    expected = sorted({int(round(zz/0.0625)) for cl in vertex_cells[0] for zz in ends[cl] if 0 < round(zz/0.0625) < 16})
    assert grain_ends[0].tolist() == expected


def test_generate_taper_moves_the_corners_of_the_ending_cells(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    blunt = ringspy.generate(SMALL, geoName="blunt", fiber_length=0.3, taper_ratio=1.0)
    tapered = ringspy.generate(SMALL, geoName="tapered", fiber_length=0.3, taper_ratio=0.1, path="t")
    assert tapered.nbeamElem == blunt.nbeamElem and tapered.mass != blunt.mass
    # same cells and fiber ends (same seed), different node positions near the cell ends only
    mesh_b = np.loadtxt(blunt.output_dir / "blunt-mesh.txt", skiprows=10)
    mesh_t = np.loadtxt(tapered.output_dir / "tapered-mesh.txt", skiprows=10)
    assert mesh_b.shape == mesh_t.shape and np.array_equal(mesh_b[:, 0:2], mesh_t[:, 0:2])
    header = (blunt.output_dir / "blunt-mesh.txt").read_text().splitlines()
    ntransverse = sum(int(header[i]) for i in (2, 4, 6))
    lengths_b = np.linalg.norm(mesh_b[:ntransverse, 5:8] - mesh_b[:ntransverse, 8:11], axis=1)
    lengths_t = np.linalg.norm(mesh_t[:ntransverse, 5:8] - mesh_t[:ntransverse, 8:11], axis=1)
    assert np.percentile(lengths_t, 5) < 0.5*np.percentile(lengths_b, 5)   # walls of tapering cells shrink
    assert np.isclose(lengths_t, lengths_b).mean() > 0.1                    # walls far from any cell end are as before
    log = (tapered.output_dir / "tapered-report.log").read_text()
    assert "taper_length:                            0.045" in log and "ncell_ends:" in log
    mesh = np.loadtxt(tapered.output_dir / "tapered-mesh.txt", skiprows=10)
    assert np.linalg.norm(mesh[:, 5:8] - mesh[:, 8:11], axis=1).min() > 0


def test_fiber_config_checks():
    with pytest.raises(ValueError, match="fiber_length"):
        ringspy.Config(fiber_length=0)
    with pytest.raises(ValueError, match="fiber_length_cov"):
        ringspy.Config(fiber_length_cov=-1)
    with pytest.raises(ValueError, match="taper_ratio"):
        ringspy.Config(taper_ratio=2)
    with pytest.raises(ValueError, match="taper_length"):
        ringspy.Config(taper_length=0)
    assert ringspy.Config(fiber_length=2.0).resolved().taper_length == pytest.approx(0.3)
    assert ringspy.Config().resolved().taper_length is None
    config = ringspy.load_config(EXAMPLES_DIR / "wood_cube_fibers.toml")
    assert config.fiber_length is not None and config.fiber_length > 0
