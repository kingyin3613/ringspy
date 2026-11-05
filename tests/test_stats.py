# -*- coding: utf-8 -*-
"""
Statistics of the cells, walls and growth rings (ringspy.stats)
"""

import csv

import numpy as np
import pytest

import ringspy
import ringspy.MeshGenTools as rpgen
from ringspy import stats

SMALL = dict(r_max=1.5, nrings=2, box_size=1.2, box_center=[0.9, 0], iter_max=300, seed=5, stlFlag=False, inpFlag=False)


def test_polygon_moments():
    # a 2 x 1 rectangle centred at (1, 0.5), rotated by 30 degrees
    rect = np.array([[-1, -0.5], [1, -0.5], [1, 0.5], [-1, 0.5]], dtype=float)
    rot = np.array([[np.cos(0.5), -np.sin(0.5)], [np.sin(0.5), np.cos(0.5)]])
    p = rect @ rot.T + [1, 0.5]
    area, centroid, cov = stats.polygon_moments(p)
    assert area == pytest.approx(2.0) and np.allclose(centroid, [1, 0.5])
    values, vectors = np.linalg.eigh(cov)
    np.testing.assert_allclose(sorted(values), [1/12, 4/12])  # second moments of a rectangle: w^2/12, h^2/12
    major = vectors[:, 1]
    assert abs(abs(np.dot(major, [np.cos(0.5), np.sin(0.5)])) - 1) < 1e-9
    # the orientation does not depend on the direction around the polygon
    area2, _, cov2 = stats.polygon_moments(p[::-1])
    assert area2 == pytest.approx(2.0) and np.allclose(cov2, cov)


def test_ridge_cells_and_polygons_on_a_square_box():
    # four sites in a unit box: four square cells, four interior walls, eight walls along the box
    sites = np.array([[0.25, 0.25], [0.75, 0.25], [0.25, 0.75], [0.75, 0.75]])
    vertices = np.array([[0.5, 0.5], [0.5, 0.0], [0.5, 1.0], [0.0, 0.5], [1.0, 0.5],
                         [0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [1.0, 1.0]], dtype=float)
    ridges = np.array([[0, 1], [0, 2], [0, 3], [0, 4],
                       [5, 1], [1, 6], [6, 4], [4, 8], [8, 2], [2, 7], [7, 3], [3, 5]])
    _, _, _, _, boundaries, _, _ = rpgen.Clipping_Box("square", (0.5, 0.5), 1.0, True)
    cells = stats.ridge_cells(sites, vertices, ridges, boundaries)
    assert np.all(cells[:4] >= 0) and np.all(cells[4:, 1] == -1)
    assert sorted(cells[0].tolist()) == [0, 1]  # the wall (0.5, 0.5)-(0.5, 0) separates the two bottom cells
    polygons, _ = stats.cell_polygons(sites, vertices, ridges, boundaries)
    for poly in polygons:
        area, _, _ = stats.polygon_moments(vertices[poly])
        assert len(poly) == 4 and area == pytest.approx(0.25)


@pytest.fixture(scope="module")
def result(tmp_path_factory):
    import os
    workdir = tmp_path_factory.mktemp("stats")
    cwd = os.getcwd()
    os.chdir(workdir)
    try:
        res = ringspy.generate(dict(SMALL, geoName="st", statsFlag=True))
    finally:
        os.chdir(cwd)
    res.output_dir = workdir / res.output_dir
    return res


def test_cell_statistics(result):
    c = result.cell_stats
    n = result.nParticles
    assert set(c) == set(stats.CELL_FIELDS) and all(len(v) == n for v in c.values())
    has = np.isfinite(c["area"])
    assert 0 < has.sum() < n  # sites outside the box have no polygon
    assert np.all(c["boundary"][~has] == 1) and np.all(c["nwalls"][has] >= 3)
    interior = has & (c["boundary"] == 0)
    assert interior.sum() > 50
    # the cells tile the box: their areas add up to the box area
    assert np.nansum(c["area"]) == pytest.approx(1.2**2, rel=5e-3)
    # interior cells: as many neighbours as walls, a few sides each, the usual shapes
    np.testing.assert_array_equal(c["nneighbours"][interior], c["nwalls"][interior])
    assert 4.5 < c["nneighbours"][interior].mean() < 6.5
    assert np.all(c["aspect_ratio"][interior] >= 1) and 1.1 < np.median(c["aspect_ratio"][interior]) < 2.5
    assert np.all(np.abs(c["orientation"][has]) <= np.pi/2 + 1e-12)
    assert np.all((c["radial_elongation"][has] >= 0) & (c["radial_elongation"][has] <= 1))
    np.testing.assert_allclose(c["equivalent_diameter"][has], np.sqrt(4*c["area"][has]/np.pi))
    assert np.all((c["solid_fraction"][has] > 0) & (c["solid_fraction"][has] <= 1))  # a tiny merged cell can be all wall
    # labels from the site positions
    band = np.searchsorted(result.radii, np.linalg.norm(result.sites[:, :2], axis=1), side="right") - 1
    np.testing.assert_array_equal(c["band"], band)
    np.testing.assert_array_equal(c["dense"], band % 2 == 0)
    assert np.all(c["type"] == 0)
    # thicknesses: the dense cells have the thicker walls
    assert np.nanmean(c["mean_thickness"][interior & (c["dense"] == 1)]) > np.nanmean(c["mean_thickness"][interior & (c["dense"] == 0)])


def test_wall_and_ring_statistics(result):
    w = result.wall_stats
    assert set(w) == set(stats.WALL_FIELDS) and all(len(v) == result.nridge for v in w.values())
    np.testing.assert_allclose(w["length"], np.linalg.norm(result.vertices[result.ridges[:, 1]] - result.vertices[result.ridges[:, 0]], axis=1))
    np.testing.assert_array_equal(w["thickness"], result.vertex_thickness[result.ridges[:, 1]])
    np.testing.assert_array_equal(w["boundary"], result.ridge_boundary)
    assert np.all(w["cell2"][w["boundary"] == 1] == -1) and np.all(w["cell2"][w["boundary"] == 0] >= 0)
    # the wall areas of the cells add up to the wall areas of the walls
    c = result.cell_stats
    assert np.nansum(c["wall_area"]) == pytest.approx(np.sum(w["length"]*w["thickness"]), rel=1e-9)
    r = result.ring_stats
    assert set(r) == set(stats.RING_FIELDS) and len(r["band"]) == len(np.unique(c["band"]))
    assert np.sum(r["ncells"]) == result.nParticles
    assert np.nansum(r["area"]) == pytest.approx(np.nansum(c["area"]))
    valid = r["area"] > 0
    assert np.all((r["porosity_2d"][valid] > 0) & (r["porosity_2d"][valid] < 1))
    # the dense bands are less porous
    dense, sparse = valid & (r["dense"] == 1), valid & (r["dense"] == 0)
    assert r["porosity_2d"][dense].mean() < r["porosity_2d"][sparse].mean()
    s = result.summary()
    assert set(s) == set(stats.SUMMARY_FIELDS) and s["ncells"] == int(np.isfinite(c["area"]).sum())
    assert s["porosity_2d_dense"] < s["porosity_2d"] < s["porosity_2d_sparse"]
    assert s["mean_wall_thickness"] == pytest.approx(np.mean(w["thickness"]))


def test_statistics_files(result):
    out = result.output_dir
    for suffix, table, fields in [("-cells.csv", result.cell_stats, stats.CELL_FIELDS),
                                  ("-rings.csv", result.ring_stats, stats.RING_FIELDS),
                                  ("-walls.csv", result.wall_stats, stats.WALL_FIELDS)]:
        with open(out / ("st" + suffix), newline="") as f:
            rows = list(csv.reader(f))
        assert rows[0] == fields and len(rows) == 1 + len(table[fields[0]])
        # numbers read back (NaN as empty cells)
        col = fields.index("area") if "area" in fields else fields.index("length")
        values = np.array([float(r[col]) if r[col] else np.nan for r in rows[1:]])
        np.testing.assert_allclose(values, table[fields[col]])
    log = (out / "st-report.log").read_text()
    assert "Cell Statistics File:" in log and "st-walls.csv" in log


def test_statistics_without_the_flag_and_with_fibers_and_hardwood():
    res = ringspy.generate(dict(SMALL, geoName="nofile"))
    assert res.cell_stats is not None and not any(p.suffix == ".csv" for p in res.files)
    assert "Cell Statistics File" not in (res.output_dir / "nofile-report.log").read_text()
    oak = ringspy.generate(dict(SMALL, geoName="oak", fiber_length=0.5, vessel_fraction=0.2, vessel_size=0.1,
                                vessel_arrangement="ring", ray_spacing=0.3, statsFlag=True))
    c = oak.cell_stats
    vessels = (c["type"] == 1) & np.isfinite(c["area"]) & (c["boundary"] == 0)
    fibers = (c["type"] == 0) & np.isfinite(c["area"]) & (c["boundary"] == 0)
    assert vessels.any() and np.median(c["area"][vessels]) > 3*np.median(c["area"][fibers])
    # the rays are continuous files: an interior ray cell touches other ray cells
    w = oak.wall_stats
    rays = (c["type"] == 2) & np.isfinite(c["area"]) & (c["boundary"] == 0)
    ray_neighbours = np.zeros(len(c["type"]), dtype=int)
    both = (w["cell2"] >= 0)
    a, b = w["cell1"][both], w["cell2"][both]
    np.add.at(ray_neighbours, a, c["type"][b] == 2)
    np.add.at(ray_neighbours, b, c["type"][a] == 2)
    assert rays.sum() > 10 and np.mean(ray_neighbours[rays] >= 1) > 0.9 and np.mean(ray_neighbours[rays] >= 2) > 0.6
    assert (oak.output_dir / "oak-rings.csv").is_file()
