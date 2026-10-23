# -*- coding: utf-8 -*-
"""
Regression tests for fixed bugs, plus end-to-end checks of the generated files
(examples/ configurations with a few parameters changed).
"""

import re
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import pytest
from scipy.spatial import Voronoi

import ringspy
import ringspy.MeshGenTools as rpgen
from helpers import EXAMPLES_DIR, OUTPUT_SUFFIXES, inside_polygon

HONEYCOMB_HEX = EXAMPLES_DIR / "honeycomb_hex.toml"


def clip(shape, merge, npts, span, box_size=1.0, center=(0.0, 0.0), seed=1):
    """Clip a random Voronoi diagram, return the vertices, ridges and the clipping polygon"""
    rng = np.random.default_rng(seed)
    sites = rng.random((npts, 2)) * span - span / 2 + np.asarray(center)
    sites = np.hstack([sites, np.full((npts, 1), 0.05)])
    vor = Voronoi(sites[:, 0:2])
    x_min, x_max, y_min, y_max, boundaries, boundary_points, _ = \
        rpgen.Clipping_Box(shape, center, box_size, "on")
    if merge:
        out = rpgen.RebuildVoronoi_merge(vor, sites, boundaries, center, x_min, x_max, y_min, y_max,
                                         center, shape, 0.01, "on")
    else:
        out = rpgen.RebuildVoronoi(vor, sites, boundaries, center, x_min, x_max, y_min, y_max,
                                   center, shape, "on")
    voronoi_vertices, voronoi_ridges = out[0], out[11]
    return voronoi_vertices, voronoi_ridges, boundary_points


def test_import_without_pkg_resources():
    # pkg_resources is gone from setuptools >= 81
    code = "import sys; sys.modules['pkg_resources'] = None; import ringspy.MeshGenTools"
    subprocess.run([sys.executable, "-c", code], check=True)


def test_cross2D():
    assert rpgen.cross2D([1, 0], [0, 1]) == 1
    np.testing.assert_allclose(rpgen.cross2D([[1, 2], [3, 4]], [[5, 6], [7, 8]]), [-4, -4])


def test_hexagon_box_centered_on_box_center():
    box_center = (1.25, -0.5)
    x_min, x_max, y_min, y_max, _, boundary_points, _ = rpgen.Clipping_Box("hexagon", box_center, 1.0, "on")
    np.testing.assert_allclose(boundary_points.mean(axis=0), box_center)
    np.testing.assert_allclose([(x_min + x_max) / 2, (y_min + y_max) / 2], box_center)


def test_lloyd_honors_generation_center():
    generation_center = (5.0, -3.0)
    sites, _ = rpgen.CellPlacement_Binary_Lloyd("lloyd", "meshes", generation_center, 1.0, 0.0, 2,
                                                0.1, 0.3, 0.1, 0.04, 0.02, 5, 1000)
    np.testing.assert_allclose(sites[:, 0:2].mean(axis=0), generation_center, atol=0.1)


def test_lloyd_ring_narrower_than_one_cell(keep_global_random_state):
    # a ring narrower than one cell has no sub-rings inside: np.concatenate of no sub-rings crashed
    np.random.seed(0)
    sites, radii = rpgen.CellPlacement_Binary_Lloyd("lloyd", "meshes", (0.0, 0.0), 1.0, 0.0, 2,
                                                    0.1, 0.3, 0.01, 0.04, 0.02, 5, 1000)
    widths = np.diff(radii)
    assert (widths[[2, 4]] < 0.02).all()  # the dense rings
    assert np.all(np.isfinite(sites))
    # the cells on the outer edge of the narrow rings are placed
    distance = np.linalg.norm(sites, axis=1)
    for i in [2, 4]:
        assert np.count_nonzero(np.isclose(distance, radii[i+1], rtol=0, atol=1e-12)) == np.ceil(2*np.pi*radii[i+1]/(2*0.02))


@pytest.mark.parametrize("width_heart", [0.003, 0.006], ids=["1_perimeter_point", "2_perimeter_points"])
def test_lloyd_heart_narrower_than_one_cell(keep_global_random_state, width_heart):
    # a heart narrower than one cell has no cells inside: the Lloyd relaxation of its 1 or 2
    # perimeter points alone crashed in Voronoi (QhullError)
    np.random.seed(0)
    sites, radii = rpgen.CellPlacement_Binary_Lloyd("lloyd", "meshes", (0.0, 0.0), 1.0, 0.0, 2,
                                                    width_heart, 0.3, 0.1, 0.04, 0.02, 5, 1000)
    assert np.ceil(2*np.pi*radii[1]/(2*0.02)) <= 2
    assert np.all(np.isfinite(sites))
    Voronoi(sites[:, 0:2])


@pytest.mark.parametrize("merge", [False, True], ids=["no_merge", "merge"])
@pytest.mark.parametrize("shape", ["square", "triangle", "hexagon"])
@pytest.mark.parametrize("npts,span,box_size,center", [
    (1500, 6.0, 1.0, (0.0, 0.0)),    # sites cover the box
    (1500, 6.0, 1.0, (1.25, -0.5)),  # box away from the origin
    (300, 1.0, 4.0, (0.0, 0.0)),     # box larger than the site cloud (infinite ridges inside the box)
], ids=["covered", "off_origin", "box_larger_than_sites"])
def test_clipped_mesh_inside_box_and_closed(shape, merge, npts, span, box_size, center):
    # used to hang (triangle) or crash (merge with triangle/hexagon, box larger than the site cloud)
    vertices, ridges, polygon = clip(shape, merge, npts, span, box_size, center)
    assert np.all(np.isfinite(vertices))
    assert inside_polygon(vertices, polygon).all()
    # every vertex closes at least two cell walls: no dangling walls at the boundary
    degree = np.bincount(ridges.ravel(), minlength=len(vertices))
    assert (degree >= 2).all()


@pytest.mark.parametrize("boundaryFlag", ["on", "off"])
def test_merge_counts_vertices_with_only_cut_ridges(keep_global_random_state, boundaryFlag):
    # with the wood_cube cells of seed 35, a vertex near a box corner has its three neighbouring
    # vertices outside the box, so all its ridges are cut. It was not counted in nvertices_in.
    np.random.seed(35)
    sites, _ = rpgen.CellPlacement_Binary((0, 0), 4, 0, 4, 0.3, 0.7, 0.3, 0.02, 0.01, 500, 10**9)
    x_min, x_max, y_min, y_max, boundaries, _, _ = rpgen.Clipping_Box("square", (0, 0), 4.0, boundaryFlag)
    out = rpgen.RebuildVoronoi_merge(Voronoi(sites[:, 0:2]), sites, boundaries, (0, 0), x_min, x_max, y_min, y_max,
                                     (0, 0), "square", 0.01, boundaryFlag)
    voronoi_vertices, finite_ridges_new, nvertex, nvertices_in, nboundary_pts = out[0], out[3], out[5], out[6], out[9]
    inside = np.all(np.abs(voronoi_vertices) < 2 - 1e-6, axis=1)
    assert np.any(inside & (np.bincount(finite_ridges_new.ravel(), minlength=nvertex) == 0))
    assert nvertices_in + nboundary_pts == nvertex == len(voronoi_vertices)


def test_wood_cube_with_vertices_with_only_cut_ridges(keep_global_random_state):
    # ~7% of unseeded wood_cube runs crashed in GenerateBeamElement (npt_per_layer was nvertex - 1)
    np.random.seed(35)
    ringspy.generate(EXAMPLES_DIR / "wood_cube.toml", seed=None, stlFlag=False)  # cells from the global random state
    out = Path("meshes", "wood_cube")
    nvertex = int((out / "wood_cube-vertex.mesh").read_text().splitlines()[2])
    ridges = np.loadtxt(out / "wood_cube-ridge.mesh", skiprows=4, usecols=(0, 1), dtype=int)
    assert set(np.unique(ridges).tolist()) == set(range(nvertex))
    assert (out / "wood_cube.inp").is_file()


def test_read_saved_sites_searches_parent_folders(tmp_path):
    np.save(tmp_path / "cube_sites.npy", np.zeros((3, 3)))
    np.save(tmp_path / "cube_radii.npy", np.ones(4))
    (tmp_path / "a" / "b").mkdir(parents=True)
    sites, radii = rpgen.ReadSavedSites(tmp_path / "a" / "b", "cube_sites.npy")
    assert sites.shape == (3, 3) and radii.shape == (4,)


def test_read_saved_sites_missing_file_raises(tmp_path):
    # used to loop until a 180 s timeout on Windows, then call exit()
    start = time.time()
    with pytest.raises(FileNotFoundError):
        rpgen.ReadSavedSites(tmp_path, "missing_sites.npy")
    assert time.time() - start < 10


def test_unknown_box_shape_raises():
    with pytest.raises(ValueError):
        rpgen.Clipping_Box("circle", (0, 0), 1.0, "on")


def test_honeycomb_hexagon_off_origin_writes_to_path():
    box_center = (0.25, 0.1)
    ringspy.generate(HONEYCOMB_HEX, geoName="hex_off_origin", path="custom_out", box_center=box_center)
    out = Path("custom_out", "hex_off_origin")
    assert not Path("meshes").exists()
    for suffix in OUTPUT_SUFFIXES + [".stl", ".inp"]:
        assert (out / ("hex_off_origin" + suffix)).is_file(), suffix
    # binary STL: 80-byte header, triangle count, 50 bytes per triangle
    stl = (out / "hex_off_origin.stl").read_bytes()
    ntriangles = int.from_bytes(stl[80:84], "little")
    assert ntriangles > 0 and len(stl) == 84 + 50 * ntriangles
    # all cell vertices lie inside the hexagon centered at box_center
    vertices = np.loadtxt(out / "hex_off_origin-vertex.mesh", skiprows=6, usecols=(0, 1))
    _, _, _, _, _, polygon, _ = rpgen.Clipping_Box("hexagon", box_center, 1.0, "on")
    assert inside_polygon(vertices, polygon, tol=1e-6).all()


def parse_inp(path):
    """Element ids per element set, section element sets, and the AllElles range of an Abaqus input file"""
    elsets, sections, current, all_elles = {}, [], None, None
    lines = Path(path).read_text().splitlines()
    for i, line in enumerate(lines):
        if line.startswith("**"):
            continue
        if line.startswith("*"):
            current = None
            keyword = line.lower().replace(" ", "")
            name = re.search(r"elset=([^,]+)", keyword)
            if keyword.startswith("*element,"):
                current = elsets.setdefault(name.group(1), [])
            elif "section," in keyword:
                sections.append(name.group(1))
            elif keyword.startswith("*elset,elset=allelles"):
                all_elles = [int(x) for x in lines[i + 1].split(",")[:2]]
        elif current is not None and line.strip():
            current.append(int(line.split(",")[0]))
    return elsets, sections, all_elles


def test_precracked_inp_is_consistent():
    nsegments = 2
    # a crack across several cells: from 0.12 to 0.72 of the box width, just above the center
    ringspy.generate(HONEYCOMB_HEX, geoName="precrack", box_shape="square", stlFlag=False, precrackFlag=True,
                     nsegments=nsegments, precrack_nodes=[[-0.38, 0.013, 0.22, 0.013]])
    elsets, sections, all_elles = parse_inp(Path("meshes", "precrack", "precrack.inp"))
    # every section refers to a non-empty element set (VisualLongConns was missing before)
    for elset in sections:
        assert elsets.get(elset), elset
    ids = [i for members in elsets.values() for i in members]
    assert len(ids) == len(set(ids))
    # AllElles covers exactly the simulation elements (it counted precracked connectors twice)
    simulation_ids = {i for name, members in elsets.items() if name.startswith("all") for i in members}
    assert simulation_ids == set(range(all_elles[0], all_elles[1] + 1))
    # the precracked cell walls are cracked in every segment, for bot, reg and top connectors
    log = Path("meshes", "precrack", "precrack-report.log").read_text()
    ncracked_walls = int(re.search(r"nconnector_t_precrack:\s+(\d+)", log).group(1))
    assert ncracked_walls > 0
    assert len(elsets["allprecracktconns"]) == ncracked_walls * 3 * nsegments
    assert len(elsets["visualprecracktconns"]) == ncracked_walls * 3 * nsegments
