# -*- coding: utf-8 -*-
"""
Clipping with any simple polygon (box_shape = 'polygon'), also non-convex ones
"""

import re

import numpy as np
import pytest
from scipy.spatial import Voronoi

import ringspy
import ringspy.MeshGenTools as rpgen
from helpers import inside_polygon


def star(n=5, r1=1.0, r2=0.45):
    angles = np.pi/2 + np.arange(2*n)*np.pi/n
    radii = np.where(np.arange(2*n) % 2 == 0, r1, r2)
    return np.column_stack((radii*np.cos(angles), radii*np.sin(angles))).tolist()


DOGBONE = [[-1.2, -0.4], [-0.6, -0.4], [-0.4, -0.2], [0.4, -0.2], [0.6, -0.4], [1.2, -0.4],
           [1.2, 0.4], [0.6, 0.4], [0.4, 0.2], [-0.4, 0.2], [-0.6, 0.4], [-1.2, 0.4]]
POLYGONS = {"star": star(), "dogbone": DOGBONE, "star_clockwise": star()[::-1]}
# a square with a V-notch from the left side to the centre, narrower than merge_tol over its last fifth
VNOTCH = [[-1, -1], [1, -1], [1, 1], [-1, 1], [-1, 0.05], [0, 0], [-1, -0.05]]


def clip_sites(sites, polygon, merge, flag=True, shape="polygon", box_size=1.0):
    vor = Voronoi(sites[:, 0:2])
    x_min, x_max, y_min, y_max, boundaries, corners, _ = \
        rpgen.Clipping_Box(shape, (0, 0), box_size, flag, polygon=polygon)
    if merge:
        out = rpgen.RebuildVoronoi_merge(vor, sites, boundaries, (0, 0), x_min, x_max, y_min, y_max, (0, 0), shape,
                                         0.01, flag)
    else:
        out = rpgen.RebuildVoronoi(vor, sites, boundaries, (0, 0), x_min, x_max, y_min, y_max, (0, 0), shape, flag)
    return out[0], out[11], corners


def clip(polygon, merge, npts=2000, span=4.0, seed=1, shape="polygon", box_size=1.0, flag=True):
    rng = np.random.default_rng(seed)
    sites = np.hstack([rng.random((npts, 2))*span - span/2, np.full((npts, 1), 0.02)])
    return clip_sites(sites, polygon, merge, flag, shape, box_size)


def walls_leaving(vertices, ridges, polygon, tol=1e-9):
    """Number of walls that cross an edge of the polygon (from one side to the other) or lie outside it"""
    def side(o, s, x): # side of the line o-s that x is on; 0 within tol of the line (walls ending on an edge)
        d = s - o
        distance = (d[..., 0]*(x[..., 1] - o[..., 1]) - d[..., 1]*(x[..., 0] - o[..., 0]))/np.linalg.norm(d, axis=-1)
        return np.where(np.abs(distance) <= tol, 0, np.sign(distance))
    p, q = vertices[ridges[:, 0]][:, np.newaxis], vertices[ridges[:, 1]][:, np.newaxis]
    a = np.asarray(polygon, dtype=float)
    b = np.roll(a, -1, axis=0)
    crossing = (side(p, q, a)*side(p, q, b) < 0) & (side(a, b, p)*side(a, b, q) < 0)
    outside = ~inside_polygon(vertices[ridges].mean(axis=1), polygon, tol=1e-9)
    return int(np.count_nonzero(crossing.any(axis=1) | outside))


def boundary_length(vertices, ridges, polygon):
    """Total length of the walls lying along the boundary: the perimeter when they follow it once, with no
    overlapping walls"""
    on = distance_to_boundary(vertices, polygon) < 1e-9
    along = on[ridges[:, 0]] & on[ridges[:, 1]] & \
        (distance_to_boundary(vertices[ridges].mean(axis=1), polygon) < 1e-9)
    return np.linalg.norm(vertices[ridges[along, 0]] - vertices[ridges[along, 1]], axis=1).sum()


def perimeter(polygon):
    a = np.asarray(polygon, dtype=float)
    return np.linalg.norm(np.roll(a, -1, axis=0) - a, axis=1).sum()


def distance_to_boundary(points, polygon):
    a = np.asarray(polygon, dtype=float)
    e = np.roll(a, -1, axis=0) - a
    rel = points[:, np.newaxis, :] - a[np.newaxis]
    u = np.clip(np.einsum('nmk,mk->nm', rel, e)/np.einsum('mk,mk->m', e, e), 0, 1)
    return np.linalg.norm(rel - u[:, :, np.newaxis]*e[np.newaxis], axis=2).min(axis=1)


def test_polygon_vertices_counterclockwise_and_closed():
    square = [[0, 0], [1, 0], [1, 1], [0, 1]]
    np.testing.assert_array_equal(rpgen.polygon_vertices(square + [[0, 0]]), square)
    np.testing.assert_array_equal(rpgen.polygon_vertices(square[::-1]), np.array(square[::-1])[::-1])


@pytest.mark.parametrize("vertices,message", [
    ([[0, 0], [1, 0]], "at least 3"),
    ([[0, 0], [1, 0], [2, 0]], "no area"),
    ([[0, 0], [1, 1], [1, 0], [0, 1]], "not simple"),
    ([[0, 0], [2, 0], [2, 2], [1, 0], [0, 2]], "not simple"),
    ([[0, 0], [1, 0], [np.nan, 1]], "pairs of numbers"),
])
def test_polygon_vertices_rejects(vertices, message):
    with pytest.raises(ValueError, match=message):
        rpgen.polygon_vertices(vertices)


@pytest.mark.parametrize("merge", [False, True], ids=["no_merge", "merge"])
@pytest.mark.parametrize("name", list(POLYGONS))
def test_clipped_mesh_inside_polygon_and_closed(name, merge):
    polygon = np.array(POLYGONS[name])
    vertices, ridges, corners = clip(polygon, merge)
    assert np.all(np.isfinite(vertices))
    assert inside_polygon(vertices, polygon, tol=1e-6).all()
    # every wall lies in the polygon (also near its reflex corners)
    assert inside_polygon(vertices[ridges].mean(axis=1), polygon, tol=1e-6).all()
    # every vertex closes at least two walls, and the polygon's corners are vertices of the mesh
    assert (np.bincount(ridges.ravel(), minlength=len(vertices)) >= 2).all()
    for corner in corners:
        assert np.linalg.norm(vertices - corner, axis=1).min() < 1e-9


def test_polygon_vertices_many_vertices():
    # a digitised outline: thousands of vertices are checked at once, not edge pair by edge pair
    angles = np.linspace(0, 2*np.pi, 3000, endpoint=False)
    circle = np.column_stack((np.cos(angles), np.sin(angles)))
    np.testing.assert_array_equal(rpgen.polygon_vertices(circle), circle)
    spike = circle.copy()
    spike[1500] = [1.5, 0.0] # reaches across the circle
    with pytest.raises(ValueError, match="not simple: edges 0 and 1499"):
        rpgen.polygon_vertices(spike)


@pytest.mark.parametrize("flag", [True, False], ids=["boundary", "no_boundary"])
@pytest.mark.parametrize("seed", [1, 2, 3])
def test_merge_does_not_bridge_a_narrow_notch(seed, flag):
    # near its tip the notch is narrower than merge_tol: a merge must not join its two sides
    vertices, ridges, corners = clip(VNOTCH, True, npts=20000, span=2.2, seed=seed, flag=flag)
    assert walls_leaving(vertices, ridges, VNOTCH) == 0
    if flag:
        for corner in corners:
            assert np.linalg.norm(vertices - corner, axis=1).min() < 1e-12
        assert boundary_length(vertices, ridges, VNOTCH) == pytest.approx(perimeter(VNOTCH), abs=1e-9)


@pytest.mark.parametrize("seed", [1, 2, 3])
def test_merge_keeps_short_polygon_edges_and_boundary(seed):
    # a U-notch drawn with edges shorter than merge_tol: the corners stay, no wall leaves the polygon, the walls
    # along the boundary follow it once (no wall lying along an edge over other boundary vertices)
    arc = np.linspace(0, -np.pi, 25) # the bottom of the notch, from right to left, edges of 0.0065
    notch = np.column_stack((0.05*np.cos(arc), 0.2 + 0.05*np.sin(arc))).tolist()
    polygon = [[-0.5, -0.5], [0.5, -0.5], [0.5, 0.5], [0.05, 0.5]] + notch + [[-0.05, 0.5], [-0.5, 0.5]]
    vertices, ridges, corners = clip(polygon, True, npts=6000, span=1.1, seed=seed)
    assert walls_leaving(vertices, ridges, polygon) == 0
    for corner in corners:
        assert np.linalg.norm(vertices - corner, axis=1).min() < 1e-12
    assert boundary_length(vertices, ridges, polygon) == pytest.approx(perimeter(polygon), abs=1e-9)
    # and the merge did merge: fewer vertices than without it
    assert len(vertices) < len(clip(polygon, False, npts=6000, span=1.1, seed=seed)[0])


def test_polygon_vertices_near_duplicates_and_edge_numbers():
    # a circle from linspace repeats its first vertex up to 1e-16
    angles = np.linspace(0, 2*np.pi, 73)
    assert len(rpgen.polygon_vertices(np.column_stack((np.cos(angles), np.sin(angles))))) == 72
    # the edges in the message are numbered on the list as given (here with a repeated vertex)
    with pytest.raises(ValueError, match="edges 1 and 3"):
        rpgen.polygon_vertices([[0, 0], [0, 0], [2, 2], [2, 0], [0, 2]])


@pytest.mark.parametrize("merge", [False, True], ids=["no_merge", "merge"])
@pytest.mark.parametrize("name", ["star", "dogbone"])
def test_no_loose_walls_without_boundary(name, merge):
    # without boundaryFlag there are no walls along the boundary, so a wall with both ends on it would be loose
    polygon = POLYGONS[name]
    vertices, ridges, _ = clip(np.array(polygon), merge, flag=False)
    on = distance_to_boundary(vertices, polygon) < 1e-9
    assert on.any()
    assert not np.any(on[ridges[:, 0]] & on[ridges[:, 1]])
    assert walls_leaving(vertices, ridges, polygon) == 0


@pytest.mark.parametrize("flag", [True, False], ids=["boundary", "no_boundary"])
@pytest.mark.parametrize("merge", [False, True], ids=["no_merge", "merge"])
def test_lattice_with_reflex_corner_on_a_voronoi_vertex(merge, flag):
    # sites of a triangular lattice: their Voronoi walls (a honeycomb) meet at the reflex corner of the L and
    # end on its edges at the same points; walls ending at one place share one vertex
    spacing = 0.05
    i, j = np.meshgrid(np.arange(-40, 41), np.arange(-40, 41))
    x, y = spacing*(i + 0.5*j), spacing*np.sqrt(3)/2*j
    keep = (np.abs(x) < 1.2) & (np.abs(y) < 1.2)
    sites = np.column_stack((x[keep], y[keep], np.full(np.count_nonzero(keep), 0.02)))
    c0 = np.array([spacing/2, spacing*np.sqrt(3)/6]) # a Voronoi vertex
    L = (c0 + np.array([[-0.5, -0.5], [0.5, -0.5], [0.5, 0], [0, 0], [0, 0.5], [-0.5, 0.5]])).tolist()
    vertices, ridges, corners = clip_sites(sites, L, merge, flag)
    lengths = np.linalg.norm(vertices[ridges[:, 0]] - vertices[ridges[:, 1]], axis=1)
    assert lengths.min() > 1e-9
    assert len(np.unique(vertices.round(9), axis=0)) == len(vertices)
    assert walls_leaving(vertices, ridges, L) == 0
    if flag:
        assert (np.bincount(ridges.ravel(), minlength=len(vertices)) >= 2).all()
        for corner in corners:
            assert np.linalg.norm(vertices - corner, axis=1).min() < 1e-12


def test_polygon_square_same_cells_as_square_box():
    square = [[-0.5, -0.5], [0.5, -0.5], [0.5, 0.5], [-0.5, 0.5]]
    box_vertices, box_ridges, _ = clip(None, False, shape="square")
    vertices, ridges, _ = clip(square, False)
    # the same vertices inside the box; the polygon clip also keeps the short walls that cut across a corner
    # of the box with both ends outside it
    inner = lambda v: v[(np.abs(v) < 0.5 - 1e-9).all(axis=1)]
    np.testing.assert_allclose(np.sort(inner(box_vertices), axis=0), np.sort(inner(vertices), axis=0))
    assert len(ridges) >= len(box_ridges)


def test_walls_between_two_boundary_points():
    # a wall can cross a corner of the polygon with both its ends outside: it becomes a wall between two points
    # on the boundary
    square = [[-0.5, -0.5], [0.5, -0.5], [0.5, 0.5], [-0.5, 0.5]]
    rng = np.random.default_rng(1)
    sites = np.hstack([rng.random((2000, 2))*4 - 2, np.full((2000, 1), 0.02)])
    *_, boundaries, _, _ = rpgen.Clipping_Box("polygon", (0, 0), 1, True, polygon=square)
    inside = rpgen._inside_clip_box("polygon", boundaries, 0, 0, 0, 0, (0, 0))
    vertices_in, points, _, ncut, _, _, cut_ridges = \
        rpgen._clip_voronoi_polygon(Voronoi(sites[:, 0:2]), sites, boundaries, (0, 0), inside, True)
    both = cut_ridges[np.all(cut_ridges >= len(vertices_in), axis=1)]
    assert len(both) > 0
    ends = points[both - len(vertices_in)]
    assert np.all(np.isclose(np.abs(ends), 0.5).any(axis=2))  # both ends on the boundary


def test_config_polygon():
    config = ringspy.Config(box_shape="polygon", box_polygon=DOGBONE)
    assert config.resolved().z_max == pytest.approx(2.4)  # largest width of the polygon
    assert ringspy.Config(box_shape="polygon", box_polygon=DOGBONE, z_max=1.0).resolved().z_max == 1.0
    with pytest.raises(ValueError, match="needs box_polygon"):
        ringspy.Config(box_shape="polygon")
    with pytest.raises(ValueError, match="only used with box_shape"):
        ringspy.Config(box_polygon=DOGBONE)
    with pytest.raises(ValueError, match="box_polygon: the polygon is not simple"):
        ringspy.Config(box_shape="polygon", box_polygon=[[0, 0], [1, 1], [1, 0], [0, 1]])
    with pytest.raises(ValueError, match="precrack_nodes"):
        ringspy.Config(box_shape="polygon", box_polygon=DOGBONE, precrackFlag=True)
    with pytest.raises(TypeError, match="box_polygon"):
        ringspy.Config(box_shape="polygon", box_polygon=[[0, 0, 1], [1, 0, 1], [0, 1, 1]])


def test_toml_polygon(tmp_path):
    toml = tmp_path / "poly.toml"
    toml.write_text('box_shape = "polygon"\nbox_polygon = [[0, 0], [2, 0], [2, 1], [0, 1]]\n')
    config = ringspy.load_config(toml)
    assert config.box_polygon == [(0, 0), (2, 0), (2, 1), (0, 1)]


def test_generate_polygon_volume_and_log():
    result = ringspy.generate(dict(geoName="dogbone_small", seed=3, r_max=1.6, nrings=2, iter_max=300,
                                   cellsize_sparse=0.04, cellsize_dense=0.02, box_shape="polygon",
                                   box_polygon=np.array(DOGBONE), z_max=0.5, stlFlag=False, inpFlag=False))
    polygon = np.array(DOGBONE)
    area = 0.5*abs(np.sum(rpgen.cross2D(polygon, np.roll(polygon, -1, axis=0))))
    assert result.bulk_volume == pytest.approx(area*0.5)
    assert 0 < result.porosity < 1
    log = (result.output_dir / "dogbone_small-report.log").read_text()
    assert re.search(r"^box_polygon:\s+\[\[-1\.2, -0\.4\]", log, re.M) # plain numbers, also from a NumPy array
