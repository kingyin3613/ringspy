# -*- coding: utf-8 -*-
"""
Statistics of the clipped 2D cell mesh: the cells as polygons (area, neighbours, shape and
orientation, wall area), the walls (length, thickness), and the growth rings (cell counts,
2D porosity), for parameter studies and for matching micrographs

author: king_yin3613
email: haoyin2022@u.northwestern.edu
"""

import csv
import math
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

from . import MeshGenTools as rpgen

CELL_FIELDS = ['site', 'x', 'y', 'type', 'band', 'ring', 'dense', 'boundary', 'nwalls', 'nneighbours', 'area',
               'perimeter', 'wall_area', 'solid_fraction', 'equivalent_diameter', 'aspect_ratio', 'orientation',
               'radial_elongation', 'mean_thickness']
RING_FIELDS = ['band', 'ring', 'dense', 'ncells', 'ncells_interior', 'area', 'wall_area', 'porosity_2d',
               'mean_cell_area', 'mean_equivalent_diameter', 'mean_aspect_ratio', 'mean_nneighbours', 'mean_thickness']
WALL_FIELDS = ['ridge', 'vertex1', 'vertex2', 'cell1', 'cell2', 'length', 'thickness', 'angle', 'ring', 'dense',
               'vessel', 'ray', 'boundary']
SUMMARY_FIELDS = ['ncells', 'mean_cell_area', 'mean_equivalent_diameter', 'mean_aspect_ratio', 'mean_nneighbours',
                  'porosity_2d', 'porosity_2d_sparse', 'porosity_2d_dense', 'mean_wall_length', 'mean_wall_thickness']


def ridge_cells(sites, vertices, ridges, boundaries=None, tol=None):
    """
    The cells (sites, by index) on the two sides of every wall (ridge) of a clipped Voronoi mesh: the
    two sites nearest to the midpoint of the wall, or one site and -1 for a wall along a side of the
    clipping box (boundaries, as Clipping_Box gives them). Array of shape (nridge, 2).
    """
    vertices = np.asarray(vertices, dtype=float)
    ridges = np.asarray(ridges)[:, 0:2].astype(int)
    sites = np.asarray(sites, dtype=float)[:, 0:2]
    mids = (vertices[ridges[:, 0]] + vertices[ridges[:, 1]])/2
    k = min(2, len(sites))
    _, index = cKDTree(sites).query(mids, k=k)
    index = np.asarray(index).reshape(len(mids), k)
    cells = -np.ones((len(mids), 2), dtype=int)
    cells[:, 0:k] = index
    if boundaries is not None and len(mids):
        starts = np.array([b[1][0] for b in boundaries], dtype=float)
        ends = np.array([b[1][1] for b in boundaries], dtype=float)
        segments = ends - starts
        lengths = np.linalg.norm(segments, axis=1)
        if tol is None:
            tol = 1e-6*lengths.max()
        along_side = np.zeros(len(mids), dtype=bool)
        for p0, seg, length in zip(starts, segments, lengths):
            rel = mids - p0
            u = np.clip((rel @ seg)/length**2, 0, 1)
            along_side |= np.linalg.norm(rel - u[:, np.newaxis]*seg, axis=1) <= tol
        cells[along_side, 1] = -1
    return cells


def cell_polygons(sites, vertices, ridges, boundaries=None):
    """
    The cells of a clipped Voronoi mesh as polygons: for every site, the indices of its vertices in
    order around the cell (an empty array for a site without walls in the mesh), and the cells of
    every ridge (see ridge_cells). The vertices of a cell are the ends of its walls, ordered by angle
    around their mean (the cells are convex).
    """
    vertices = np.asarray(vertices, dtype=float)
    ridges = np.asarray(ridges)[:, 0:2].astype(int)
    cells = ridge_cells(sites, vertices, ridges, boundaries)
    members = [[] for _ in range(len(sites))]
    for r in range(len(ridges)):
        for cell in cells[r]:
            if cell >= 0:
                members[cell].append(r)
    polygons = []
    for walls in members:
        if not walls:
            polygons.append(np.zeros(0, dtype=int))
            continue
        idx = np.unique(ridges[walls].reshape(-1))
        p = vertices[idx]
        center = p.mean(axis=0)
        order = np.argsort(np.arctan2(p[:, 1] - center[1], p[:, 0] - center[0]))
        polygons.append(idx[order])
    return polygons, cells


def polygon_moments(p):
    """Area, centroid and covariance matrix (second central moments per unit area) of a polygon (counter-clockwise)"""
    x, y = p[:, 0], p[:, 1]
    x1, y1 = np.roll(x, -1), np.roll(y, -1)
    cross = x*y1 - x1*y
    area = 0.5*np.sum(cross)
    if abs(area) < 1e-300:
        return 0.0, p.mean(axis=0), np.zeros((2, 2))
    cx = np.sum((x + x1)*cross)/(6*area)
    cy = np.sum((y + y1)*cross)/(6*area)
    ixx = np.sum((y*y + y*y1 + y1*y1)*cross)/12  # second moments about the origin
    iyy = np.sum((x*x + x*x1 + x1*x1)*cross)/12
    ixy = np.sum((x*y1 + 2*x*y + 2*x1*y1 + x1*y)*cross)/24
    cov = np.array([[iyy/area - cx*cx, ixy/area - cx*cy], [ixy/area - cx*cy, ixx/area - cy*cy]])
    return abs(area), np.array([cx, cy]), cov


def cell_statistics(result, boundaries=None):
    """
    One value per cell site of a Result (dict of arrays, keys CELL_FIELDS): site index, coordinates,
    type (0 fiber, 1 vessel, 2 ray), band, growth ring and dense flag (from the ring coordinate of
    the site), boundary (True for a cell cut by the clipping box, or without a closed polygon), the
    number of walls and of neighbouring cells, area, perimeter, wall area (half of every shared wall,
    the whole of a wall along the box), solid fraction (wall area / area), equivalent diameter, aspect
    ratio and orientation (radians, of the major axis, from the second moments of the polygon),
    radial elongation (|cos| of the angle between the major axis and the direction from the
    generation center, 1 for a cell elongated radially) and the mean thickness of its walls. Cells
    without a closed polygon (sites outside the box) have NaN for the geometry; their wall area and
    mean thickness are kept whenever they own walls, so that the wall areas of all cells add up.
    """
    sites = np.asarray(result.sites, dtype=float)
    vertices = np.asarray(result.vertices, dtype=float)
    ridges = np.asarray(result.ridges, dtype=int)
    polygons, cells = cell_polygons(sites, vertices, ridges, boundaries)
    n = len(sites)
    thickness = np.asarray(result.vertex_thickness, dtype=float)[ridges[:, 1]]
    length = np.linalg.norm(vertices[ridges[:, 1]] - vertices[ridges[:, 0]], axis=1)
    shared = cells[:, 1] >= 0
    out = {name: np.full(n, np.nan) for name in CELL_FIELDS}
    out['site'] = np.arange(n)
    out['x'], out['y'] = sites[:, 0], sites[:, 1]
    out['type'] = np.asarray(result.site_type, dtype=int) if len(result.site_type) == n else np.zeros(n, dtype=int)
    shape = result.ring_shape if result.ring_shape is not None else rpgen.RingShape(result.config.generation_center)
    band = np.searchsorted(result.radii, shape.rho(sites[:, 0:2]), side='right') - 1
    ring, dense = rpgen.RingLabels(band)
    out['band'], out['ring'], out['dense'] = band, ring, dense.astype(int)
    out['nwalls'] = np.zeros(n, dtype=int)
    out['nneighbours'] = np.zeros(n, dtype=int)
    out['boundary'] = np.ones(n, dtype=int)
    wall_area = np.zeros(n)
    thick_sum = np.zeros(n)
    for r in range(len(ridges)):
        a, b = cells[r]
        out['nwalls'][a] += 1
        thick_sum[a] += thickness[r]
        if b >= 0:
            out['nwalls'][b] += 1
            out['nneighbours'][a] += 1
            out['nneighbours'][b] += 1
            wall_area[a] += length[r]*thickness[r]/2
            wall_area[b] += length[r]*thickness[r]/2
            thick_sum[b] += thickness[r]
        else:
            wall_area[a] += length[r]*thickness[r]
    center = np.asarray(result.config.generation_center, dtype=float)
    vertex_boundary = np.asarray(result.vertex_boundary, dtype=bool)
    has_walls = out['nwalls'] > 0
    out['wall_area'][has_walls] = wall_area[has_walls]  # also for the cells without a closed polygon
    out['mean_thickness'][has_walls] = thick_sum[has_walls]/out['nwalls'][has_walls]
    for i, poly in enumerate(polygons):
        if len(poly) < 3:
            continue
        p = vertices[poly]
        area, centroid, cov = polygon_moments(p)
        if area <= 0:
            continue
        out['boundary'][i] = int(vertex_boundary[poly].any())
        out['area'][i] = area
        out['perimeter'][i] = np.sum(np.linalg.norm(np.roll(p, -1, axis=0) - p, axis=1))
        out['solid_fraction'][i] = min(wall_area[i]/area, 1.0)
        out['equivalent_diameter'][i] = math.sqrt(4*area/math.pi)
        values, vectors = np.linalg.eigh(cov)
        lo, hi = max(values[0], 1e-300), max(values[1], 1e-300)
        out['aspect_ratio'][i] = math.sqrt(hi/lo)
        major = vectors[:, 1]
        orientation = math.atan2(major[1], major[0])
        if orientation <= -math.pi/2:
            orientation += math.pi
        elif orientation > math.pi/2:
            orientation -= math.pi
        out['orientation'][i] = orientation
        radial = math.atan2(sites[i, 1] - center[1], sites[i, 0] - center[0])
        out['radial_elongation'][i] = abs(math.cos(orientation - radial))
    return out


def wall_statistics(result, boundaries=None):
    """One value per wall (ridge) of a Result (dict of arrays, keys WALL_FIELDS): its vertices, the cells on
    its two sides (-1 along the box), length, thickness, angle, growth ring, dense, vessel, ray, boundary"""
    vertices = np.asarray(result.vertices, dtype=float)
    ridges = np.asarray(result.ridges, dtype=int)
    cells = ridge_cells(result.sites, vertices, ridges, boundaries)
    vector = vertices[ridges[:, 1]] - vertices[ridges[:, 0]]
    return {'ridge': np.arange(len(ridges)), 'vertex1': ridges[:, 0], 'vertex2': ridges[:, 1],
            'cell1': cells[:, 0], 'cell2': cells[:, 1],
            'length': np.linalg.norm(vector, axis=1),
            'thickness': np.asarray(result.vertex_thickness, dtype=float)[ridges[:, 1]],
            'angle': np.arctan2(vector[:, 1], vector[:, 0]),
            'ring': np.asarray(result.ridge_ring, dtype=int), 'dense': np.asarray(result.ridge_dense, dtype=int),
            'vessel': np.asarray(result.ridge_vessel, dtype=int), 'ray': np.asarray(result.ridge_ray, dtype=int),
            'boundary': np.asarray(result.ridge_boundary, dtype=int)}


def ring_statistics(cells):
    """
    One row per band of the growth rings that has cells (dict of arrays, keys RING_FIELDS), from the
    cell statistics: the band, growth ring and dense flag, the number of cells (all, and the interior
    ones not cut by the box), the area of the cells and of their walls, the 2D porosity
    (1 - wall area / area), and the mean area, equivalent diameter, aspect ratio, number of neighbours
    and wall thickness of the interior cells.
    """
    bands = np.unique(cells['band'])
    out = {name: [] for name in RING_FIELDS}
    for b in bands:
        sel = (cells['band'] == b) & np.isfinite(cells['area'])
        interior = sel & (cells['boundary'] == 0)
        ring, dense = rpgen.RingLabels(np.array([b]))
        out['band'].append(int(b))
        out['ring'].append(int(ring[0]))
        out['dense'].append(int(dense[0]))
        out['ncells'].append(int(np.sum(cells['band'] == b)))
        out['ncells_interior'].append(int(interior.sum()))
        area = float(np.sum(cells['area'][sel]))
        wall = float(np.sum(cells['wall_area'][sel]))
        out['area'].append(area)
        out['wall_area'].append(wall)
        out['porosity_2d'].append(1 - min(wall/area, 1.0) if area > 0 else np.nan)
        for name, key in [('mean_cell_area', 'area'), ('mean_equivalent_diameter', 'equivalent_diameter'),
                          ('mean_aspect_ratio', 'aspect_ratio'), ('mean_nneighbours', 'nneighbours'),
                          ('mean_thickness', 'mean_thickness')]:
            out[name].append(float(np.mean(cells[key][interior])) if interior.any() else np.nan)
    return {name: np.asarray(values) for name, values in out.items()}


def summary(cells, rings, walls):
    """Descriptors of the whole mesh for the batch table (dict, keys SUMMARY_FIELDS): the number of cells, the
    mean area, equivalent diameter, aspect ratio and number of neighbours of the interior cells, the 2D porosity
    of all cells and of the sparse and the dense rings, and the mean wall length and thickness"""
    interior = (cells['boundary'] == 0) & np.isfinite(cells['area'])
    finite = np.isfinite(cells['area'])
    mean = lambda key, sel: float(np.mean(cells[key][sel])) if sel.any() else np.nan
    area, wall = float(np.sum(cells['area'][finite])), float(np.sum(cells['wall_area'][finite]))
    def porosity(sel):
        a, w = float(np.sum(cells['area'][sel])), float(np.sum(cells['wall_area'][sel]))
        return 1 - min(w/a, 1.0) if a > 0 else np.nan
    return {'ncells': int(finite.sum()), 'mean_cell_area': mean('area', interior),
            'mean_equivalent_diameter': mean('equivalent_diameter', interior),
            'mean_aspect_ratio': mean('aspect_ratio', interior), 'mean_nneighbours': mean('nneighbours', interior),
            'porosity_2d': 1 - min(wall/area, 1.0) if area > 0 else np.nan,
            'porosity_2d_sparse': porosity(finite & (cells['dense'] == 0)),
            'porosity_2d_dense': porosity(finite & (cells['dense'] == 1)),
            'mean_wall_length': float(np.mean(walls['length'])) if len(walls['length']) else np.nan,
            'mean_wall_thickness': float(np.mean(walls['thickness'])) if len(walls['thickness']) else np.nan}


def write_table(file, table, fields):
    """Write a dict of equal-length arrays as a CSV file with the given columns"""
    with open(file, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(fields)
        columns = [np.asarray(table[name]) for name in fields]
        for row in zip(*columns):
            writer.writerow([_text(v) for v in row])


def _text(value):
    if isinstance(value, (np.integer, int)):
        return int(value)
    if isinstance(value, (np.floating, float)):
        return '' if np.isnan(value) else repr(float(value))
    return value


def write_statistics(geoName, cells, rings, walls, path='meshes'):
    """Write geoName-cells.csv, geoName-rings.csv and geoName-walls.csv (see cell_statistics,
    ring_statistics and wall_statistics)"""
    folder = Path(path, geoName)
    write_table(folder / (geoName + '-cells.csv'), cells, CELL_FIELDS)
    write_table(folder / (geoName + '-rings.csv'), rings, RING_FIELDS)
    write_table(folder / (geoName + '-walls.csv'), walls, WALL_FIELDS)
