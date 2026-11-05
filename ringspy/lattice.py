# -*- coding: utf-8 -*-
"""
Lattices with fibers of various lengths

In the regular extrusion, every grain (beam along a cell corner) is cut into the same segments,
so the longitudinal connectors (the joints between the fibers of a grain, where fibers end) lie
on common planes, like knife cuts through the specimen. Here every grain is its own sequence of
fibers of random length, ending at different heights, as the cells of a fibrous material do.

The specimen height is divided into a grid of `nelem` quadratic Bezier beam elements of length h.
A grain's fiber ends lie on interior grid lines; the fiber above and the fiber below a fiber end
have separate nodes there (connector_l_length apart, symmetric about the grid line), joined by a
longitudinal connector per wing of the grain. On every grid line, every cell wall (ridge) gets a
transverse connector between the nodes of its two grains: one of height h where both grains are
continuous (a "reg" connector), or two of height h/2, one below the grid line between the nodes
of the fibers ending there or continuing ("top"), one above between the nodes of the fibers
starting there or continuing ("bot"), as at the bottom and top faces of the specimen. With the
same fiber ends on every grain the lattice is the regular one, with 2 elements per segment.

author: king_yin3613
email: haoyin2022@u.northwestern.edu
"""

import math
from dataclasses import dataclass
from pathlib import Path
from typing import List

import numpy as np
from scipy.spatial import cKDTree

from . import MeshGenTools as rpgen

NURBS_DEGREE = 2


def fiber_end_layers(ngrain, nelem, h, fiber_length, fiber_length_cov=0.3, rng=None):
    """
    Grid lines (1 .. nelem-1) on which the fibers of each grain end.

    The fiber lengths are drawn from a lognormal distribution with mean fiber_length (mm) and
    coefficient of variation fiber_length_cov, independently for every grain; the first fiber of a
    grain crosses the bottom face of the specimen at a random point of its length, so the ends are
    not aligned. Each end is rounded to the nearest of the nelem grid lines of spacing h.
    rng :: None (global NumPy random state), an int seed, or a np.random.Generator.
    Returns a list of ngrain sorted integer arrays.
    """
    rng = rpgen._rng(rng)
    if fiber_length <= 0:
        raise ValueError("fiber_length must be positive")
    if fiber_length_cov < 0:
        raise ValueError("fiber_length_cov must not be negative")
    sigma = math.sqrt(math.log(1 + fiber_length_cov**2))
    mu = math.log(fiber_length) - sigma**2/2
    height = nelem*h
    ends = []
    for _ in range(ngrain):
        length = rng.lognormal(mu, sigma)
        z = length*(1 - rng.random())  # end of the fiber crossing the bottom face
        ks = []
        while z < height:
            k = int(round(z/h))
            if 0 < k < nelem:
                ks.append(k)
            z += rng.lognormal(mu, sigma)
        ends.append(np.unique(np.array(ks, dtype=int)))
    return ends


def cell_end_heights(cells, z_min, z_max, fiber_length, fiber_length_cov=0.3, rng=None):
    """
    Heights (mm) inside the specimen at which each cell ends and the next cell begins: cell lengths from
    a lognormal distribution of mean fiber_length and coefficient of variation fiber_length_cov, the first
    cell of each crossing the bottom face at a random point of its length. `cells` is the number of
    cells or a list of their ids; returns {cell id: array of heights}.
    rng :: None (global NumPy random state), an int seed, or a np.random.Generator.
    """
    rng = rpgen._rng(rng)
    if fiber_length <= 0:
        raise ValueError("fiber_length must be positive")
    if fiber_length_cov < 0:
        raise ValueError("fiber_length_cov must not be negative")
    sigma = math.sqrt(math.log(1 + fiber_length_cov**2))
    mu = math.log(fiber_length) - sigma**2/2
    ids = range(cells) if isinstance(cells, (int, np.integer)) else list(cells)
    ends = {}
    for cell in ids:
        z = z_min + rng.lognormal(mu, sigma)*(1 - rng.random())
        heights = []
        while z < z_max:
            heights.append(z)
            z += rng.lognormal(mu, sigma)
        ends[cell] = np.array(heights)
    return ends


def grain_ends_from_cells(vertex_cells, cell_ends, z_min, h, nelem):
    """
    Grid lines (1 .. nelem-1, spacing h from z_min) where the fibers of each grain end: wherever one of
    the cells meeting at the grain ends (vertex_cells: the cell ids of each vertex; cell_ends:
    {cell id: heights}, e.g. from cell_end_heights). Returns a list of sorted integer arrays.
    """
    ends = []
    for cells in vertex_cells:
        ks = set()
        for cell in cells:
            for z in cell_ends.get(cell, ()):
                k = int(round((z - z_min)/h))
                if 0 < k < nelem:
                    ks.add(k)
        ends.append(np.array(sorted(ks), dtype=int))
    return ends


def cells_of_vertices(vertices, ridges, sites, boundaries=None, tol=None):
    """
    The cells (sites of the Voronoi tessellation, by index) meeting at each vertex of a clipped Voronoi
    mesh: the two sites nearest to the midpoint of a ridge are the cells the ridge separates; a ridge
    along a side of the clipping box (boundaries, as Clipping_Box gives them) belongs to one cell.
    Returns a list of sets of site indices.
    """
    vertices = np.asarray(vertices, dtype=float)
    ridges = np.asarray(ridges)[:, 0:2].astype(int)
    sites = np.asarray(sites, dtype=float)[:, 0:2]
    mids = (vertices[ridges[:, 0]] + vertices[ridges[:, 1]])/2
    k = min(2, len(sites))
    distance, index = cKDTree(sites).query(mids, k=k)
    index = np.asarray(index).reshape(len(mids), k)
    along_side = np.zeros(len(mids), dtype=bool)
    if boundaries is not None:
        starts = np.array([b[1][0] for b in boundaries], dtype=float)
        ends = np.array([b[1][1] for b in boundaries], dtype=float)
        segments = ends - starts
        lengths = np.linalg.norm(segments, axis=1)
        if tol is None:
            tol = 1e-6*lengths.max()
        for p0, p1, seg, length in zip(starts, ends, segments, lengths):
            rel = mids - p0
            u = np.clip((rel @ seg)/length**2, 0, 1)
            along_side |= np.linalg.norm(rel - u[:, np.newaxis]*seg, axis=1) <= tol
    cells = [set() for _ in range(len(vertices))]
    for r, (a, b) in enumerate(ridges):
        members = index[r, 0:1] if along_side[r] or k < 2 else index[r]
        for cell in members.tolist():
            cells[a].add(cell)
            cells[b].add(cell)
    return cells


def taper_vertices(vertex_layers, z_layers, vertex_cells, cell_ends, taper_length, taper_ratio=0.15, kind=None,
                   boundaries=None):
    """
    Taper the cells at their ends: within taper_length of an end of a cell (cell_ends: {cell id: heights}),
    the vertices of the cell (vertex_cells: the cell ids of each vertex) move toward the center of the cell
    (the mean of its vertices in that layer), so that its cross-section shrinks, with a cosine profile,
    to taper_ratio times its size at the end, and the neighbouring cells fill the space; the cell
    beginning above the end widens the same way. The displacements of the cells ending near each other
    add up. Vertices of kind 2 (corners of the box) stay, vertices of kind 1 (boundary points) stay on
    their side of the box (boundaries, as Clipping_Box gives them).

    vertex_layers :: (nlayers, nvertex, 2) in-plane coordinates of the vertices on every layer at the
        heights z_layers. Returns the tapered copy.
    """
    V = np.array(vertex_layers, dtype=float, copy=True)
    z = np.asarray(z_layers, dtype=float)
    if V.ndim != 3 or V.shape[0] != len(z):
        raise ValueError("vertex_layers must have the shape (len(z_layers), nvertex, 2)")
    if not 0 <= taper_ratio <= 1:
        raise ValueError("taper_ratio must be between 0 (the cell closes) and 1 (no taper)")
    if taper_length <= 0 or taper_ratio >= 1:
        return V
    members = {}
    for v, cells in enumerate(vertex_cells):
        for cell in cells:
            if len(cell_ends.get(cell, ())):
                members.setdefault(cell, []).append(v)
    displacement = np.zeros_like(V)
    for cell, verts in members.items():
        verts = np.array(sorted(set(verts)), dtype=int)
        if len(verts) < 2:
            continue
        for z_end in cell_ends[cell]:
            t = (z - z_end)/taper_length
            layers = np.nonzero(np.abs(t) < 1)[0]
            if len(layers) == 0:
                continue
            w = 0.5*(1 + np.cos(math.pi*t[layers]))*(1 - taper_ratio)
            p = V[layers][:, verts]
            center = p.mean(axis=1, keepdims=True)
            displacement[np.ix_(layers, verts)] -= w[:, np.newaxis, np.newaxis]*(p - center)
    tapered = V + displacement
    if kind is not None:
        kind = np.asarray(kind)
        tapered[:, kind == 2] = V[:, kind == 2]
        if boundaries is not None and np.any(kind == 1):
            rows = np.nonzero(kind == 1)[0]
            for m in range(len(z)):
                tapered[m, rows] = snap_to_boundary(tapered[m, rows], V[m, rows], boundaries)
    return tapered


def snap_to_boundary(points, reference, boundaries):
    """Project points onto the boundary line (segment of `boundaries`, as Clipping_Box gives them) closest
    to their reference positions"""
    points = np.asarray(points, dtype=float).reshape(-1, 2)
    reference = np.asarray(reference, dtype=float).reshape(-1, 2)
    starts = np.array([b[1][0] for b in boundaries], dtype=float)
    ends = np.array([b[1][1] for b in boundaries], dtype=float)
    segments = ends - starts
    lengths2 = np.einsum('mk,mk->m', segments, segments)
    rel = reference[:, np.newaxis, :] - starts[np.newaxis, :, :]
    u = np.clip(np.einsum('nmk,mk->nm', rel, segments)/lengths2, 0, 1)
    nearest = np.linalg.norm(rel - u[:, :, np.newaxis]*segments[np.newaxis, :, :], axis=2).argmin(axis=1)
    seg = segments[nearest]
    u = np.clip(np.einsum('nk,nk->n', points - starts[nearest], seg)/lengths2[nearest], 0, 1)
    return starts[nearest] + u[:, np.newaxis]*seg


def fiber_statistics(fiber_ends, nelem, h):
    """Numbers and lengths of the fibers of a lattice (dict, for the report log)"""
    lengths = []
    nends = 0
    for ks in fiber_ends:
        bounds = np.concatenate(([0], np.asarray(ks, dtype=int), [nelem]))
        lengths.extend((np.diff(bounds)*h).tolist())
        nends += len(ks)
    lengths = np.asarray(lengths, dtype=float)
    return {'nfibers': int(len(lengths)), 'nfiber_ends': int(nends),
            'mean fiber length (cut by the faces)': float(lengths.mean()) if len(lengths) else 0.0,
            'longest fiber': float(lengths.max()) if len(lengths) else 0.0,
            'grains without a fiber end': int(sum(1 for ks in fiber_ends if len(ks) == 0))}


@dataclass
class Lattice:
    """The nodes, beams and connectors of a lattice with fibers of various lengths"""
    nodes: np.ndarray                # (nnode, 3) coordinates of the control points (nodes)
    node_grain: np.ndarray           # grain (2D vertex) of each node
    node_layer: np.ndarray           # half-layer m (0 .. 2*nelem) of each node, z = z_min + m*h/2 (+- connector_l_length/2 at a fiber end)
    patches: List[np.ndarray]        # node indices (0-based) of each fiber: an odd number, 2 elements... per fiber
    patch_grain: np.ndarray
    element_nodes: np.ndarray        # (ngrain, nelem, 3) node indices (0-based) of the element of each grain on each grid interval
    beam_connectivity: np.ndarray    # (nbeamElem, 3) 1-based, one row per Bezier element
    connector_t_bot: np.ndarray      # (n, 2) 1-based node pairs
    connector_t_reg: np.ndarray
    connector_t_top: np.ndarray
    connector_t_ridge: tuple         # ridge of each transverse connector (bot, reg, top arrays)
    connector_t_line: tuple          # grid line k of each transverse connector (bot, reg, top arrays)
    connector_l: np.ndarray          # (n, 2) 1-based node pairs
    connector_l_wing: np.ndarray     # wing (row of flattened_all_vertices_2D) of each longitudinal connector
    Meshdata: np.ndarray             # (nconnector, 13): [inode jnode centerx centery centerz dx1 dy1 dz1 dx2 dy2 dz2 width height]
    conn_t_tangents: np.ndarray      # (nconnector_t, 3) unit tangents along the grains
    conn_l_tangents: np.ndarray      # (nconnector_l, 3) in-plane directions of the wings
    fiber_ends: list                 # grid lines of the fiber ends of each grain
    z_grid: np.ndarray               # (nelem+1,) heights of the grid lines
    h: float                         # length of the beam elements
    connector_l_length: float

    @property
    def nconnector_t(self):
        return len(self.connector_t_bot) + len(self.connector_t_reg) + len(self.connector_t_top)

    @property
    def nconnector_l(self):
        return len(self.connector_l)

    @property
    def nbeamElem(self):
        return len(self.beam_connectivity)

    @property
    def npatch(self):
        return len(self.patches)


def build_lattice(all_ridges, nvertex, nridge, all_vertices_2D, max_wings, flattened_all_vertices_2D, layer_points,
                  z_min, z_max, nelem, long_connector_ratio, fiber_ends, ridge_widths=None):
    """
    Build the lattice of a 2D mesh with fibers ending on the given grid lines.

    all_ridges, nvertex, nridge, all_vertices_2D, max_wings, flattened_all_vertices_2D :: the 2D mesh
        (RidgeMidQuarterPts, VertexandRidgeinfo).
    layer_points :: in-plane coordinates of all points of the 2D mesh (vertices, ridge mid and quarter
        points) on every half-layer m = 0 .. 2*nelem, shape (2*nelem+1, npt_per_layer_vtk, 2): the
        2D mesh rotated by the twist, or the cell corners followed through a 3D image.
    z_min, z_max, nelem :: the grid: nelem beam elements of length h = (z_max - z_min)/nelem.
    long_connector_ratio :: length of the longitudinal connectors = ratio * 2h (2h: the segment length
        of the regular extrusion).
    fiber_ends :: for each grain, the grid lines (1 .. nelem-1) where its fibers end (fiber_end_layers).
    ridge_widths :: wall thickness of each ridge, shape (nridge,), or on each grid line, shape
        (nelem+1, nridge); default: the widths of the 2D mesh (all_vertices_2D).
    """
    layer_points = np.asarray(layer_points, dtype=float)
    all_ridges = np.asarray(all_ridges)
    nhalf = 2*nelem + 1
    if layer_points.shape[0] != nhalf or layer_points.shape[1] < nvertex:
        raise ValueError("layer_points must have the shape (2*nelem+1, npoints >= nvertex, 2), not {}".format(layer_points.shape))
    if nelem < 1:
        raise ValueError("nelem must be at least 1")
    h = (z_max - z_min)/nelem
    z_half = z_min + np.arange(nhalf)*h/2
    z_grid = z_half[0::2]
    delta = 2*h*long_connector_ratio
    ridges = all_ridges[:, 0:2].astype(int)
    if ridge_widths is None:
        ridge_widths = all_vertices_2D[ridges[:, 1], 3 + max_wings*3]
    ridge_widths = np.asarray(ridge_widths, dtype=float)
    if ridge_widths.ndim == 1:
        ridge_widths = np.broadcast_to(ridge_widths.reshape(1, -1), (nelem + 1, nridge))
    if ridge_widths.shape != (nelem + 1, nridge):
        raise ValueError("ridge_widths must have the shape (nridge,) or (nelem+1, nridge), not {}".format(ridge_widths.shape))

    is_end = np.zeros((nvertex, nelem + 1), dtype=bool)
    for g, ks in enumerate(fiber_ends):
        ks = np.asarray(ks, dtype=int)
        if np.any((ks < 1) | (ks > nelem - 1)):
            raise ValueError("fiber ends must lie on the interior grid lines 1 .. {:d}".format(nelem - 1))
        is_end[g, ks] = True

    # nodes, fiber by fiber: a fiber end has a node below and a node above the grid line
    coords = []
    node_grain = []
    node_layer = []
    lower = np.full((nvertex, nhalf), -1, dtype=int)  # node of grain g on half-layer m seen from below
    upper = np.full((nvertex, nhalf), -1, dtype=int)  # ... seen from above
    patches = []
    patch_grain = []
    for g in range(nvertex):
        bounds = np.concatenate(([0], np.nonzero(is_end[g])[0], [nelem]))
        for k0, k1 in zip(bounds[:-1], bounds[1:]):
            patch = []
            for m in range(2*k0, 2*k1 + 1):
                z = z_half[m]
                if m == 2*k0 and k0 > 0:
                    z += delta/2
                if m == 2*k1 and k1 < nelem:
                    z -= delta/2
                idx = len(coords)
                coords.append((layer_points[m, g, 0], layer_points[m, g, 1], z))
                node_grain.append(g)
                node_layer.append(m)
                patch.append(idx)
                if m == 2*k0:
                    upper[g, m] = idx
                elif m == 2*k1:
                    lower[g, m] = idx
                else:
                    lower[g, m] = upper[g, m] = idx
            patches.append(np.array(patch, dtype=int))
            patch_grain.append(g)
    nodes = np.array(coords, dtype=float)
    node_grain = np.array(node_grain, dtype=int)
    node_layer = np.array(node_layer, dtype=int)

    # beam elements (3 control points each), and the element of each grain on each grid interval
    beam = np.vstack([np.column_stack((patch[0:-2:2], patch[1:-1:2], patch[2::2])) for patch in patches])
    element_nodes = np.stack((upper[:, 0:-2:2], upper[:, 1:-1:2], lower[:, 2::2]), axis=-1)  # (ngrain, nelem, 3)

    # transverse connectors on every grid line
    bot, reg, top = [], [], []
    bot_r, reg_r, top_r = [], [], []
    bot_k, reg_k, top_k = [], [], []
    a, b = ridges[:, 0], ridges[:, 1]
    for k in range(nelem + 1):
        m = 2*k
        if k == 0:
            rows = np.arange(nridge)
            bot.append(np.column_stack((upper[a, m], upper[b, m])))
            bot_r.append(rows)
            bot_k.append(np.full(nridge, k))
        elif k == nelem:
            rows = np.arange(nridge)
            top.append(np.column_stack((lower[a, m], lower[b, m])))
            top_r.append(rows)
            top_k.append(np.full(nridge, k))
        else:
            split = is_end[a, k] | is_end[b, k]
            rows = np.nonzero(split)[0]
            if len(rows):
                top.append(np.column_stack((lower[a[rows], m], lower[b[rows], m])))
                top_r.append(rows)
                top_k.append(np.full(len(rows), k))
                bot.append(np.column_stack((upper[a[rows], m], upper[b[rows], m])))
                bot_r.append(rows)
                bot_k.append(np.full(len(rows), k))
            rows = np.nonzero(~split)[0]
            if len(rows):
                reg.append(np.column_stack((lower[a[rows], m], lower[b[rows], m])))
                reg_r.append(rows)
                reg_k.append(np.full(len(rows), k))
    stack = lambda parts, n: np.vstack(parts) if parts else np.zeros((0, n), dtype=int)  # noqa: E731
    cat = lambda parts: np.concatenate(parts) if parts else np.zeros(0, dtype=int)  # noqa: E731
    connector_t_bot, connector_t_reg, connector_t_top = stack(bot, 2), stack(reg, 2), stack(top, 2)
    ridge_of = (cat(bot_r), cat(reg_r), cat(top_r))
    line_of = (cat(bot_k), cat(reg_k), cat(top_k))

    # longitudinal connectors: one per wing of the grain at every fiber end
    nwings = all_vertices_2D[:, 2].astype(int)
    wing_first = np.concatenate(([0], np.cumsum(nwings)))
    wing_other = flattened_all_vertices_2D[1, :].astype(int)
    wing_ridge = flattened_all_vertices_2D[0, :].astype(int)
    conn_l, conn_l_wing, conn_l_line = [], [], []
    for g in range(nvertex):
        for k in np.nonzero(is_end[g])[0]:
            wings = np.arange(wing_first[g], wing_first[g + 1])
            conn_l.append(np.column_stack((np.full(len(wings), lower[g, 2*k]), np.full(len(wings), upper[g, 2*k]))))
            conn_l_wing.append(wings)
            conn_l_line.append(np.full(len(wings), k))
    connector_l = stack(conn_l, 2)
    connector_l_wing = cat(conn_l_wing)
    connector_l_line = cat(conn_l_line)

    # directions of the grains at their nodes, and the tangents of the transverse connectors
    node_dir = np.zeros((len(nodes), 3))
    for patch in patches:
        before = patch[np.maximum(np.arange(len(patch)) - 1, 0)]
        after = patch[np.minimum(np.arange(len(patch)) + 1, len(patch) - 1)]
        node_dir[patch] = nodes[after] - nodes[before]
    node_dir /= np.maximum(np.linalg.norm(node_dir, axis=1, keepdims=True), 1e-300)
    connector_t = np.vstack((connector_t_bot, connector_t_reg, connector_t_top))
    conn_t_tangents = node_dir[connector_t[:, 0]] + node_dir[connector_t[:, 1]]
    conn_t_tangents /= np.maximum(np.linalg.norm(conn_t_tangents, axis=1, keepdims=True), 1e-300)

    # geometry of the longitudinal connectors on their half-layer: wing directions and half lengths
    wing_vertex = np.repeat(np.arange(nvertex), nwings)
    m_l = 2*connector_l_line
    wing_vec = layer_points[m_l, wing_other[connector_l_wing]] - layer_points[m_l, wing_vertex[connector_l_wing]]
    wing_length = np.linalg.norm(wing_vec, axis=1)
    conn_l_tangents = np.zeros((len(connector_l), 3))
    conn_l_tangents[:, 0:2] = wing_vec/np.maximum(wing_length, 1e-300)[:, np.newaxis]

    # connector data, as ConnectorMeshFile: [inode jnode centerx centery centerz dx1 dy1 dz1 dx2 dy2 dz2 width height]
    nt = len(connector_t)
    connectivity = np.vstack((connector_t, connector_l))
    n = len(connectivity)
    Meshdata = np.zeros((n, 13))
    Meshdata[:, 0:2] = connectivity + 1
    p1 = nodes[connectivity[:, 0]]
    p2 = nodes[connectivity[:, 1]]
    Meshdata[:, 2:5] = (p1 + p2)/2
    Meshdata[:, 5:8] = Meshdata[:, 2:5] - p1
    Meshdata[:, 8:11] = Meshdata[:, 2:5] - p2
    heights_t = np.concatenate((np.full(len(connector_t_bot), h/2), np.full(len(connector_t_reg), h),
                                np.full(len(connector_t_top), h/2)))
    Meshdata[:nt, 12] = heights_t
    Meshdata[nt:, 12] = wing_length/2
    Meshdata[:nt, 11] = ridge_widths[np.concatenate(line_of), np.concatenate(ridge_of)]
    Meshdata[nt:, 11] = ridge_widths[connector_l_line, wing_ridge[connector_l_wing]]
    Meshdata[nt:, 2:5] += conn_l_tangents*Meshdata[nt:, 12][:, np.newaxis]/2  # eccentricity along the wing

    return Lattice(nodes=nodes, node_grain=node_grain, node_layer=node_layer, patches=patches,
                   patch_grain=np.array(patch_grain, dtype=int), element_nodes=element_nodes,
                   beam_connectivity=beam + 1, connector_t_bot=connector_t_bot + 1, connector_t_reg=connector_t_reg + 1,
                   connector_t_top=connector_t_top + 1, connector_t_ridge=ridge_of, connector_t_line=line_of,
                   connector_l=connector_l + 1, connector_l_wing=connector_l_wing, Meshdata=Meshdata,
                   conn_t_tangents=conn_t_tangents, conn_l_tangents=conn_l_tangents,
                   fiber_ends=[np.asarray(ks, dtype=int) for ks in fiber_ends], z_grid=z_grid, h=h,
                   connector_l_length=delta)


def write_connector_file(geoName, lattice, path='meshes'):
    """The connector data file geoName-mesh.txt, as ConnectorMeshFile writes it"""
    rpgen._savetxt(Path(path, geoName, geoName + '-mesh.txt'), lattice.Meshdata, fmt='%.16g', delimiter=' ',
                   header='# Connector Data Generated with RingsPy Mesh Generation Tool\n'
                   'Number of bot connectors\n' + str(len(lattice.connector_t_bot)) + '\n'
                   'Number of reg connectors\n' + str(len(lattice.connector_t_reg)) + '\n'
                   'Number of top connectors\n' + str(len(lattice.connector_t_top)) + '\n'
                   'Number of long connectors\n' + str(len(lattice.connector_l)) + '\n'
                   '[inode jnode centerx centery centerz dx1 dy1 dz1 dx2 dy2 dz2 width height]', comments='')


def write_iga_file(geoName, lattice, path='meshes'):
    """
    The NURBS beam file geoNameIGA.txt, as BezierBeamFile writes it. The patches (fibers) have
    various numbers of elements: the header gives the largest numbers of control points and elements
    of a patch, and each patch's knot vector tells its own: as in BezierBeamFile, a patch of n
    elements has the knot vector [0, 0, 0, 1/n, ..., 1, 1, 1] of n + 5 knots and 2n + 1 control points.
    """
    sizes = [len(patch) for patch in lattice.patches]
    nelem_patch = [(size - 1)//NURBS_DEGREE for size in sizes]
    with open(Path(path, geoName, geoName + 'IGA.txt'), 'w') as txtfile:
        txtfile.write('# Dimension of beam elements \n')
        txtfile.write('1 \n')
        txtfile.write('# Order of basis function \n')
        txtfile.write('{:d} \n'.format(NURBS_DEGREE))
        txtfile.write('# Number of control points per patch \n')
        txtfile.write('{:d} \n'.format(max(sizes)))
        txtfile.write('# Number of elements per patch \n')
        txtfile.write('{:d} \n'.format(max(nelem_patch)))
        txtfile.write('# Number of Patches \n')
        txtfile.write('{:d} \n'.format(len(sizes)))
        for ipatch, ne in enumerate(nelem_patch):
            knots = np.concatenate((np.zeros(NURBS_DEGREE), np.linspace(0, 1, ne + 1), np.ones(NURBS_DEGREE)))
            txtfile.write('{:s} \n'.format('PATCH-' + str(ipatch + 1)))
            txtfile.write('Size of knot vectors \n')
            txtfile.write('{:d} \n'.format(len(knots)))
            txtfile.write('knot vectors \n')
            txtfile.write(''.join('{:f} '.format(knot) for knot in knots))
            txtfile.write('\n')


def _connector_cross_sections(Meshdata, tangents, nconnector_t, nbot, ntop):
    """Corners of the cross-sections (4 per connector) and of the ends (8 per connector) of the connectors,
    as VisualizationFiles draws them: the bot/top connectors are shifted off their grid line by half their height"""
    data = Meshdata.copy()
    n = len(data)
    data[nconnector_t:, 5:7] = 0  # longitudinal connectors are drawn vertical
    data[nconnector_t:, 8:10] = 0
    center = data[:, 2:5].copy()
    normal1 = data[:, 5:8]
    normal2 = data[:, 8:11]
    length1 = np.linalg.norm(normal1, axis=1)
    length2 = np.linalg.norm(normal2, axis=1)
    degenerate = length1 == 0
    normal = normal1/np.where(degenerate, 1.0, length1)[:, np.newaxis]
    normal[degenerate] = [1.0, 0.0, 0.0]
    width = data[:, 11]
    height = data[:, 12]
    tangent = tangents/np.maximum(np.linalg.norm(tangents, axis=1, keepdims=True), 1e-300)
    bitangent = np.cross(normal, tangent)
    center[0:nbot] += tangent[0:nbot]*height[0:nbot, np.newaxis]/2
    center[nconnector_t - ntop:nconnector_t] -= tangent[nconnector_t - ntop:nconnector_t]*height[nconnector_t - ntop:nconnector_t, np.newaxis]/2
    t = tangent*height[:, np.newaxis]/2
    bt = bitangent*width[:, np.newaxis]/2
    corners = [center - t - bt, center - t + bt, center + t + bt, center + t - bt]
    quad = np.vstack(corners)
    hexa = np.vstack([c - normal*length1[:, np.newaxis] for c in corners] + [c + normal*length2[:, np.newaxis] for c in corners])
    hex1 = np.arange(8*n).reshape(8, -1).T
    hex2 = hex1.copy()
    hex2[:, 4:8] += 4*n
    return quad, hexa, np.vstack((hex1, hex2))


def element_bands(lattice, vertex_band, ridge_band, wing_ridge):
    """
    Ring bands (see MeshGenTools.RingBands) of the elements of a lattice, for the ring element sets of
    AbaqusFile: {'Beams': band of each beam element (that of its grain), 'TConns': of each transverse
    connector (bot, reg, top in a row: that of its ridge), 'LConns': of each longitudinal connector
    (that of the ridge of its wing)}. vertex_band, ridge_band: the bands of the vertices and ridges of the
    2D mesh; wing_ridge: the ridge of each wing (flattened_all_vertices_2D[0, :]).
    """
    vertex_band = np.asarray(vertex_band).reshape(-1)
    ridge_band = np.asarray(ridge_band).reshape(-1)
    wing_ridge = np.asarray(wing_ridge, dtype=int).reshape(-1)
    beams = lattice.beam_connectivity - 1
    return {'Beams': vertex_band[lattice.node_grain[beams[:, 0]]],
            'TConns': ridge_band[np.concatenate(lattice.connector_t_ridge)],
            'LConns': ridge_band[wing_ridge[lattice.connector_l_wing]]}


def write_vtu_files(geoName, lattice, layer_points, all_ridges, nvertex, path='meshes', vtk_format='ascii',
                    vertex_band=None, ridge_band=None, extra_labels=None):
    """
    The Paraview files geoName_vertices.vtu, geoName_beams.vtu, geoName_conns.vtu and
    geoName_conns_vol.vtu of a lattice, as VisualizationFiles writes them: the nodes (with the IGA
    collocation flag), the beams with the wings of the cell walls, the connectors as lines and
    cross-sections, and the connector volumes. With vertex_band and ridge_band (the ring bands of
    the vertices and ridges of the 2D mesh, see MeshGenTools.RingBands), every cell gets the cell
    data 'Ring' and 'Dense' as VisualizationFiles writes them; extra_labels ({name: (values per vertex,
    values per ridge)}) adds more integer cell data in the same way. Returns the connector volumes
    ({'points', 'hexahedra'}) for the STL file, as VisualizationFiles does.
    """
    if vtk_format not in ['ascii', 'binary']:
        raise ValueError("vtk_format must be 'ascii' or 'binary', not {!r}".format(vtk_format))
    layer_points = np.asarray(layer_points, dtype=float)
    all_ridges = np.asarray(all_ridges)
    nhalf, npt_vtk = layer_points.shape[0:2]
    nextra = npt_vtk - nvertex
    nnode = len(lattice.nodes)
    # points: the nodes, then the ridge mid and quarter points of every half-layer
    z_half = lattice.z_grid[0] + np.arange(nhalf)*lattice.h/2
    extra = np.concatenate((layer_points[:, nvertex:, :], np.broadcast_to(z_half[:, np.newaxis, np.newaxis], (nhalf, nextra, 1))), axis=2)
    points = np.vstack((lattice.nodes, extra.reshape(-1, 3)))
    npoints = len(points)

    def grid_point(m, p):  # index of mesh point p (>= nvertex) on half-layer m
        return nnode + m*nextra + (p - nvertex)

    collocation = np.zeros(npoints)
    collocation[:nnode] = (lattice.node_layer % 2 == 0).astype(float)

    # beams: quadratic edges, and quadratic quads for the two wings of every cell wall on every element
    beam_vtk = lattice.beam_connectivity[:, [0, 2, 1]] - 1
    nelem = lattice.element_nodes.shape[1]
    quads = []
    for wing in (0, 1):
        vertex = all_ridges[:, wing]
        mid = all_ridges[:, 2 + wing]
        qrt = all_ridges[:, 4 + wing]
        elements = lattice.element_nodes[vertex]  # (nridge, nelem, 3)
        m0 = np.broadcast_to(2*np.arange(nelem)[np.newaxis, :], elements.shape[0:2])
        quads.append(np.stack((elements[:, :, 0], elements[:, :, 2],
                               grid_point(m0 + 2, mid[:, np.newaxis]), grid_point(m0, mid[:, np.newaxis]),
                               elements[:, :, 1], grid_point(m0 + 2, qrt[:, np.newaxis]),
                               grid_point(m0 + 1, mid[:, np.newaxis]), grid_point(m0, qrt[:, np.newaxis])), axis=-1).reshape(-1, 8))
    quad_connectivity = np.vstack(quads)

    # connectors
    nt, nl = lattice.nconnector_t, lattice.nconnector_l
    tangents = np.vstack((lattice.conn_t_tangents, lattice.conn_l_tangents))
    quad_pts, hex_pts, hex_conn = _connector_cross_sections(lattice.Meshdata, tangents, nt, len(lattice.connector_t_bot),
                                                            len(lattice.connector_t_top))
    conn_lines = (lattice.Meshdata[:, 0:2] - 1).astype(np.int64)
    quad_conn = np.arange(4*(nt + nl)).reshape(4, -1).T + npoints
    hex_conn = hex_conn + npoints
    widths = np.tile(lattice.Meshdata[:, 11], 2)

    # growth ring and sparse/dense labels (and extra labels) of the cells of each file, in the order of the cells
    # above: the cells of a vertex (beams, wing quads, nodes) take the value of the vertex, the connectors (and the
    # mid and quarter points) that of their ridge
    def cell_values(vertex_values, ridge_values):
        vertex_values = np.asarray(vertex_values).reshape(-1).astype(np.int32)
        ridge_values = np.asarray(ridge_values).reshape(-1).astype(np.int32)
        extra_values = np.empty(nextra, dtype=np.int32)  # the mid and quarter points take the value of their ridge
        extra_values[all_ridges[:, 2:6].reshape(-1) - nvertex] = np.repeat(ridge_values, 4)
        point = np.concatenate((vertex_values[lattice.node_grain], np.tile(extra_values, nhalf)))
        quad = np.concatenate([np.repeat(vertex_values[all_ridges[:, wing]], nelem) for wing in (0, 1)])
        order = np.argsort(np.asarray(all_ridges)[:, 0:2].reshape(-1), kind='stable')
        wing_ridge = order // 2  # ridge of each wing, wings sorted by vertex
        conn = np.concatenate((ridge_values[np.concatenate(lattice.connector_t_ridge)],
                               ridge_values[wing_ridge[lattice.connector_l_wing]]))
        return point, np.concatenate((vertex_values[lattice.node_grain[beam_vtk[:, 0]]], quad)), np.tile(conn, 2)
    label_arrays = []  # (name, values per vertex, values per ridge)
    if vertex_band is not None and ridge_band is not None:
        ring_v, dense_v = rpgen.RingLabels(np.asarray(vertex_band).reshape(-1))
        ring_r, dense_r = rpgen.RingLabels(np.asarray(ridge_band).reshape(-1))
        label_arrays += [('Ring', ring_v, ring_r), ('Dense', dense_v, dense_r)]
    for name, (vertex_values, ridge_values) in (extra_labels or {}).items():
        label_arrays.append((name, vertex_values, ridge_values))
    ring_vertices, ring_beams, ring_conns = [], [], []
    for name, vertex_values, ridge_values in label_arrays:
        values_vertices, values_beams, values_conns = cell_values(vertex_values, ridge_values)
        ring_vertices.append((name, values_vertices))
        ring_beams.append((name, values_beams))
        ring_conns.append((name, values_conns))

    files = [
        ('_vertices', [points], [np.arange(npoints).reshape(-1, 1)], np.ones(npoints, dtype=int),
         [('IGAcollocation_flag', collocation)], ring_vertices),
        ('_beams', [points], [beam_vtk, quad_connectivity],
         np.concatenate((21*np.ones(len(beam_vtk)), 23*np.ones(len(quad_connectivity)))).astype(int), [], ring_beams),
        ('_conns', [points, quad_pts], [conn_lines[:nt], conn_lines[nt:], quad_conn],
         np.concatenate((3*np.ones(nt + nl), 9*np.ones(nt + nl))).astype(int), [],
         [('Connector_width', widths)] + ring_conns),
        ('_conns_vol', [points, quad_pts, hex_pts], [hex_conn], (12*np.ones(2*(nt + nl))).astype(int), [],
         [('Connector_width', widths)] + ring_conns),
    ]
    if vtk_format == 'binary':
        for suffix, pts, cells, types, point_data, cell_data in files:
            rpgen._write_vtu_appended(Path(path, geoName, geoName + suffix + '.vtu'), np.vstack(pts), cells, types,
                                      point_data, cell_data)
    else:
        rpgen._write_vtu_files_ascii(path, geoName, files)
    # the connector volumes (hexahedra), for the STL file (MeshGenTools.StlModelFile)
    return {'points': np.vstack((points, quad_pts, hex_pts)), 'hexahedra': hex_conn}
