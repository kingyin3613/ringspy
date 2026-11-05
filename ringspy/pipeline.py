# -*- coding: utf-8 -*-
"""
The full mesh generation: cell placement, Voronoi tessellation and clipping,
longitudinal extrusion, and the output files

author: king_yin3613
email: haoyin2022@u.northwestern.edu
"""

import json
import os
import shutil
import time
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional

import matplotlib.patches as patches
import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial import Voronoi, voronoi_plot_2d

from . import MeshGenTools as rpgen
from . import lattice as rplat
from .config import Config

NURBS_DEGREE = 2


@dataclass
class Result:
    """What generate() made: the output folder and files, the cell sites, the clipped 2D cell mesh
    with its growth ring labels, and model properties"""
    config: Config               # with the defaults filled in (Config.resolved())
    output_dir: Path
    files: List[Path]
    sites: np.ndarray            # cell sites, one row (x, y, cell size) per cell
    radii: np.ndarray            # radii of the growth rings: 0, the heart, then the sparse and dense part of each ring
    nParticles: int
    nbeamElem: int
    mass: float
    bulk_volume: float
    bulk_density: float
    porosity: float
    vertices: np.ndarray         # (nvertex, 2) vertices of the clipped (and merged) 2D cell mesh: the grains, one beam each
    ridges: np.ndarray           # (nridge, 2) vertex indices of the cell walls: the transverse connectors
    vertex_ring: np.ndarray      # growth ring of each vertex: 0 the heart, 1 ... nrings (nrings+1 beyond the last ring)
    vertex_dense: np.ndarray     # True for the vertices in the heart and in the dense (latewood) part of a ring
    vertex_boundary: np.ndarray  # True for the vertices on the clipping boundary
    vertex_thickness: np.ndarray # cell wall thickness at each vertex (the thickness of its walls)
    ridge_ring: np.ndarray       # growth ring of each wall: that of its second vertex, whose thickness it takes
                                 # (the connectors of a wall, transverse and longitudinal, have this ring and thickness)
    ridge_dense: np.ndarray      # True for the walls in the heart and in the dense part of a ring
    ridge_boundary: np.ndarray   # True for the walls along the clipping boundary (both ends on it)
    site_type: np.ndarray        # type of each cell site: 0 fiber, 1 vessel, 2 ray cell (MeshGenTools.CELL_TYPES)
    vertex_vessel: np.ndarray    # True for the vertices of a vessel (one of their walls bounds a vessel)
    vertex_ray: np.ndarray       # True for the vertices of a ray cell (and not of a vessel)
    ridge_vessel: np.ndarray     # True for the walls of the vessels
    ridge_ray: np.ndarray        # True for the walls of the ray cells that are not vessel walls
    ring_shape: object = None    # the shape of the growth rings (MeshGenTools.RingShape): ring_shape.rho(points) is the
                                 # ring coordinate, ring_shape.boundary(result.radii[k]) the outline of a ring boundary
    graph: Optional[dict] = None # the arrays of geoName-graph.npz when graphFlag is on, see graph_arrays()

    @property
    def nvertex(self):
        return len(self.vertices)

    @property
    def nridge(self):
        return len(self.ridges)


def generate(config, **changes):
    """Generate a mesh, with the same steps and output files as the template scripts of RingsPy <= 1.1.

    config :: a Config, a TOML file name, or a dict of parameters.
    changes :: parameters to override, e.g. ``generate('wood_cube.toml', seed=1, stlFlag=False)``.

    The files go to config.path/config.geoName/, which is deleted first if it exists.
    """
    if isinstance(config, (str, os.PathLike)):
        config = Config.from_toml(config)
    elif isinstance(config, dict):
        config = Config.from_dict(config)
    elif not isinstance(config, Config):
        raise TypeError("config must be a Config, a TOML file name or a dict, not {}".format(type(config).__name__))
    if changes:
        config = config.replace(**changes)
    c = config.resolved()

    # close the preview figure afterwards, also on errors (generate() may run many times in one session)
    open_figures = set(plt.get_fignums())
    try:
        return _generate(c)
    finally:
        for num in set(plt.get_fignums()) - open_figures:
            plt.close(num)


def _generate(c):
    startTime = time.time()
    geoName, path = c.geoName, c.path
    rng = None if c.seed is None else np.random.default_rng(c.seed)

    # ==================================================================
    # Remove directory/files if exists already
    output_dir = Path(path, geoName)
    shutil.rmtree(output_dir, ignore_errors=True)
    os.makedirs(output_dir, exist_ok=True)

    # ==================================================================
    # Place cells with a specific radial growth pattern
    sites, radii, site_type, hardwood_counts = place_cells(c, rng, with_types=True)
    hardwood = bool(np.any(site_type != rpgen.CELL_TYPES['fiber']))
    shape = c.ring_shape()

    placementTime = time.time()
    nParticles = len(sites)
    print('{:d} particles/cells placed in {:.3f} seconds'.format(nParticles, (placementTime - startTime)))

    if c.save_sites:
        np.save(Path(output_dir, geoName + '_sites.npy'), sites)
        np.save(Path(output_dir, geoName + '_radii.npy'), radii)
        if hardwood:
            np.save(Path(output_dir, geoName + '_types.npy'), site_type)

    # ==================================================================
    # Clipping box (boundaries) of the model
    x_min, x_max, y_min, y_max, boundaries, boundary_points, boundarylines = \
        rpgen.Clipping_Box(c.box_shape, c.box_center, c.box_size, c.boundaryFlag, polygon=c.box_polygon)
    box_corners = boundary_points

    # Visualize the original 2D Voronoi diagram
    vor = Voronoi(sites[:, 0:2])
    voronoi_plot_2d(vor, show_vertices=False, line_width=0.5, line_alpha=0.6, point_size=2)
    fig = plt.gcf()
    plt.xlim(x_min-0.1*abs(x_max-x_min), x_max+0.1*abs(x_max-x_min))
    plt.ylim(y_min-0.1*abs(y_max-y_min), y_max+0.1*abs(y_max-y_min))
    plt.plot(boundary_points[:, 0], boundary_points[:, 1], 'bo')

    voronoiTime = time.time()
    print('Original Voronoi tessellation generated in {:.3f} seconds'.format(voronoiTime - placementTime))

    ax = plt.gca()
    boundarylines = patches.Polygon(boundary_points, closed=True, linewidth=2, edgecolor='k', facecolor='none')
    ax.set_aspect('equal', adjustable='box')
    ax.add_patch(boundarylines)

    # ==================================================================
    # Rebuild the Voronoi mesh
    if c.merge_operation:
        [voronoi_vertices, finite_ridges, boundary_points, finite_ridges_new,
         boundary_ridges_new, nvertex, nvertices_in, nfinite_ridge, nboundary_ridge,
         nboundary_pts, nboundary_pts_featured, voronoi_ridges, nridge] = \
            rpgen.RebuildVoronoi_merge(vor, sites, boundaries, c.generation_center, x_min, x_max, y_min, y_max,
                                       c.box_center, c.box_shape, c.merge_tol, c.boundaryFlag)
        RebuildvorTime = time.time()
        print('Voronoi tessellation rebuilt and merged in {:.3f} seconds'.format(RebuildvorTime - voronoiTime))
    else:
        [voronoi_vertices, finite_ridges, boundary_points, finite_ridges_new,
         boundary_ridges_new, nvertex, nvertices_in, nfinite_ridge, nboundary_ridge,
         nboundary_pts, nboundary_pts_featured, voronoi_ridges, nridge] = \
            rpgen.RebuildVoronoi(vor, sites, boundaries, c.generation_center, x_min, x_max, y_min, y_max,
                                 c.box_center, c.box_shape, c.boundaryFlag)
        RebuildvorTime = time.time()
        print('Voronoi tessellation rebuilt in {:.3f} seconds'.format(RebuildvorTime - voronoiTime))

    # ==================================================================
    # Insert mid and quarter points on the Voronoi ridges (potential failure positions on cell walls)
    [all_pts_2D, all_ridges, npt_per_layer, npt_per_layer_normal, npt_per_layer_vtk] = \
        rpgen.RidgeMidQuarterPts(voronoi_vertices, nvertex, nvertices_in, voronoi_ridges,
                                 finite_ridges_new, boundary_ridges_new, nfinite_ridge,
                                 nboundary_ridge, nboundary_pts, nboundary_pts_featured)

    # ==================================================================
    # Generate a file for vertices and ridges info
    [all_vertices_2D, max_wings, flattened_all_vertices_2D, all_ridges] = \
        rpgen.VertexandRidgeinfo(all_pts_2D, all_ridges, npt_per_layer, npt_per_layer_normal,
                                 npt_per_layer_vtk, nridge, geoName, radii, c.generation_center,
                                 c.cellwallthickness_sparse, c.cellwallthickness_dense, path=path,
                                 transition=c.transition, transition_width=c.transition_width, shape=shape)

    # ==================================================================
    # Growth ring of every vertex (grain) and ridge (cell wall), for the ParaView and Abaqus files: a wall
    # is labelled by its second vertex, whose thickness it takes, and a connector by its wall (ConnectorMeshFile)
    vertex_band = rpgen.RingBands(all_pts_2D[0:npt_per_layer], radii, c.generation_center, shape=shape)
    ridge_band = vertex_band[np.asarray(all_ridges)[0:nridge, 1]]

    # the walls and vertices of the vessels and rays of a hardwood (all False for a softwood)
    ridge_vessel, ridge_ray, vertex_vessel, vertex_ray = \
        rpgen.CellTypeLabels(sites, site_type, all_pts_2D, all_ridges, nridge, npt_per_layer)
    cell_labels = {'site_type': site_type, 'vertex_vessel': vertex_vessel, 'vertex_ray': vertex_ray,
                   'ridge_vessel': ridge_vessel, 'ridge_ray': ridge_ray, 'ring_shape': shape}
    extra_labels = {'Vessel': (vertex_vessel, ridge_vessel), 'Ray': (vertex_ray, ridge_ray)} if hardwood else None
    hardwood_log = None
    if hardwood or c.vessel_fraction > 0 or c.ray_spacing is not None:
        hardwood_log = {'vessel_fraction': c.vessel_fraction, 'vessel_size': c.vessel_size,
                        'vessel_arrangement': c.vessel_arrangement, 'vessel_zone': c.vessel_zone,
                        'ray_spacing': c.ray_spacing, 'ray_width': c.ray_width, 'ray_cell_length': c.ray_cell_length}
        hardwood_log.update(hardwood_counts)
        hardwood_log.update({'nvessel_walls': int(ridge_vessel.sum()), 'nray_walls': int(ridge_ray.sum())})

    if c.fiber_length is not None:
        # fibers of various lengths: every grain its own sequence of fibers (lattice module)
        return _generate_fibers(c, rng, output_dir, fig, sites, radii, nParticles, x_min, x_max, y_min, y_max,
                                boundary_points, box_corners, boundaries, all_pts_2D, all_ridges, npt_per_layer,
                                npt_per_layer_normal, npt_per_layer_vtk, nvertex, nvertices_in, nboundary_pts_featured,
                                nridge, all_vertices_2D, max_wings, flattened_all_vertices_2D, vertex_band, ridge_band,
                                [startTime, placementTime, voronoiTime, RebuildvorTime], cell_labels, extra_labels,
                                hardwood_log)

    # ==================================================================
    # Extrude in the parallel-to-grain (longitudinal) direction
    [IGAvertices, vertex_connectivity, beam_connectivity_original, nbeam_total,
     beam_connectivity, nbeamElem, nctrlpt_per_beam, nlayers, segment_length, connector_t_connectivity,
     connector_t_bot_connectivity, connector_t_top_connectivity,
     connector_t_reg_connectivity, connector_l_connectivity, nconnector_t_per_beam,
     nconnector_t_per_grain, nconnector_t, nconnector_l, nconnector_total,
     theta, z_coord, connector_l_vertex_dict] = \
        rpgen.GenerateBeamElement(NURBS_DEGREE, c.nsegments, c.theta_min, c.theta_max, c.z_min, c.z_max,
                                  c.long_connector_ratio, npt_per_layer, voronoi_vertices,
                                  nvertex, voronoi_ridges, nridge, c.generation_center,
                                  all_vertices_2D, max_wings, flattened_all_vertices_2D, all_ridges)

    BeamTime = time.time()
    print('{:d} beam elements generated in {:.3f} seconds'.format(nbeamElem, (BeamTime - RebuildvorTime)))

    # ==================================================================
    # Insert precracks
    if c.precrackFlag:
        precrack_nodes = np.array(c.precrack_nodes if c.precrack_nodes is not None
                                  else default_precrack_nodes(c, x_min))
        [precrack_elem, nconnector_t_precrack, nconnector_l_precrack] = \
            rpgen.insert_precracks(all_pts_2D, all_ridges, nridge, npt_per_layer,
                                   npt_per_layer_normal, npt_per_layer_vtk,
                                   nlayers, precrack_nodes, c.precrack_widths,
                                   float(np.max(c.cellsize_sparse)))
    else:
        precrack_elem = []
        nconnector_t_precrack = 0
        nconnector_l_precrack = 0

    # ==================================================================
    # Calculate mesh info
    [ConnMeshData, conn_l_tangents] = \
        rpgen.ConnectorMeshFile(geoName, IGAvertices, connector_t_bot_connectivity,
                                connector_t_reg_connectivity, connector_t_top_connectivity,
                                connector_l_connectivity, all_vertices_2D,
                                max_wings, flattened_all_vertices_2D, c.nsegments, segment_length,
                                nctrlpt_per_beam, theta, nridge, connector_l_vertex_dict, path=path)

    # ==================================================================
    # Calculate model properties
    [mass, bulk_volume, bulk_density, porosity] = \
        rpgen.ModelInfo(c.box_shape, boundary_points, c.z_min, c.z_max, c.skeleton_density, ConnMeshData,
                        polygon=box_corners if c.box_shape in rpgen.POLYGON_SHAPES else None)

    # ==================================================================
    # Bezier extraction
    knotVec = rpgen.BezierExtraction(NURBS_DEGREE, nbeam_total)
    npatch = beam_connectivity_original.shape[0]

    rpgen.BezierBeamFile(geoName, NURBS_DEGREE, nctrlpt_per_beam, nconnector_t_per_beam, npatch, knotVec, path=path)

    # ==================================================================
    # Generate Paraview visulization files
    rpgen.VisualizationFiles(geoName, NURBS_DEGREE, nlayers, npt_per_layer_vtk, all_pts_2D,
                             segment_length, theta, z_coord, c.nsegments, nridge,
                             voronoi_ridges, c.generation_center, all_ridges, nvertex, nconnector_t,
                             nconnector_l, nctrlpt_per_beam, ConnMeshData, conn_l_tangents,
                             all_vertices_2D, max_wings, flattened_all_vertices_2D, path=path,
                             vertex_band=vertex_band, ridge_band=ridge_band, extra_labels=extra_labels)

    fig.savefig(Path(output_dir, geoName + '.png'))

    # ==================================================================
    # Generate 3D model files
    if c.stlFlag:
        rpgen.StlModelFile(geoName, path=path)

    # ==================================================================
    # Generate input files for numerical simulations (Config checks that inpType is Abaqus)
    if c.inpFlag:
        rpgen.AbaqusFile(geoName, NURBS_DEGREE, npatch, IGAvertices, beam_connectivity,
                         connector_t_bot_connectivity, connector_t_reg_connectivity,
                         connector_t_top_connectivity, connector_l_connectivity, nbeamElem,
                         nconnector_t, nconnector_l, nconnector_t_precrack, nconnector_l_precrack,
                         segment_length, c.long_connector_ratio,
                         x_max, x_min, y_max, y_min, z_coord, c.box_shape, c.box_size,
                         float(np.mean(c.cellwallthickness_sparse)), float(np.mean(c.cellwallthickness_dense)),
                         c.merge_operation, c.merge_tol,
                         c.precrackFlag, precrack_elem, path=path, vertex_band=vertex_band, ridge_band=ridge_band,
                         wing_ridge=flattened_all_vertices_2D[0, :].astype(int), extra_labels=extra_labels)

    # ==================================================================
    # The 2D cell mesh with its ring labels, and as a graph
    vertices = np.array(all_pts_2D[0:npt_per_layer], dtype=float)
    ridges = np.array(all_ridges, dtype=int)[0:nridge, 0:2]
    vertex_ring, vertex_dense = rpgen.RingLabels(vertex_band)
    ridge_ring, ridge_dense = rpgen.RingLabels(ridge_band)
    vertex_boundary = rpgen.OnBoundary(vertices, box_corners)
    result = Result(config=c, output_dir=output_dir, files=[], sites=sites, radii=radii,
                    nParticles=nParticles, nbeamElem=nbeamElem, mass=mass, bulk_volume=bulk_volume,
                    bulk_density=bulk_density, porosity=porosity, vertices=vertices, ridges=ridges,
                    vertex_ring=vertex_ring, vertex_dense=vertex_dense, vertex_boundary=vertex_boundary,
                    vertex_thickness=np.array(all_vertices_2D[:, 3 + 3*max_wings], dtype=float),  # width of the first wing
                    ridge_ring=ridge_ring, ridge_dense=ridge_dense,
                    ridge_boundary=vertex_boundary[ridges[:, 0]] & vertex_boundary[ridges[:, 1]], **cell_labels)
    if c.graphFlag:
        result.graph = graph_arrays(result)
        rpgen.GraphFile(geoName, result.graph, path=path)

    FileTime = time.time()
    print('Files generated in {:.3f} seconds'.format(FileTime - BeamTime))

    # ==================================================================
    # Generate log file for the generation
    rpgen.LogFile(geoName, c.iter_max, c.r_min, c.r_max, c.nrings, c.width_heart, c.width_sparse, c.width_dense,
                  c.generation_center, c.box_shape, c.box_center, c.box_size, x_min, x_max, y_min, y_max,
                  c.cellsize_sparse, c.cellsize_dense, c.cellwallthickness_sparse, c.cellwallthickness_dense,
                  c.merge_operation, c.merge_tol, c.precrackFlag, c.precrack_widths, c.boundaryFlag,
                  segment_length, c.theta_min, c.theta_max, c.z_min, c.z_max, c.long_connector_ratio,
                  NURBS_DEGREE, nctrlpt_per_beam, nconnector_t_precrack, nconnector_l_precrack,
                  nParticles, nbeamElem, c.skeleton_density, mass, bulk_volume, bulk_density, porosity,
                  c.stlFlag, c.inpFlag, c.inpType, c.radial_growth_rule,
                  startTime, placementTime, voronoiTime, RebuildvorTime, BeamTime, FileTime,
                  path=path, seed=c.seed, box_polygon=c.box_polygon,
                  transition=c.transition, transition_width=c.transition_width, radii=radii, graphFlag=c.graphFlag,
                  hardwood=hardwood_log, ring_shape=shape)

    result.files = sorted(output_dir.iterdir())
    return result


NODE_FEATURES = ['ring', 'dense', 'boundary', 'thickness', 'degree', 'radius', 'vessel', 'ray']
EDGE_FEATURES = ['length', 'angle', 'thickness', 'ring', 'dense', 'boundary', 'vessel', 'ray']


def graph_arrays(result):
    """The 2D cell mesh of a Result as a graph, in the layout of a PyTorch Geometric ``Data`` object
    (the arrays of the geoName-graph.npz file written with graphFlag):

    - ``pos`` (nvertex, 2): coordinates of the nodes, the cell vertices (one beam each in 3D)
    - ``x`` (nvertex, 8): node features, in the order of ``node_features``: growth ring (0: heart),
      dense (1 in the heart and the dense rings), boundary (1 on the clipping boundary), cell wall
      thickness at the vertex, degree (number of walls), distance to the generation center, vessel
      and ray (1 for the vertices of a vessel or of a ray cell of a hardwood)
    - ``edge_index`` (2, 2*nridge): the cell walls as edges, each in both directions
    - ``edge_attr`` (2*nridge, 8): edge features, in the order of ``edge_features``: length, angle of
      the wall from its first to its second node (radians, counter-clockwise from the x axis),
      thickness, growth ring, dense, boundary (1 for walls along the clipping boundary), vessel and
      ray (1 for the walls of a vessel or of a ray cell)
    - ``sites`` (nParticles, 3), ``radii``, ``generation_center``: the cells and growth rings behind the mesh
    - ``nParticles``, ``nbeamElem``, ``mass``, ``bulk_volume``, ``bulk_density``, ``porosity``: model properties
    - ``config``: the parameters of the run as a JSON string

    To load it: ``g = np.load(file); Data(x=torch.from_numpy(g['x']), edge_index=torch.from_numpy(g['edge_index']),
    edge_attr=torch.from_numpy(g['edge_attr']), pos=torch.from_numpy(g['pos']))``.
    """
    pos = np.asarray(result.vertices, dtype=float)
    ridges = np.asarray(result.ridges, dtype=np.int64)
    degree = np.bincount(ridges.ravel(), minlength=len(pos))
    radius = np.linalg.norm(pos - np.asarray(result.config.generation_center, dtype=float), axis=1)
    x = np.column_stack((result.vertex_ring, result.vertex_dense, result.vertex_boundary,
                         result.vertex_thickness, degree, radius, result.vertex_vessel, result.vertex_ray)).astype(float)
    # every wall in both directions: (i, j) then (j, i)
    edge_index = np.hstack((ridges.T, ridges[:, ::-1].T))
    vector = pos[edge_index[1]] - pos[edge_index[0]]
    ridge_thickness = result.vertex_thickness[ridges[:, 1]]
    edge_attr = np.column_stack((np.linalg.norm(vector, axis=1), np.arctan2(vector[:, 1], vector[:, 0]),
                                 np.tile(ridge_thickness, 2), np.tile(result.ridge_ring, 2),
                                 np.tile(result.ridge_dense, 2), np.tile(result.ridge_boundary, 2),
                                 np.tile(result.ridge_vessel, 2), np.tile(result.ridge_ray, 2))).astype(float)
    return {'pos': pos, 'x': x, 'edge_index': edge_index, 'edge_attr': edge_attr,
            'node_features': np.array(NODE_FEATURES), 'edge_features': np.array(EDGE_FEATURES),
            'sites': np.asarray(result.sites, dtype=float), 'radii': np.asarray(result.radii, dtype=float),
            'generation_center': np.asarray(result.config.generation_center, dtype=float),
            'nParticles': np.int64(result.nParticles), 'nbeamElem': np.int64(result.nbeamElem),
            'mass': np.float64(result.mass), 'bulk_volume': np.float64(result.bulk_volume),
            'bulk_density': np.float64(result.bulk_density), 'porosity': np.float64(result.porosity),
            'config': np.array(json.dumps(result.config.to_dict()))}


def _generate_fibers(c, rng, output_dir, fig, sites, radii, nParticles, x_min, x_max, y_min, y_max, boundary_points,
                     box_corners, boundaries, all_pts_2D, all_ridges, npt_per_layer, npt_per_layer_normal,
                     npt_per_layer_vtk, nvertex, nvertices_in, nboundary_pts_featured, nridge, all_vertices_2D,
                     max_wings, flattened_all_vertices_2D, vertex_band, ridge_band, times, cell_labels, extra_labels,
                     hardwood_log):
    """The rest of _generate for fiber_length: the cells end at random heights, tapering to their ends, and
    the 2D mesh is extruded into fibers of various lengths (one fiber end on a grain per end of its cells)"""
    startTime, placementTime, voronoiTime, RebuildvorTime = times
    geoName, path = c.geoName, c.path
    nelem = 2*c.nsegments
    nhalf = 2*nelem + 1
    h = (c.z_max - c.z_min)/nelem
    segment_length = 2*h
    theta_half = np.linspace(c.theta_min, c.theta_max, nhalf)
    z_half = np.linspace(c.z_min, c.z_max, nhalf)
    kind = np.zeros(nvertex, dtype=int)
    kind[nvertices_in:] = 1
    kind[nvertex - nboundary_pts_featured:] = 2
    # the cells meeting at each vertex, the heights at which the cells end, the fiber ends of the grains
    vertex_cells = rplat.cells_of_vertices(all_pts_2D[:nvertex], all_ridges, sites, boundaries)
    cells = sorted(set().union(*vertex_cells))
    cell_ends = rplat.cell_end_heights(cells, c.z_min, c.z_max, c.fiber_length, c.fiber_length_cov, rng=rng)
    fiber_ends = rplat.grain_ends_from_cells(vertex_cells, cell_ends, c.z_min, h, nelem)
    # the vertices on every half-layer: twisted, and drawn in toward the cells tapering to their ends
    vertex_layers = rpgen._rotate_layers(all_pts_2D[:nvertex], theta_half, z_half, c.generation_center)[:, :, 0:2]
    vertex_layers = rplat.taper_vertices(vertex_layers, z_half, vertex_cells, cell_ends, c.taper_length, c.taper_ratio,
                                         kind=kind, boundaries=boundaries)
    layer_points = rpgen.LayerPoints(vertex_layers, all_ridges, npt_per_layer_vtk)
    lat = rplat.build_lattice(all_ridges, nvertex, nridge, all_vertices_2D, max_wings, flattened_all_vertices_2D,
                              layer_points, c.z_min, c.z_max, nelem, c.long_connector_ratio, fiber_ends)
    nbeamElem = lat.nbeamElem
    BeamTime = time.time()
    print('{:d} beam elements of {:d} fibers generated in {:.3f} seconds'.format(nbeamElem, lat.npatch, BeamTime - RebuildvorTime))

    if c.precrackFlag:
        precrack_nodes = np.array(c.precrack_nodes if c.precrack_nodes is not None
                                  else default_precrack_nodes(c, x_min))
        [precrack_elem, nconnector_t_precrack, nconnector_l_precrack] = \
            rpgen.insert_precracks(all_pts_2D, all_ridges, nridge, npt_per_layer, npt_per_layer_normal,
                                   npt_per_layer_vtk, nhalf, precrack_nodes, c.precrack_widths,
                                   float(np.max(c.cellsize_sparse)))
    else:
        precrack_elem, nconnector_t_precrack, nconnector_l_precrack = [], 0, 0

    rplat.write_connector_file(geoName, lat, path=path)
    [mass, bulk_volume, bulk_density, porosity] = \
        rpgen.ModelInfo(c.box_shape, boundary_points, c.z_min, c.z_max, c.skeleton_density, lat.Meshdata,
                        polygon=box_corners if c.box_shape in rpgen.POLYGON_SHAPES else None)
    rplat.write_iga_file(geoName, lat, path=path)
    wing_ridge = flattened_all_vertices_2D[0, :].astype(int)
    rplat.write_vtu_files(geoName, lat, layer_points, all_ridges, nvertex, path=path, vertex_band=vertex_band,
                          ridge_band=ridge_band, extra_labels=extra_labels)
    fig.savefig(Path(output_dir, geoName + '.png'))
    if c.stlFlag:
        rpgen.StlModelFile(geoName, path=path)
    if c.inpFlag:
        rpgen.AbaqusFile(geoName, NURBS_DEGREE, lat.npatch, lat.nodes, lat.beam_connectivity,
                         lat.connector_t_bot, lat.connector_t_reg, lat.connector_t_top, lat.connector_l, nbeamElem,
                         lat.nconnector_t, lat.nconnector_l, nconnector_t_precrack, nconnector_l_precrack,
                         segment_length, c.long_connector_ratio, x_max, x_min, y_max, y_min, lat.nodes[:, 2],
                         c.box_shape, c.box_size, float(np.mean(c.cellwallthickness_sparse)),
                         float(np.mean(c.cellwallthickness_dense)), c.merge_operation, c.merge_tol, c.precrackFlag,
                         precrack_elem, path=path, connector_t_ridge=lat.connector_t_ridge,
                         element_bands=rplat.element_bands(lat, vertex_band, ridge_band, wing_ridge),
                         element_labels={name: rplat.element_bands(lat, vv, rv, wing_ridge)
                                         for name, (vv, rv) in (extra_labels or {}).items()})

    # the 2D cell mesh with its ring labels, and as a graph
    vertices = np.array(all_pts_2D[0:npt_per_layer], dtype=float)
    ridges = np.array(all_ridges, dtype=int)[0:nridge, 0:2]
    vertex_ring, vertex_dense = rpgen.RingLabels(vertex_band)
    ridge_ring, ridge_dense = rpgen.RingLabels(ridge_band)
    vertex_boundary = rpgen.OnBoundary(vertices, box_corners)
    result = Result(config=c, output_dir=output_dir, files=[], sites=sites, radii=radii,
                    nParticles=nParticles, nbeamElem=nbeamElem, mass=mass, bulk_volume=bulk_volume,
                    bulk_density=bulk_density, porosity=porosity, vertices=vertices, ridges=ridges,
                    vertex_ring=vertex_ring, vertex_dense=vertex_dense, vertex_boundary=vertex_boundary,
                    vertex_thickness=np.array(all_vertices_2D[:, 3 + 3*max_wings], dtype=float),
                    ridge_ring=ridge_ring, ridge_dense=ridge_dense,
                    ridge_boundary=vertex_boundary[ridges[:, 0]] & vertex_boundary[ridges[:, 1]], **cell_labels)
    if c.graphFlag:
        result.graph = graph_arrays(result)
        rpgen.GraphFile(geoName, result.graph, path=path)
    FileTime = time.time()
    print('Files generated in {:.3f} seconds'.format(FileTime - BeamTime))

    fibers = {'fiber_length': c.fiber_length, 'fiber_length_cov': c.fiber_length_cov,
              'taper_length': c.taper_length, 'taper_ratio': c.taper_ratio, 'nelem': nelem,
              'ncells': len(cells), 'ncell_ends': int(sum(len(v) for v in cell_ends.values()))}
    fibers.update(rplat.fiber_statistics(fiber_ends, nelem, h))
    rpgen.LogFile(geoName, c.iter_max, c.r_min, c.r_max, c.nrings, c.width_heart, c.width_sparse, c.width_dense,
                  c.generation_center, c.box_shape, c.box_center, c.box_size, x_min, x_max, y_min, y_max,
                  c.cellsize_sparse, c.cellsize_dense, c.cellwallthickness_sparse, c.cellwallthickness_dense,
                  c.merge_operation, c.merge_tol, c.precrackFlag, c.precrack_widths, c.boundaryFlag,
                  segment_length, c.theta_min, c.theta_max, c.z_min, c.z_max, c.long_connector_ratio,
                  NURBS_DEGREE, 2*NURBS_DEGREE + 1, nconnector_t_precrack, nconnector_l_precrack,
                  nParticles, nbeamElem, c.skeleton_density, mass, bulk_volume, bulk_density, porosity,
                  c.stlFlag, c.inpFlag, c.inpType, c.radial_growth_rule,
                  startTime, placementTime, voronoiTime, RebuildvorTime, BeamTime, FileTime,
                  path=path, seed=c.seed, box_polygon=c.box_polygon,
                  transition=c.transition, transition_width=c.transition_width, radii=radii, graphFlag=c.graphFlag,
                  fibers=fibers, hardwood=hardwood_log, ring_shape=cell_labels['ring_shape'])
    result.files = sorted(output_dir.iterdir())
    return result


def place_cells(c, rng=None, with_types=False):
    """The cell sites (one row (x, y, cell size) per cell) and growth ring radii of a resolved Config:
    placed with its radial growth rule, or loaded from its *_sites.npy file.

    with_types: also return the type of every site (0 fiber, 1 vessel, 2 ray cell, see
    MeshGenTools.CELL_TYPES) and a dict of the numbers of vessels, rays and ray cells placed (empty
    without hardwood features)."""
    geoName, path = c.geoName, c.path
    site_type, counts = None, {}
    shape = c.ring_shape()
    if c.radial_growth_rule == 'binary' and (c.vessel_fraction > 0 or c.ray_spacing is not None):
        # hardwood: the vessels and rays first, then the fibers packed around them
        sites, radii, site_type, counts = rpgen.CellPlacement_Hardwood(
            c.generation_center, c.r_max, c.r_min, c.nrings, c.width_heart, c.width_sparse, c.width_dense,
            c.cellsize_sparse, c.cellsize_dense, c.iter_max, c.print_interval, rng=rng,
            transition=c.transition, transition_width=c.transition_width, vessel_fraction=c.vessel_fraction,
            vessel_size=c.vessel_size, vessel_arrangement=c.vessel_arrangement, vessel_zone=c.vessel_zone,
            ray_spacing=c.ray_spacing, ray_width=c.ray_width, ray_cell_length=c.ray_cell_length, shape=shape)
    elif c.radial_growth_rule == 'binary':
        # binary radial growth pattern (e.g. wood microstructure with earlywood-latewood alternations)
        sites, radii = rpgen.CellPlacement_Binary(c.generation_center, c.r_max, c.r_min, c.nrings, c.width_heart,
                                                  c.width_sparse, c.width_dense, c.cellsize_sparse, c.cellsize_dense,
                                                  c.iter_max, c.print_interval, rng=rng,
                                                  transition=c.transition, transition_width=c.transition_width,
                                                  shape=shape)
    elif c.radial_growth_rule == 'binary_lloyd':
        # binary with Lloyd's algorithm (more regular cell shapes)
        sites, radii = rpgen.CellPlacement_Binary_Lloyd(geoName, path, c.generation_center, c.r_max, c.r_min,
                                                        c.nrings, c.width_heart, c.width_sparse, c.width_dense,
                                                        c.cellsize_sparse, c.cellsize_dense, c.iter_max,
                                                        c.print_interval, rng=rng,
                                                        transition=c.transition, transition_width=c.transition_width,
                                                        shape=shape)
    elif c.radial_growth_rule == 'regular_hexagonal':
        # hexagonal honeycomb-like geometry
        sites, radii = rpgen.CellPlacement_Honeycomb(c.generation_center, c.r_max, c.r_min, c.nrings,
                                                     c.box_center, c.box_size, c.width_heart,
                                                     c.width_sparse, c.width_dense,
                                                     c.cellsize_sparse, c.cellsize_dense,
                                                     c.cellwallthickness_sparse, c.cellwallthickness_dense,
                                                     c.iter_max, c.print_interval, rng=rng)
    else:
        # load saved cell sites and radii data (and the cell types of a hardwood, if saved)
        sites, radii, site_type = rpgen.ReadSavedSites(Path(c.sites_dir or '.'), c.radial_growth_rule, with_types=True)
    if not with_types:
        return sites, radii
    if site_type is None:
        site_type = np.zeros(len(sites), dtype=int)
    return sites, radii, np.asarray(site_type, dtype=int), counts


def default_precrack_nodes(c, x_min):
    """The precrack of the template scripts: a horizontal line through box_center, starting
    0.12*box_size from the left edge of the box and 0.1*box_size long"""
    x_indent = x_min + c.box_size*0.120
    x_precrack = x_indent + c.box_size*0.1
    y_precrack = c.box_center[1]
    return [[x_indent, y_precrack, x_precrack, y_precrack]]
