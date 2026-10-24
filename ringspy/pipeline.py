# -*- coding: utf-8 -*-
"""
The full mesh generation: cell placement, Voronoi tessellation and clipping,
longitudinal extrusion, and the output files

author: king_yin3613
email: haoyin2022@u.northwestern.edu
"""

import os
import shutil
import time
from dataclasses import dataclass
from pathlib import Path
from typing import List

import matplotlib.patches as patches
import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial import Voronoi, voronoi_plot_2d

from . import MeshGenTools as rpgen
from .config import Config

NURBS_DEGREE = 2


@dataclass
class Result:
    """What generate() made: the output folder and files, the cell sites, and model properties"""
    config: Config               # with the defaults filled in (Config.resolved())
    output_dir: Path
    files: List[Path]
    sites: np.ndarray            # cell sites, one row (x, y, cell size) per cell
    radii: np.ndarray            # radii of the growth rings
    nParticles: int
    nbeamElem: int
    mass: float
    bulk_volume: float
    bulk_density: float
    porosity: float


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
    if c.radial_growth_rule == 'binary':
        # binary radial growth pattern (e.g. wood microstructure with earlywood-latewood alternations)
        sites, radii = rpgen.CellPlacement_Binary(c.generation_center, c.r_max, c.r_min, c.nrings, c.width_heart,
                                                  c.width_sparse, c.width_dense, c.cellsize_sparse, c.cellsize_dense,
                                                  c.iter_max, c.print_interval, rng=rng)
    elif c.radial_growth_rule == 'binary_lloyd':
        # binary with Lloyd's algorithm (more regular cell shapes)
        sites, radii = rpgen.CellPlacement_Binary_Lloyd(geoName, path, c.generation_center, c.r_max, c.r_min,
                                                        c.nrings, c.width_heart, c.width_sparse, c.width_dense,
                                                        c.cellsize_sparse, c.cellsize_dense, c.iter_max,
                                                        c.print_interval, rng=rng)
    elif c.radial_growth_rule == 'regular_hexagonal':
        # hexagonal honeycomb-like geometry
        sites, radii = rpgen.CellPlacement_Honeycomb(c.generation_center, c.r_max, c.r_min, c.nrings,
                                                     c.box_center, c.box_size, c.width_heart,
                                                     c.width_sparse, c.width_dense,
                                                     c.cellsize_sparse, c.cellsize_dense,
                                                     c.cellwallthickness_sparse, c.cellwallthickness_dense,
                                                     c.iter_max, c.print_interval, rng=rng)
    else:
        # load saved cell sites and radii data
        sites, radii = rpgen.ReadSavedSites(Path(c.sites_dir or '.'), c.radial_growth_rule)

    placementTime = time.time()
    nParticles = len(sites)
    print('{:d} particles/cells placed in {:.3f} seconds'.format(nParticles, (placementTime - startTime)))

    if c.save_sites:
        np.save(Path(output_dir, geoName + '_sites.npy'), sites)
        np.save(Path(output_dir, geoName + '_radii.npy'), radii)

    # ==================================================================
    # Clipping box (boundaries) of the model
    x_min, x_max, y_min, y_max, boundaries, boundary_points, boundarylines = \
        rpgen.Clipping_Box(c.box_shape, c.box_center, c.box_size, c.boundaryFlag)

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
                                 c.cellwallthickness_sparse, c.cellwallthickness_dense, path=path)

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
                                   c.cellsize_sparse)
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
        rpgen.ModelInfo(c.box_shape, boundary_points, c.z_min, c.z_max, c.skeleton_density, ConnMeshData)

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
                             all_vertices_2D, max_wings, flattened_all_vertices_2D, path=path)

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
                         c.cellwallthickness_sparse, c.cellwallthickness_dense,
                         c.merge_operation, c.merge_tol,
                         c.precrackFlag, precrack_elem, path=path)

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
                  path=path, seed=c.seed)

    return Result(config=c, output_dir=output_dir, files=sorted(output_dir.iterdir()), sites=sites, radii=radii,
                  nParticles=nParticles, nbeamElem=nbeamElem, mass=mass, bulk_volume=bulk_volume,
                  bulk_density=bulk_density, porosity=porosity)


def default_precrack_nodes(c, x_min):
    """The precrack of the template scripts: a horizontal line through box_center, starting
    0.12*box_size from the left edge of the box and 0.1*box_size long"""
    x_indent = x_min + c.box_size*0.120
    x_precrack = x_indent + c.box_size*0.1
    y_precrack = c.box_center[1]
    return [[x_indent, y_precrack, x_precrack, y_precrack]]
