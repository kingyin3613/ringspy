# -*- coding: utf-8 -*-
"""
RingsPy
a Python Voronoi diagrams-based geometric generation tool that generates 
3D meshes of prismatic cellular solids with radial growth rules.

author: king_yin3613
email: haoyin2022@u.northwestern.edu
"""

import numpy as np
from scipy.spatial import Voronoi, voronoi_plot_2d
from scipy.spatial import Delaunay, delaunay_plot_2d
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.tri as mtri
import os
import time
import RingsPy.MeshGenTools as RPgen
import shutil
from pathlib import Path

def test_hexalattice_tri():
    # performance check only
    startTime = time.time()

    # ==================================================================
    # Input parameters
    geoName = 'honeycomb_tri'
    path = 'meshes'

    iter_max = 100
    print_interval = 500

    # length unit: mm
    r_min = 0   # inner radius of wood log
    r_max = 2   # outer radius of wood log
    nrings = 4 # number of rings
    width_heart = 0.3*(r_max-r_min)/nrings # heart wood ring width
    width_early = 0.7*(r_max-r_min)/nrings # early wood ring width
    width_late = 0.3*(r_max-r_min)/nrings # late wood ring width
    log_center = (0,1) # coordinates of log center in the global system of reference
    box_center = (0,0) # coordinates of box center in the global system of reference
    box_size   = 1 # side length of the clipping box
    
    cellsize_early = 0.2
    cellsize_late = 0.2
    cellwallthickness_early = 0.010
    cellwallthickness_late = 0.006
    cellangle = 0.0

    merge_operation = 'off'
    merge_tol = 0.01

    precrackFlag = 'off'
    precrack_widths = 0.1

    boundaryFlag = 'on'
    STLFlag = 'on'

    # ==================================================================
    # Clipping box (Boundaries) of the model

    # Triangle box
    #----------------------------------
    box_shape = 'triangle'
    x_0 = box_center[0] - box_size
    x_1 = box_center[0]
    x_2 = box_center[0] + box_size
    y_0 = box_center[1]
    y_1 = box_center[1] + box_size*np.sqrt(3)
    boundary_points = np.array([[x_0, y_0], [x_2, y_0], [x_1, y_1]])
    
    l0 = [(x_0, y_0), (x_2, y_0)]
    l1 = [(x_2, y_0), (x_1, y_1)]
    l2 = [(x_1, y_1), (x_0, y_0)]
    
    boundaries = [('bottom', l0), ('top-right', l1), ('top-left', l2)]
    
    if boundaryFlag in ['on','On','Y','y','Yes','yes']:
        boundarylines = patches.RegularPolygon((box_center[0]+0,box_center[1]+np.sqrt(3)/3*box_size),numVertices=3,radius=2*np.sqrt(3)/3*box_size,orientation=0,linewidth=2,edgecolor='k',facecolor='none')
    else:
        boundarylines = patches.RegularPolygon((box_center[0]+0,box_center[1]+np.sqrt(3)/3*box_size),numVertices=3,radius=2*np.sqrt(3)/3*box_size,orientation=0,linewidth=2,linestyle='-.',edgecolor='k',facecolor='none')
    
    x_min = x_0
    x_max = x_2
    y_min = y_0
    y_max = y_1
    # ----------------------------------
    
    # if precracked
    x_indent_size = box_size*0.120
    y_indent_size = box_size*0.125
    x_precrack_size = box_size*0.1
    y_precrack_size = box_size*0.02
    x_indent = x_min + x_indent_size
    y_indent_min = box_center[1] - y_indent_size/2
    y_indent_max = box_center[1] + y_indent_size/2
    x_precrack = x_indent + x_precrack_size
    y_precrack = box_center[1]
    
    # ==================================================================
    # Remove directory/files if exists already
    try:
        shutil.rmtree(Path(path, geoName))
    except:
        pass

    if not os.path.exists(Path(path, geoName)):
        os.makedirs(Path(path, geoName))
        
    # ==================================================================
    # Random cell site placement with a radial growth pattern

    # ----------------------------------
    # hexagonal honeycomb lattice
    sites,radius = RPgen.CellPlacement_Honeycomb(log_center,r_max,r_min,nrings,\
                                    box_center,box_size,width_heart,\
                                    width_early,width_late,\
                                    cellsize_early,cellsize_late,\
                                    cellwallthickness_early,cellwallthickness_late,\
                                    cellangle,iter_max,print_interval)

    placementTime = time.time()
    nParticles = len(sites)

    print('{:d} particles/cells placed in {:.3f} seconds'.format(nParticles, (placementTime - startTime)))
                    
    # Visualize the original Voronoi diagram generated with the wood mesh sites    
    vor = Voronoi(sites[:,0:2])
    voronoi_plot_2d(vor, show_vertices=False,line_width=0.5, line_alpha=0.6, point_size=2)
    plt.xlim(x_min-0.1*abs(x_max-x_min), x_max+0.1*abs(x_max-x_min))
    plt.ylim(y_min-0.1*abs(y_max-y_min), y_max+0.1*abs(y_max-y_min))
    plt.plot(boundary_points[:, 0], boundary_points[:, 1], 'bo')
    # plt.show()

    voronoiTime = time.time() 
    print('Original Voronoi tessellation generated in {:.3f} seconds'.format(voronoiTime - placementTime))

    ax = plt.gca()
    boundarylines = patches.Polygon(boundary_points,closed=True,linewidth=2,edgecolor='k',facecolor='none')
    ax.set_aspect('equal', adjustable='box')
    ax.add_patch(boundarylines)

    # ==================================================================
    # Rebuild the Voronoi mesh
    if merge_operation in ['on','On','Y','y','Yes','yes']:
        [voronoi_vertices,finite_ridges,boundary_points,finite_ridges_new,\
         boundary_ridges_new,nvertex,nvertices_in,nfinite_ridge,nboundary_ridge,\
         nboundary_pts,nboundary_pts_featured,voronoi_ridges,nridge] = \
            RPgen.RebuildVoronoi_new(vor,sites,boundaries,log_center,x_min,x_max,y_min,y_max,box_center,merge_tol,boundaryFlag)
        
        RebuildvorTime = time.time() 
        print('Voronoi tessellation rebuilt and merged in {:.3f} seconds'.format(RebuildvorTime - voronoiTime))
    else:
        [voronoi_vertices,finite_ridges,boundary_points,finite_ridges_new,\
          boundary_ridges_new,nvertex,nvertices_in,nfinite_ridge,nboundary_ridge,\
          nboundary_pts,nboundary_pts_featured,voronoi_ridges,nridge] = \
            RPgen.RebuildVoronoi(vor,sites,boundaries,log_center,x_min,x_max,y_min,y_max,box_center,boundaryFlag)
        
        RebuildvorTime = time.time() 
        print('Voronoi tessellation rebuilt in {:.3f} seconds'.format(RebuildvorTime - voronoiTime))

    # ===============================================
    # Insert mid and quarter points on the Voronoi ridges (can be used as potential failure positions on cell walls)
    [all_pts_2D,all_ridges,npt_per_layer,npt_per_layer_normal,npt_per_layer_vtk] = \
        RPgen.RidgeMidQuarterPts(voronoi_vertices,nvertex,nvertices_in,voronoi_ridges,\
                           finite_ridges_new,boundary_ridges_new,nfinite_ridge,\
                           nboundary_ridge,nboundary_pts,nboundary_pts_featured)

    # ==================================================================        
    # Generate a file for the geometry info for vertices and ridges
    [all_vertices_2D, max_wings, flattened_all_vertices_2D, all_ridges] = \
        RPgen.VertexandRidgeinfo(all_pts_2D,all_ridges,\
                           npt_per_layer,npt_per_layer_normal,\
                           npt_per_layer_vtk,nridge,geoName,radius,log_center,\
                           cellwallthickness_early,cellwallthickness_late)
        
    ###############################################################################
    # Extrusion in the parallel-to-grain (longitudinal) direction
    fiberlength = 0.5*box_size
    theta_min = 0 # radian
    theta_max = 0 # radian
    z_min = 0
    z_max = box_size

    long_connector_ratio = 0.02

    NURBS_degree = 2
    nctrlpt_per_beam = 5

    [woodIGAvertices,vertex_connectivity,beam_connectivity_original,nbeam_total,\
     beam_connectivity,nbeamElem,nlayer,connector_t_connectivity,\
     connector_t_bot_connectivity,connector_t_top_connectivity,\
     connector_t_reg_connectivity,connector_l_connectivity,nconnector_t_per_beam,\
     nconnector_t_per_grain,nconnector_t,nconnector_l,nconnector_total,\
     theta,z_coord,nbeam_per_grain,connector_l_vertex_dict] = \
        RPgen.GenerateBeamElement(NURBS_degree,nctrlpt_per_beam,\
                            fiberlength,theta_min,theta_max,z_min,z_max,\
                            long_connector_ratio,npt_per_layer,voronoi_vertices,\
                            nvertex,voronoi_ridges,nridge,log_center,\
                            all_vertices_2D,max_wings,flattened_all_vertices_2D,all_ridges)

    BeamTime = time.time() 
    print('{:d} beam elements generated in {:.3f} seconds'.format(nbeamElem, (BeamTime - RebuildvorTime)))

    # ==================================================================
    # Insertion of precracks
    if precrackFlag in ['on','On','Y','y','Yes','yes']:
        # pre-crack node list: [c1x1 c1y1 c1x2 c1y2 ...]
        precrack_nodes = np.array([[x_indent, y_precrack, x_precrack, y_precrack]])
        [precracked_elem,nconnector_t_precrack,nconnector_l_precrack] = \
            RPgen.insert_precracks(all_pts_2D,all_ridges,nridge,npt_per_layer,\
                                     npt_per_layer_normal,npt_per_layer_vtk,\
                                     nlayer,precrack_nodes,precrack_widths,\
                                     cellsize_early)
    else:
        precrack_nodes = []
        precracked_elem = []
        nconnector_t_precrack = 0
        nconnector_l_precrack = 0

    # ==================================================================
    # Connector Calculations
    height_connector_t = fiberlength/4

    ConnMeshData = RPgen.ConnectorMeshFile(geoName,woodIGAvertices,connector_t_bot_connectivity,\
                     connector_t_reg_connectivity,connector_t_top_connectivity,\
                     height_connector_t,connector_l_connectivity,all_vertices_2D,\
                     max_wings,flattened_all_vertices_2D,nbeam_per_grain,nridge,\
                     connector_l_vertex_dict)

    # ==================================================================
    # Bezier extraction 
    knotVec = RPgen.BezierExtraction(NURBS_degree,nctrlpt_per_beam,nbeam_total)
    npatch = beam_connectivity_original.shape[0]

    mkBezierBeamFile = RPgen.BezierBeamFile(geoName,NURBS_degree,nctrlpt_per_beam,\
                       nconnector_t_per_beam,npatch,knotVec)

    # ==================================================================
    # Generate visulization data 
    RPgen.VisualizationFiles(geoName,NURBS_degree,nlayer,npt_per_layer_vtk,all_pts_2D,\
                       fiberlength,theta,z_coord,nbeam_per_grain,nridge,\
                       voronoi_ridges,log_center,all_ridges,nvertex,nconnector_t,\
                       nconnector_l,nctrlpt_per_beam,ConnMeshData,all_vertices_2D,\
                       max_wings,flattened_all_vertices_2D)
        
    plt.savefig(Path('meshes/' + geoName + '/' + geoName + '.png'))

    FileTime = time.time() 
    print('Files generated in {:.3f} seconds'.format(FileTime - BeamTime))

    # ==================================================================
    # Generate 3D printing data
    if STLFlag in ['on','On','Y','y','Yes','yes']:
        RPgen.StlModelFile(geoName)

    # ==================================================================
    # Generate log file for the mesh generation
    RPgen.LogFile(geoName,iter_max,r_min,r_max,nrings,width_heart,width_early,width_late,\
            log_center,box_shape,box_center,box_size,x_min,x_max,y_min,y_max,
            cellsize_early,cellsize_late,cellwallthickness_early,cellwallthickness_late,\
            merge_operation,merge_tol,precrackFlag,precrack_widths,boundaryFlag,\
            fiberlength,theta_min,theta_max,z_min,z_max,long_connector_ratio,\
            NURBS_degree,nctrlpt_per_beam,nconnector_t_precrack,nconnector_l_precrack,\
            nParticles,nbeamElem,\
            STLFlag,\
            startTime,placementTime,voronoiTime,RebuildvorTime,BeamTime,FileTime)

if __name__ == "__main__":
    test_hexalattice_tri()