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
import ringspy.MeshGenTools as rpgen
import shutil
from pathlib import Path

def test_read_existing_sites():
    # performance check only
    startTime = time.time()
    
    # ==================================================================
    # Input parameters
    geoName = 'test_cube'
    path = 'meshes'
    
    radial_growth_rule = 'test_sites.npy'   
    iter_max = 500 # increase this number to achieve a more regular geometry
    print_interval = 500 # interval for printing prgress info
    
    # Radial cell growth parameters
    # length unit: mm
    r_min = 0   # inner radius of generation domain
    r_max = 4   # outer radius of generation domain
    nrings = 4 # number of rings
    width_heart = 0.3*(r_max-r_min)/nrings # ring width for the innermost ring
    width_sparse = 0.7*(r_max-r_min)/nrings # ring width for rings with sparse cells
    width_dense = 0.3*(r_max-r_min)/nrings # ring width for rings with dense cells
    generation_center = (0,0) # coordinates of generation domain center
    
    cellsize_sparse = 0.02
    cellsize_dense = 0.01
    cellwallthickness_sparse = 0.010
    cellwallthickness_dense = 0.006
    
    # clipping box parameters
    boundaryFlag = 'on'
    box_shape = 'square'
    box_center = (1.25,0) # coordinates of clipping box center
    box_size = 1.5 # side length
        
    # longitudinal direction parameters
    segment_length = 0.5*box_size
    theta_min = 0 # unit: radian
    theta_max = 0 # unit: radian
    z_min = 0
    z_max = box_size
    long_connector_ratio = 0.02 # longitudinal joint length = ratio * segment_length
    
    # material parameters
    skeleton_density = 1.5e-9 # unit: tonne/mm3
    
    # generation parameters
    merge_operation = 'on'
    merge_tol = 0.01
    
    precrackFlag = 'off'
    precrack_widths = 0.1
    
    stlFlag = 'on'
    
    inpFlag = 'on'
    inpType = 'Abaqus'
    
    # ==================================================================
    # Remove directory/files if exists already
    try:
        shutil.rmtree(Path(path, geoName))
    except:
        pass
    
    if not os.path.exists(Path(path, geoName)):
        os.makedirs(Path(path, geoName))
        
    # ==================================================================
    # Place cells with a specific radial growth pattern
    
    if radial_growth_rule == 'binary':
        # ---------------------------------------------
        # binary radial growth pattern (e.g. wood microstructure with earlywood-latewood alternations)
        sites,radii = rpgen.CellPlacement_Binary(generation_center,r_max,r_min,nrings,width_heart,\
                            width_sparse,width_dense,cellsize_sparse,cellsize_dense,\
                            iter_max,print_interval)
    elif radial_growth_rule == 'binary_lloyd':
        # ---------------------------------------------
        # binary with Lloyd's algorithm (e.g. wood microstructure with earlywood-latewood alternations, but more regular cell shapes)
        sites,radii = rpgen.CellPlacement_Binary_Lloyd(geoName,path,generation_center,r_max,r_min,\
                                                        nrings,width_heart,width_sparse,width_dense,\
                                                        cellsize_sparse,cellsize_dense,iter_max,\
                                                        print_interval)
            
    elif radial_growth_rule == 'regular_hexagonal':
        # ----------------------------------
        # hexagonal honeycomb-like geometry
        sites,radii = rpgen.CellPlacement_Honeycomb(generation_center,r_max,r_min,nrings,\
                            box_center,box_size,width_heart,\
                            width_sparse,width_dense,\
                            cellsize_sparse,cellsize_dense,\
                            cellwallthickness_sparse,cellwallthickness_dense,\
                            iter_max,print_interval)
    elif os.path.splitext(radial_growth_rule)[1] == '.npy':
        # ----------------------------------
        # load saved cell sites and radii data
        sites_path = Path(os.path.dirname(os.path.abspath(__file__)))
        sites,radii = rpgen.ReadSavedSites(sites_path,radial_growth_rule)
    
    else:
        print('Growth rule: {:s} is not supported for the current version, please check the README for more details.'.format(radial_growth_rule))
        print('Now exitting...')
        exit()
        
    placementTime = time.time()
    nParticles = len(sites)
    
    print('{:d} particles/cells placed in {:.3f} seconds'.format(nParticles, (placementTime - startTime)))
    
    # ==================================================================
    # Clipping box (boundaries) of the model
    x_min,x_max,y_min,y_max,boundaries,boundary_points,boundarylines = \
        rpgen.Clipping_Box(box_shape,box_center,box_size,boundaryFlag)
    
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
    
    # Visualize the original 2D Voronoi diagram
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
            rpgen.RebuildVoronoi_merge(vor,sites,boundaries,generation_center,x_min,x_max,y_min,y_max,box_center,box_shape,merge_tol,boundaryFlag)
        
        RebuildvorTime = time.time() 
        print('Voronoi tessellation rebuilt and merged in {:.3f} seconds'.format(RebuildvorTime - voronoiTime))
    else:
        [voronoi_vertices,finite_ridges,boundary_points,finite_ridges_new,\
          boundary_ridges_new,nvertex,nvertices_in,nfinite_ridge,nboundary_ridge,\
          nboundary_pts,nboundary_pts_featured,voronoi_ridges,nridge] = \
            rpgen.RebuildVoronoi(vor,sites,boundaries,generation_center,x_min,x_max,y_min,y_max,box_center,box_shape,boundaryFlag)
        
        RebuildvorTime = time.time() 
        print('Voronoi tessellation rebuilt in {:.3f} seconds'.format(RebuildvorTime - voronoiTime))
    
    # ===============================================
    # Insert mid and quarter points on the Voronoi ridges (can be used as potential failure positions on cell walls)
    [all_pts_2D,all_ridges,npt_per_layer,npt_per_layer_normal,npt_per_layer_vtk] = \
        rpgen.RidgeMidQuarterPts(voronoi_vertices,nvertex,nvertices_in,voronoi_ridges,\
                           finite_ridges_new,boundary_ridges_new,nfinite_ridge,\
                           nboundary_ridge,nboundary_pts,nboundary_pts_featured)
    
    # ==================================================================        
    # Generate a file for vertices and ridges info
    [all_vertices_2D, max_wings, flattened_all_vertices_2D, all_ridges] = \
        rpgen.VertexandRidgeinfo(all_pts_2D,all_ridges,\
                           npt_per_layer,npt_per_layer_normal,\
                           npt_per_layer_vtk,nridge,geoName,radii,generation_center,\
                           cellwallthickness_sparse,cellwallthickness_dense)
        
    ###############################################################################
    # Extrude in the parallel-to-grain (longitudinal) direction
    NURBS_degree = 2
    nctrlpt_per_beam = 5
    
    [IGAvertices,vertex_connectivity,beam_connectivity_original,nbeam_total,\
     beam_connectivity,nbeamElem,nlayer,connector_t_connectivity,\
     connector_t_bot_connectivity,connector_t_top_connectivity,\
     connector_t_reg_connectivity,connector_l_connectivity,nconnector_t_per_beam,\
     nconnector_t_per_grain,nconnector_t,nconnector_l,nconnector_total,\
     theta,z_coord,nbeam_per_grain,connector_l_vertex_dict] = \
        rpgen.GenerateBeamElement(NURBS_degree,nctrlpt_per_beam,\
                            segment_length,theta_min,theta_max,z_min,z_max,\
                            long_connector_ratio,npt_per_layer,voronoi_vertices,\
                            nvertex,voronoi_ridges,nridge,generation_center,\
                            all_vertices_2D,max_wings,flattened_all_vertices_2D,all_ridges)
    
    BeamTime = time.time() 
    print('{:d} beam elements generated in {:.3f} seconds'.format(nbeamElem, (BeamTime - RebuildvorTime)))
    
    # ==================================================================
    # Insert precracks
    if precrackFlag in ['on','On','Y','y','Yes','yes']:
        precrack_nodes = np.array([[x_indent, y_precrack, x_precrack, y_precrack]])
        [precrack_elem,nconnector_t_precrack,nconnector_l_precrack] = \
            rpgen.insert_precracks(all_pts_2D,all_ridges,nridge,npt_per_layer,\
                                     npt_per_layer_normal,npt_per_layer_vtk,\
                                     nlayer,precrack_nodes,precrack_widths,\
                                     cellsize_sparse)
    else:
        precrack_nodes = []
        precrack_elem = []
        nconnector_t_precrack = 0
        nconnector_l_precrack = 0
    
    # ==================================================================
    # Calculate mesh info
    height_connector_t = segment_length/4
    
    ConnMeshData = rpgen.ConnectorMeshFile(geoName,IGAvertices,connector_t_bot_connectivity,\
                     connector_t_reg_connectivity,connector_t_top_connectivity,\
                     height_connector_t,connector_l_connectivity,all_vertices_2D,\
                     max_wings,flattened_all_vertices_2D,nbeam_per_grain,nridge,\
                     connector_l_vertex_dict)

    # ==================================================================
    # Calculate model properties
    [mass,bulk_volume,bulk_density,porosity] = \
        rpgen.ModelInfo(box_shape,boundary_points,z_min,z_max,skeleton_density,ConnMeshData)
    
    # ==================================================================
    # Bezier extraction 
    knotVec = rpgen.BezierExtraction(NURBS_degree,nctrlpt_per_beam,nbeam_total)
    npatch = beam_connectivity_original.shape[0]
    
    mkBezierBeamFile = rpgen.BezierBeamFile(geoName,NURBS_degree,nctrlpt_per_beam,\
                       nconnector_t_per_beam,npatch,knotVec)
    
    # ==================================================================
    # Generate Paraview visulization files
    rpgen.VisualizationFiles(geoName,NURBS_degree,nlayer,npt_per_layer_vtk,all_pts_2D,\
                       segment_length,theta,z_coord,nbeam_per_grain,nridge,\
                       voronoi_ridges,generation_center,all_ridges,nvertex,nconnector_t,\
                       nconnector_l,nctrlpt_per_beam,ConnMeshData,all_vertices_2D,\
                       max_wings,flattened_all_vertices_2D)
        
    plt.savefig(Path('meshes/' + geoName + '/' + geoName + '.png'))
    
    # ==================================================================
    # Generate 3D model files
    if stlFlag in ['on','On','Y','y','Yes','yes']:
        rpgen.StlModelFile(geoName)
    
    # ==================================================================
    # Generate input files for numerical simulations
    if inpFlag in ['on','On','Y','y','Yes','yes']:
        if inpType in ['abaqus','Abaqus','ABQ','abq','ABAQUS','Abq']:
            rpgen.AbaqusFile(geoName,NURBS_degree,npatch,nbeam_per_grain,IGAvertices,beam_connectivity,\
                            connector_t_bot_connectivity,connector_t_reg_connectivity,\
                            connector_t_top_connectivity,connector_l_connectivity,nbeamElem,\
                            nconnector_t,nconnector_l,nconnector_t_precrack,nconnector_l_precrack,\
                            segment_length,height_connector_t,long_connector_ratio,\
                            x_max,x_min,y_max,y_min,z_coord,box_shape,box_size,\
                            cellwallthickness_sparse,cellwallthickness_dense,\
                            merge_operation,merge_tol,\
                            precrackFlag,precrack_elem)
        else:
            np.save(Path('meshes/' + geoName + '/' + geoName + '_sites.npy'),sites)
            np.save(Path('meshes/' + geoName + '/' + geoName + '_radii.npy'),radii)
            print('Input files type: {:s} is not supported for the current version, please check the README for more details.'.format(inpType))
            print('Generated cells and rings info has been saved.')
            print('Now exitting...')
            exit()
    
    FileTime = time.time() 
    print('Files generated in {:.3f} seconds'.format(FileTime - BeamTime))
    
    # ==================================================================
    # Generate log file for the generation
    rpgen.LogFile(geoName,iter_max,r_min,r_max,nrings,width_heart,width_sparse,width_dense,\
            generation_center,box_shape,box_center,box_size,x_min,x_max,y_min,y_max,
            cellsize_sparse,cellsize_dense,cellwallthickness_sparse,cellwallthickness_dense,\
            merge_operation,merge_tol,precrackFlag,precrack_widths,boundaryFlag,\
            segment_length,theta_min,theta_max,z_min,z_max,long_connector_ratio,\
            NURBS_degree,nctrlpt_per_beam,nconnector_t_precrack,nconnector_l_precrack,\
            nParticles,nbeamElem,skeleton_density,mass,bulk_volume,bulk_density,porosity,\
            stlFlag,inpFlag,inpType,radial_growth_rule,\
            startTime,placementTime,voronoiTime,RebuildvorTime,BeamTime,FileTime)

if __name__ == "__main__":
    test_read_existing_sites()