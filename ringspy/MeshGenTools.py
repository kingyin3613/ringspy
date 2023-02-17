# -*- coding: utf-8 -*-
"""
RingsPy utilities and functions

author: king_yin3613
email: haoyin2022@u.northwestern.edu
"""

import math
import bisect
import numpy as np
from scipy.spatial.distance import cdist
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from scipy.spatial import cKDTree
from scipy.spatial import ConvexHull
from scipy.spatial import Voronoi
from pathlib import Path
import datetime
import pkg_resources


def intersect2D(a, b):
  """
  Find row intersection between 2D numpy arrays, a and b.
  Returns another numpy array with shared rows
  credit to: https://gist.github.com/Robaina/b742f44f489a07cd26b49222f6063ef7
  """
  return np.array([x for x in set(tuple(x) for x in a) & set(tuple(x) for x in b)])

def TriangleArea2D(x1, y1, x2, y2, x3, y3):
    """
    A utility function to calculate area
    of triangle formed by (x1, y1),
    (x2, y2) and (x3, y3)
    """
    return abs((x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2)) / 2.0)

def check_isinside_boundcircle(circC,generation_center,rmin,rmax):
    """ Make sure this circle does not protrude outside the ring's bounding circle """
    x = circC[0]
    y = circC[1]
    w = circC[2]
    within = ( math.sqrt((x-generation_center[0])**2 + (y-generation_center[1])**2) > rmin+w ) and ( math.sqrt((x-generation_center[0])**2 + (y-generation_center[1])**2)<rmax-w )
    
    return within

def check_isinside_boundTriangle2D(circC,boundaries):
    """ 
    A function to check whether point P(x, y)
    lies inside the triangle formed by
    A(x1, y1), B(x2, y2) and C(x3, y3)
    credit to: https://www.geeksforgeeks.org/check-whether-a-given-point-lies-inside-a-triangle-or-not/
    """
    x = circC[0]
    y = circC[1]
    w = circC[2]
    
    x1 = boundaries[0][1][0][0]
    y1 = boundaries[0][1][0][1]
    x2 = boundaries[1][1][0][0]
    y2 = boundaries[1][1][0][1]
    x3 = boundaries[2][1][0][0]
    y3 = boundaries[2][1][0][1]
    
    # Calculate area of triangle ABC
    A = TriangleArea2D(x1, y1, x2, y2, x3, y3)
 
    # Calculate area of triangle PBC
    A1 = TriangleArea2D(x, y, x2, y2, x3, y3)
     
    # Calculate area of triangle PAC
    A2 = TriangleArea2D(x1, y1, x, y, x3, y3)
     
    # Calculate area of triangle PAB
    A3 = TriangleArea2D(x1, y1, x2, y2, x, y)
     
    # Check if sum of A1, A2 and A3 is same as A
    # within = (A == A1 + A2 + A3)
    # with tolerance
    within = (np.abs(A - (A1 + A2 + A3)) <= A*1e-3)

    return within

def check_isinside_boundbox2D(circC,x_min,x_max,y_min,y_max):
    """ Make sure this circle does not protrude outside the bounding box """
    x = circC[0]
    y = circC[1]
    w = circC[2]
    within = ( x > x_min+w ) and ( x < x_max-w ) and ( y > y_min+w ) and ( y < y_max-w )
    
    return within

def check_isinside_boundbox2Dindent(circC,x_min,x_max,y_min,y_max,x_indent,y_indent_min,y_indent_max):
    """ Make sure this circle does not protrude outside the indent bounding box """
    x = circC[0]
    y = circC[1]
    w = circC[2]
    within = ( x > x_min+w ) and ( x < x_max-w ) and ( y > y_min+w ) and ( y < y_max-w )
    within_indent = ( x > x_min+w ) and ( x < x_indent-w ) and ( y > y_indent_min+w ) and ( y < y_indent_max-w )
    within = np.subtract(within,within_indent,dtype=np.float32)
    return within

def check_isinside_boundbox2Dprecrack(circC,x_min,x_max,y_min,y_max,x_indent,y_indent,x_precrack,y_precrack):
    """ Make sure this circle does not protrude outside the ring's bounding box """
# Neither/nor set operation, ref:https://math.stackexchange.com/questions/1257097/how-to-represent-a-neither-nor-set-operation
    x = circC[0]
    y = circC[1]
    w = circC[2]
    within = ( x > x_min+w ) and ( x < x_max-w ) and ( y > y_min+w ) and ( y < y_max-w )
    within_indent = ( x > x_min+w ) and ( x < x_indent-w ) and ( y > -y_indent+w ) and ( y < y_indent-w )
    within_precrack = ( x > x_indent+w ) and ( x < x_precrack-w ) and ( y > -y_precrack+w ) and ( y < y_precrack-w )
    within_both = within_indent or within_precrack
    within = np.subtract(within,within_both,dtype=np.float32)
    return within

def check_isinside_boundHexagon2D(circC,L,box_center):
    """ Make sure this circle does not protrude outside the bounding hexagon 
    credit to: http://www.playchilla.com/how-to-check-if-a-point-is-inside-a-hexagon
    """
    x = np.abs(circC[0] - box_center[0])
    y = np.abs(circC[1] - box_center[1])
    v = L/2
    h = L/np.sqrt(3)/2
    within = ( x < 2*h ) and ( y < v ) and ((2*v*h - x*v - h*y) > 0)
    return within

def check_overlap(circles,circC):
    """ Make sure the distance between the current circle's center and all
        other circle centers is greater than or equal to the circle's perimeter (2r)
    """
    x = circC[0]
    y = circC[1]
    w = circC[2]
    nooverlap = all( (x - c[0])**2 + (y - c[1])**2 >= (w + c[2])**2 for c in circles )
    return nooverlap

def check_iscollinear(p1,p2,boundaries):
    """ Make sure the line segment x1-x2 is collinear with 
        any line segment of the boundary polygon
    """
    collinear = 0
    if (p2[0]-p1[0]) == 0: # inf slope
        k1 = 99999999 
    else:
        k1 = (p2[1]-p1[1])/(p2[0]-p1[0]) # slope of the new line segment
        
    for boundary in boundaries: # loop over boundary lines
        bp1 = boundary[1][0]
        bp2 = boundary[1][1]
        
        if (bp2[0]-bp1[0]) == 0: # inf slope
            k2 = 99999999
        else:
            k2 = (bp2[1]-bp1[1])/(bp2[0]-bp1[0]) # slope of the boundary line segment
            
        if math.isclose(k1,k2): # check if slopes are equal (with a tolerance)
            p3 = (p1+p2)/2 # mid point of the new line
            if k2 == 99999999: # inf slope
                p3_on = (p3[0] == bp1[0])
            else:
                p3_on = math.isclose((p3[1] - bp1[1]),k2*(p3[0] - bp1[0]))
            p3_between = (min(bp1[0], bp2[0]) <= p3[0] <= max(bp1[0], bp2[0])) and (min(bp1[1], bp2[1]) <= p3[1] <= max(bp1[1], bp2[1]))
            if (p3_on and p3_between): # check if mid point of the new line is on the boundary line segment
                collinear += 1
                
    return collinear

def rotate_around_point_highperf(xy, radians, origin=(0, 0)):
    """Rotate a point around a given point.
    
    I call this the "high performance" version since we're caching some
    values that are needed >1 time. It's less readable than the previous
    function but it's faster.
    """
    x, y = xy
    offset_x, offset_y = origin
    adjusted_x = (x - offset_x)
    adjusted_y = (y - offset_y)
    cos_rad = math.cos(radians)
    sin_rad = math.sin(radians)
    qx = offset_x + cos_rad * adjusted_x + sin_rad * adjusted_y
    qy = offset_y + -sin_rad * adjusted_x + cos_rad * adjusted_y

    return qx, qy


def Clipping_Box(box_shape,box_center,box_size,boundaryFlag):
    """
    clipping box for the 2D Voronoi diagram
    """
    
    if box_shape in ['square','Square','SQUARE', 'cube','s']: # sqaure box
        x_min = box_center[0] - box_size/2
        x_max = box_center[0] + box_size/2 
        y_min = box_center[1] - box_size/2
        y_max = box_center[1] + box_size/2
        
        boundary_points = np.array([[x_min, y_min], [x_max, y_min], [x_max, y_max], [x_min, y_max]])
        l0 = [(x_min, y_min), (x_max, y_min)]
        l1 = [(x_max, y_min), (x_max, y_max)]
        l2 = [(x_max, y_max), (x_min, y_max)]
        l3 = [(x_min, y_max), (x_min, y_min)]
        boundaries = [('bottom', l0), ('right', l1), ('top', l2), ('left', l3)]

        if boundaryFlag in ['on','On','Y','y','Yes','yes']:
            boundarylines = patches.Rectangle(boundaries[0][1][0], x_max-x_min, y_max-y_min, linewidth=2, edgecolor='k',facecolor='none')
        else:
            boundarylines = patches.Rectangle(boundaries[0][1][0], x_max-x_min, y_max-y_min, linewidth=2, linestyle='-.', edgecolor='k',facecolor='none')
            
    elif box_shape in ['triangle','Triangle','TRIANGLE','triangular','tri','^']: # Triangle box
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
        
    elif box_shape in ['hexagon','Hexagon','HEXAGON','hexagonal','hex','h']: # hexagonal box
        x_0 = box_center[1] - box_size
        x_1 = box_center[1] - box_size/2.0
        x_2 = box_center[1] + box_size/2.0
        x_3 = box_center[1] + box_size
        y_0 = box_center[0] - box_size/2*np.sqrt(3)
        y_1 = box_center[0]
        y_2 = box_center[0] + box_size/2*np.sqrt(3)
        
        boundary_points = np.array([[x_1, y_0], [x_2, y_0], [x_3, y_1], [x_2, y_2], [x_1, y_2], [x_0, y_1]])
        
        l0 = [(x_1, y_0), (x_2, y_0)]
        l1 = [(x_2, y_0), (x_3, y_1)]
        l2 = [(x_3, y_1), (x_2, y_2)]
        l3 = [(x_2, y_2), (x_1, y_2)]
        l4 = [(x_1, y_2), (x_0, y_1)]
        l5 = [(x_0, y_1), (x_1, y_0)]
        
        boundaries = [('bottom', l0), ('bottom-right', l1), ('top-right', l2), ('top', l3), ('top-left', l4), ('bottom-left', l5)]
        
        if boundaryFlag in ['on','On','Y','y','Yes','yes']:
            boundarylines = patches.RegularPolygon(box_center,numVertices=6,radius=box_size,orientation=np.pi/6,linewidth=2,edgecolor='k',facecolor='none')
        else:
            boundarylines = patches.RegularPolygon(box_center,numVertices=6,radius=box_size,orientation=np.pi/6,linewidth=2,linestyle='-.',edgecolor='k',facecolor='none')
        
        x_min = x_0
        x_max = x_3
        y_min = y_0
        y_max = y_2
        
    else:
        print('box_shape: {:s} is not supported for the current version, please check the README for more details.'.format(box_shape))
        print('Now exitting...')
        exit()
        
    return x_min,x_max,y_min,y_max,boundaries,boundary_points,boundarylines


def relax_points(vor,omega):
    """
    Moves each point to the centroid of its cell in the Voronoi map to "relax"
    the points (i.e. jitter them so as to spread them out within the space).
    """

    filtered_regions = []
    
    nonempty_regions = list(filter(None, vor.regions))
    for region in nonempty_regions:
        if not any(x < 0 for x in region):
            filtered_regions.append(region)
    
    centroids = np.zeros((len(filtered_regions),2))
    for i in range(0,len(filtered_regions)):   
        vertices = vor.vertices[filtered_regions[i] + [filtered_regions[i][0]], :]
        # plt.plot(vertices[:,0],vertices[:,1], 'go', markersize=2.0)
        centroid = find_centroid(vertices,omega) # get the centroid of these verts
        
        # plt.plot(centroid[:,0],centroid[:,1], 'ro', markersize=2.0)
        centroids[i,:] = centroid
        
    return centroids # store the centroids as the new point positions


def find_centroid(vertices,omega):
    """
    Find the centroid of a Voroni region described by `vertices`, and return a
    np array with the x and y coords of that centroid.
    The equation for the method used here to find the centroid of a 2D polygon
    is given here: https://en.wikipedia.org/wiki/Centroid#Of_a_polygon
    @params: np.array `vertices` a numpy array with shape n,2
    @returns np.array a numpy array that defines the x, y coords
      of the centroid described by `vertices`
    """
    area = 0
    centroid_x = 0
    centroid_y = 0
    for i in range(len(vertices)-1):
      step = (vertices[i, 0] * vertices[i+1, 1]) - (vertices[i+1, 0] * vertices[i, 1])
      
      step = step*omega # omega -> relaxation factor (>1 for over-relaxation, faster convergence)
      
      area += step
      centroid_x += (vertices[i, 0] + vertices[i+1, 0]) * step
      centroid_y += (vertices[i, 1] + vertices[i+1, 1]) * step
    area /= 2
    centroid_x = (1.0/(6.0*area)) * centroid_x
    centroid_y = (1.0/(6.0*area)) * centroid_y
    return np.array([[centroid_x, centroid_y]])


def CellPlacement_Binary(generation_center,r_max,r_min,nrings,width_heart,
                         width_sparse,width_dense,cellsize_sparse,cellsize_dense,\
                         iter_max,print_interval):
    """
    packing cells by randomly placing new cells in the generation rings
    """
    
    # generate radii for rings
    width = np.concatenate(([width_heart],np.tile([width_sparse,width_dense],nrings)))
    noise = np.random.normal(1,0.25,len(width))
    width = np.multiply(width,noise)
    radii = np.concatenate(([0],np.cumsum(width)))
    
    # Place circular cells in each ring (can be parallelized in the future)
    # circles: list = list()
    circles = []
    for ibin in range(0,nrings*2+1):
        iter = 0 
        while iter < iter_max:
            r = (radii[ibin+1]-radii[ibin])*math.sqrt(np.random.random()) + radii[ibin]  # randomly generate a point in the domain
            t = np.random.random()*2*math.pi 
            
            x = r*math.cos(t) + generation_center[0]
            y = r*math.sin(t) + generation_center[1]
    
            if (ibin % 2) == 0: # if even, dense cells
                w = cellsize_dense
            else:
                w = cellsize_sparse
    
            circC = [x, y, w]
            
            if check_overlap(circles, circC) and check_isinside_boundcircle(circC,generation_center,radii[ibin],radii[ibin+1]):
                circles.append(circC)
                
                # regularly print 
                if (len(circles) % print_interval == 0):
                    if (nrings*2+1-ibin == 1):
                        print('{:d} ring remaining; {:d} cells/particles placed.'.format(nrings*2+1-ibin,len(circles)))
                    else:
                        print('{:d} rings remaining; {:d} cells/particles placed.'.format(nrings*2+1-ibin,len(circles)))
            iter += 1 
        
    sites = np.array(circles)

    # out-of-boundary detection for Voronoi vertices
    outofbound = []
    for i in range(0,sites.shape[0]):
        if( ((sites[i,0]-generation_center[0])**2+(sites[i,1]-generation_center[1])**2) > (1.2*r_max)**2 ):
            outofbound.append(i)
    sites = np.delete(sites,outofbound,0)

    return sites, radii


def CellPlacement_Binary_Lloyd(geoName,path,generation_center,r_max,r_min,\
                               nrings,width_heart,width_sparse,width_dense,\
                               cellsize_sparse,cellsize_dense,iter_max,\
                               print_interval,omega=1.0):
    """
    packing cells by firstly placing new cells and then performing Lloyd's relaxation in the generation rings
    """
    # generate radii for rings
    radii = np.concatenate(([width_heart],np.tile([width_sparse,width_dense],nrings)))
    noise = np.random.normal(1,0.25,len(radii))
    radii = np.multiply(radii,noise)
    radii = np.concatenate(([0],np.cumsum(radii)))

    # generate perimeter points for heart region
    PerimeterPoints = []
    npoint = int(np.ceil(2*np.pi*radii[1]/(2*cellsize_dense)))
    t = np.linspace(0, 2*np.pi, npoint, endpoint=False)
    x = radii[1] * np.cos(t)
    y = radii[1] * np.sin(t)
    w = cellsize_dense*np.ones(x.shape)
    PerimeterPoints.append(np.c_[x, y])
    PerimeterPointsSites = np.concatenate(PerimeterPoints, axis=0)
    
    #generate internal points for heart region
    n_nonoverlapped_cells = int(2*np.floor((radii[1]/cellsize_dense)**2))
    inside_cells = 1e-3*(np.random.rand(n_nonoverlapped_cells,2) - [0.5,0.5]) # add cell points in a very small area
    
    sites = np.vstack((PerimeterPointsSites,inside_cells))
    
    filenames = []
    frame_iter = 0
    
    lloyd_iter = 0
    while lloyd_iter < iter_max:
        vor = Voronoi(sites)
        sites = relax_points(vor,omega)
        
        sites = np.vstack((PerimeterPointsSites,sites))
        lloyd_iter += 1

    
    existing_sites = np.copy(sites[npoint:,:])
    
    # generate sites for each ring    
    for i in range(1,len(radii)-1):
        
        OuterPerimeterPoints = []
        OuterPerimeter_radii = radii[i+1]
        
        if (i % 2) == 0: # if even, dense cells
            cellsize = cellsize_dense
        else:
            cellsize = cellsize_sparse

        OuterPerimeter_npoint = int(np.ceil(2*np.pi*OuterPerimeter_radii/(2*cellsize)))
        
        t = np.linspace(0, 2*np.pi, OuterPerimeter_npoint, endpoint=False)
        x = OuterPerimeter_radii * np.cos(t)
        y = OuterPerimeter_radii * np.sin(t)
        w = cellsize*np.ones(x.shape)
        OuterPerimeterPoints.append(np.c_[x, y])
        OuterPerimeterPointsSites = np.concatenate(OuterPerimeterPoints, axis=0)
        
        # generate internal points
        nsubrings = int(0.5*np.ceil((radii[i+1]-radii[i])/cellsize))
        subradii = np.linspace(radii[i], radii[i+1], nsubrings, endpoint=False)
        
        subnpoints =np.ceil(2*np.pi*subradii/(2*cellsize)).astype(int)
        
        circles = []
        
        for subr, subnpoint in zip(subradii, subnpoints):
            t = np.linspace(0, 2*np.pi, subnpoint, endpoint=False)
            x = subr * np.cos(t)
            y = subr * np.sin(t)
            circles.append(np.c_[x, y])
            
        inside_cells = np.concatenate(circles, axis=0)
        
        noise = (np.random.rand(len(inside_cells),2) - [0.5,0.5])*cellsize
        inside_cells = inside_cells + 0.25*noise

        sites = np.vstack((OuterPerimeterPointsSites,inside_cells))
                    
        existing_sites = np.vstack((sites,existing_sites))
        PerimeterPointsSites = np.copy(OuterPerimeterPointsSites)
    sites = existing_sites

    return sites, radii


def CellPlacement_Honeycomb(generation_center,r_max,r_min,nrings,box_center,box_size,\
                            width_heart,width_sparse,width_dense,\
                            cellsize_sparse,cellsize_dense,\
                            cellwallthickness_sparse,cellwallthickness_dense,\
                            iter_max,print_interval):
    from hexalattice.hexalattice import create_hex_grid
    """
    packing cells in hexagonal grids
    """
    
    # generate radii for rings
    width = np.concatenate(([width_heart],np.tile([width_sparse,width_dense],nrings)))
    noise = np.random.normal(1,0.25,len(width))
    width = np.multiply(width,noise)
    radii = np.concatenate(([0],np.cumsum(width)))
    
    cellangle = 0.0
    cellsize = (cellsize_sparse+cellsize_dense)/2
    nx = int(2*r_max/cellsize)
    ny = int(2*r_max/cellsize)
    # The hexagonal lattice sites are generated via hexalattice: https://pypi.org/project/hexalattice/
    hex_centers, _ = create_hex_grid(nx=nx,
                                     ny=ny,
                                     min_diam=cellsize,
                                     rotate_deg=cellangle,
                                     do_plot=False)
    sites = np.hstack((hex_centers+generation_center,np.ones([hex_centers.shape[0],1])))
    sites = np.asarray(sites)

    # out-of-boundary detection for Voronoi vertices
    outofbound = []
    for i in range(0,sites.shape[0]):
        if( ((sites[i,0]-generation_center[0])**2+(sites[i,1]-generation_center[1])**2) > (1.2*r_max)**2 ):
            outofbound.append(i)
    sites = np.delete(sites,outofbound,0)

    return sites, radii


def RebuildVoronoi(vor,circles,boundaries,generation_center,x_min,x_max,y_min,y_max,box_center,box_shape,boundaryFlag):
    """Clip Voronoi mesh by the boundaries, rebuild the new Voronoi mesh"""
    # Store indices of Voronoi vertices for each finite ridge
    finite_ridges = []
    finite_ridges_pointid = []
    infinite_ridges = []
    infinite_ridges_pointid = []
    boundary_points = []
    
    # Form two new Voronoi vertices lists with and without out-of-bound vertices
    voronoi_vertices_in = []
    voronoi_vertices_out = []
    
    if box_shape == 'triangle':
        for pointidx, simplex in zip(vor.ridge_points, vor.ridge_vertices):
            simplex = np.asarray(simplex)
            simplex = np.where(simplex == -1, -9999999999, simplex) # use a large negative number instead of -1 to represent the infinite ridge vertices 
            if np.all(simplex >= 0) and check_isinside_boundTriangle2D(np.append(vor.vertices[simplex[0]],0),boundaries) and check_isinside_boundTriangle2D(np.append(vor.vertices[simplex[1]],0),boundaries):
                voronoi_vertices_in.append(vor.vertices[simplex[0]])
                voronoi_vertices_in.append(vor.vertices[simplex[1]])
                finite_ridges.append(simplex)
                finite_ridges_pointid.append(pointidx)
                plt.plot(vor.vertices[simplex, 0], vor.vertices[simplex, 1], 'k-',linewidth=1)
            else: # infinite Voronoi ridges and finite ridges with out-of-bound vetices
                if np.all(simplex >= 0) and (check_isinside_boundTriangle2D(np.append(vor.vertices[simplex[0]],0),boundaries) or check_isinside_boundTriangle2D(np.append(vor.vertices[simplex[1]],0),boundaries)): # finite Voronoi ridge but one vertex is far
        
                    if (check_isinside_boundTriangle2D(np.append(vor.vertices[simplex[0]],0),boundaries) == 0):# if vertex is out-of-bound, inverse its index
                        voronoi_vertices_out.append(vor.vertices[simplex[0]])
                        if simplex[0] == 0:
                            simplex[0] = -9999999998 # since inverse of 0 is still zero, we use another large number to represent -0
                        else:
                            simplex[0] = -simplex[0]
                    elif (check_isinside_boundTriangle2D(np.append(vor.vertices[simplex[1]],0),boundaries) == 0):
                        voronoi_vertices_out.append(vor.vertices[simplex[1]])
                        if simplex[1] == 0:
                            simplex[1] = -9999999998 # since inverse of 0 is still zero, we use another large number to represent -0
                        else:
                            simplex[1] = -simplex[1]
                    infinite_ridges.append(simplex)
                    infinite_ridges_pointid.append(pointidx) 
                elif check_isinside_boundTriangle2D(np.append(vor.vertices[simplex[simplex >= 0][0]],0),boundaries): # index of finite end Voronoi vertex of the infinite ridge
                    infinite_ridges.append(simplex)
                    infinite_ridges_pointid.append(pointidx) 
    elif box_shape == 'square':
        for pointidx, simplex in zip(vor.ridge_points, vor.ridge_vertices):
            simplex = np.asarray(simplex)
            simplex = np.where(simplex == -1, -9999999999, simplex) # use a large negative number instead of -1 to represent the infinite ridge vertices 
            if np.all(simplex >= 0) and check_isinside_boundbox2D(np.append(vor.vertices[simplex[0]],0),x_min,x_max,y_min,y_max) and check_isinside_boundbox2D(np.append(vor.vertices[simplex[1]],0),x_min,x_max,y_min,y_max):
                voronoi_vertices_in.append(vor.vertices[simplex[0]])
                voronoi_vertices_in.append(vor.vertices[simplex[1]])
                finite_ridges.append(simplex)
                finite_ridges_pointid.append(pointidx)
                plt.plot(vor.vertices[simplex, 0], vor.vertices[simplex, 1], 'k-',linewidth=1)
            else: # infinite Voronoi ridges and finite ridges with out-of-bound vetices
                if np.all(simplex >= 0) and (check_isinside_boundbox2D(np.append(vor.vertices[simplex[0]],0),x_min,x_max,y_min,y_max) or check_isinside_boundbox2D(np.append(vor.vertices[simplex[1]],0),x_min,x_max,y_min,y_max)): # finite Voronoi ridge but one vertex is far
        
                    if (check_isinside_boundbox2D(np.append(vor.vertices[simplex[0]],0),x_min,x_max,y_min,y_max) == 0):# if vertex is out-of-bound, inverse its index
                        voronoi_vertices_out.append(vor.vertices[simplex[0]])
                        if simplex[0] == 0:
                            simplex[0] = -9999999998 # since inverse of 0 is still zero, we use another large number to represent -0
                        else:
                            simplex[0] = -simplex[0]
                    elif (check_isinside_boundbox2D(np.append(vor.vertices[simplex[1]],0),x_min,x_max,y_min,y_max) == 0):
                        voronoi_vertices_out.append(vor.vertices[simplex[1]])
                        if simplex[1] == 0:
                            simplex[1] = -9999999998 # since inverse of 0 is still zero, we use another large number to represent -0
                        else:
                            simplex[1] = -simplex[1]
                    infinite_ridges.append(simplex)
                    infinite_ridges_pointid.append(pointidx) 
                elif check_isinside_boundbox2D(np.append(vor.vertices[simplex[simplex >= 0][0]],0),x_min,x_max,y_min,y_max): # index of finite end Voronoi vertex of the infinite ridge
                    infinite_ridges.append(simplex)
                    infinite_ridges_pointid.append(pointidx) 
    elif box_shape == 'hexagon':
        L = y_max - y_min
        box_center = generation_center
        for pointidx, simplex in zip(vor.ridge_points, vor.ridge_vertices):
            simplex = np.asarray(simplex)
            simplex = np.where(simplex == -1, -9999999999, simplex) # use a large negative number instead of -1 to represent the infinite ridge vertices 
            if np.all(simplex >= 0) and check_isinside_boundHexagon2D(np.append(vor.vertices[simplex[0]],0),L,box_center) and check_isinside_boundHexagon2D(np.append(vor.vertices[simplex[1]],0),L,box_center):
                voronoi_vertices_in.append(vor.vertices[simplex[0]])
                voronoi_vertices_in.append(vor.vertices[simplex[1]])
                
                finite_ridges.append(simplex)
                finite_ridges_pointid.append(pointidx)
                plt.plot(vor.vertices[simplex, 0], vor.vertices[simplex, 1], 'k-',linewidth=1)
            else: # infinite Voronoi ridges and finite ridges with out-of-bound vetices
                
                if np.all(simplex >= 0) and (check_isinside_boundHexagon2D(np.append(vor.vertices[simplex[0]],0),L,box_center) or check_isinside_boundHexagon2D(np.append(vor.vertices[simplex[1]],0),L,box_center)): # finite Voronoi ridge but one vertex is far
        
                    if (check_isinside_boundHexagon2D(np.append(vor.vertices[simplex[0]],0),L,box_center) == 0):# if vertex is out-of-bound, inverse its index
                        voronoi_vertices_out.append(vor.vertices[simplex[0]])
                        if simplex[0] == 0:
                            simplex[0] = -9999999998 # since inverse of 0 is still zero, we use another large number to represent -0
                        else:
                            simplex[0] = -simplex[0]
                    elif (check_isinside_boundHexagon2D(np.append(vor.vertices[simplex[1]],0),L,box_center) == 0):
                        voronoi_vertices_out.append(vor.vertices[simplex[1]])
                        if simplex[1] == 0:
                            simplex[1] = -9999999998 # since inverse of 0 is still zero, we use another large number to represent -0
                        else:
                            simplex[1] = -simplex[1]
                    infinite_ridges.append(simplex)
                    infinite_ridges_pointid.append(pointidx) 
                elif check_isinside_boundHexagon2D(np.append(vor.vertices[simplex[simplex >= 0][0]],0),L,box_center): # index of finite end Voronoi vertex of the infinite ridge
                    infinite_ridges.append(simplex)
                    infinite_ridges_pointid.append(pointidx) 
    elif box_shape == 'notched_square': # notched
        x_min = boundaries[0][1][0][0]
        x_max = boundaries[1][1][0][0]
        y_min = boundaries[0][1][0][1]
        y_max = boundaries[2][1][0][1]
        x_indent = boundaries[5][1][0][0]
        y_indent_min = boundaries[6][1][0][1]
        y_indent_max = boundaries[4][1][0][1]
        
        
        for pointidx, simplex in zip(vor.ridge_points, vor.ridge_vertices):
            simplex = np.asarray(simplex)
            simplex = np.where(simplex == -1, -9999999999, simplex) # use a large negative number instead of -1 to represent the infinite ridge vertices 
            if np.all(simplex >= 0) and check_isinside_boundbox2Dindent(np.append(vor.vertices[simplex[0]],0),x_min,x_max,y_min,y_max,x_indent,y_indent_min,y_indent_max) and check_isinside_boundbox2Dindent(np.append(vor.vertices[simplex[1]],0),x_min,x_max,y_min,y_max,x_indent,y_indent_min,y_indent_max):
                voronoi_vertices_in.append(vor.vertices[simplex[0]])
                voronoi_vertices_in.append(vor.vertices[simplex[1]])
                finite_ridges.append(simplex)
                finite_ridges_pointid.append(pointidx)
                plt.plot(vor.vertices[simplex, 0], vor.vertices[simplex, 1], 'k-')
            else: # infinite Voronoi ridges and finite ridges with out-of-bound vetices
                if np.all(simplex >= 0) and (check_isinside_boundbox2Dindent(np.append(vor.vertices[simplex[0]],0),x_min,x_max,y_min,y_max,x_indent,y_indent_min,y_indent_max) or check_isinside_boundbox2Dindent(np.append(vor.vertices[simplex[1]],0),x_min,x_max,y_min,y_max,x_indent,y_indent_min,y_indent_max)): # finite Voronoi ridge but one vertex is far

                    if (check_isinside_boundbox2Dindent(np.append(vor.vertices[simplex[0]],0),x_min,x_max,y_min,y_max,x_indent,y_indent_min,y_indent_max) == 0):# if vertex is out-of-bound, inverse its index
                        voronoi_vertices_out.append(vor.vertices[simplex[0]])
                        if simplex[0] == 0:
                            simplex[0] = -9999999998 # since inverse of 0 is still zero, we use another large number to represent -0
                        else:
                            simplex[0] = -simplex[0]
                    elif (check_isinside_boundbox2Dindent(np.append(vor.vertices[simplex[1]],0),x_min,x_max,y_min,y_max,x_indent,y_indent_min,y_indent_max) == 0):
                        voronoi_vertices_out.append(vor.vertices[simplex[1]])
                        if simplex[1] == 0:
                            simplex[1] = -9999999998 # since inverse of 0 is still zero, we use another large number to represent -0
                        else:
                            simplex[1] = -simplex[1]
                    infinite_ridges.append(simplex)
                    infinite_ridges_pointid.append(pointidx) 
                elif check_isinside_boundbox2Dindent(np.append(vor.vertices[simplex[simplex >= 0][0]],0),x_min,x_max,y_min,y_max,x_indent,y_indent_min,y_indent_max): # index of finite end Voronoi vertex of the infinite ridge
                    infinite_ridges.append(simplex)
                    infinite_ridges_pointid.append(pointidx) 
    else:
        return
                    
    voronoi_vertices_in = np.unique(voronoi_vertices_in,axis=0)
    voronoi_vertices_out = np.unique(voronoi_vertices_out,axis=0)
    # plt.plot(voronoi_vertices_in[:, 0], voronoi_vertices_in[:, 1], 'ro')
    
    
            
    # Find the intersect point of infinite ridges (rays) with the boundary lines
    # Ref: https://stackoverflow.com/questions/14307158/how-do-you-check-for-intersection-between-a-line-segment-and-a-line-ray-emanatin#:~:text=Let%20r%20%3D%20(cos%20%CE%B8%2C,0%20%E2%89%A4%20u%20%E2%89%A4%201).&text=Then%20your%20line%20segment%20intersects,0%20%E2%89%A4%20u%20%E2%89%A4%201.
    for pointidx, simplex in zip(infinite_ridges_pointid, infinite_ridges):
        simplex = np.asarray(simplex)
            
        i = simplex[simplex >= 0][0] # index of finite end Voronoi vertex of the infinite ridge
        p = vor.vertices[i]
        s = circles[pointidx[1],0:2] - circles[pointidx[0],0:2]
        if simplex[simplex < 0][0] == -9999999999:
            tangent = s / np.linalg.norm(s)
            normal = np.array([-tangent[1], tangent[0]]) # normal
            midpoint = circles[:,0:2][pointidx].mean(axis=0) # midpoint between site points
            for boundary in boundaries: # loop over boundary lines
                q = np.asarray(boundary[1][0])
                s = np.asarray(boundary[1][1]) - np.asarray(boundary[1][0])
                t = np.cross((q-p),s)/np.cross(normal,s)
                u = np.cross((q-p),normal)/np.cross(normal,s)
                
                if (u >= 0) and (u <= 1):
                    if np.sign(np.dot(midpoint - generation_center, normal)) == 1: # facing outwards
                        if (t >= 0) and math.isfinite(t):
                            t_final = t
                    elif np.sign(np.dot(midpoint - generation_center, normal)) == -1: # facing inwards
                        if (t < 0) and math.isfinite(t):
                            t_final = t
                            
        elif simplex[simplex < 0][0] == -9999999998:
            t_final = np.inf
            normal = (vor.vertices[0] - p)/np.linalg.norm(vor.vertices[0] - p)
            for boundary in boundaries: # loop over boundary lines
                q = np.asarray(boundary[1][0])
                s = np.asarray(boundary[1][1]) - np.asarray(boundary[1][0])
                t = np.cross((q-p),s)/np.cross(normal,s)
                u = np.cross((q-p),normal)/np.cross(normal,s)
                
                if (u >= 0) and (u <= 1):
                    if (t >= 0) and (t < t_final):
                        t_final = t
        else:
            t_final = np.inf
            normal = (vor.vertices[-simplex[simplex < 0][0]] - p)/np.linalg.norm(vor.vertices[-simplex[simplex < 0][0]] - p)
            for boundary in boundaries: # loop over boundary lines
                q = np.asarray(boundary[1][0])
                s = np.asarray(boundary[1][1]) - np.asarray(boundary[1][0])
                t = np.cross((q-p),s)/np.cross(normal,s)
                u = np.cross((q-p),normal)/np.cross(normal,s)
                
                if (u >= 0) and (u <= 1):
                    if (t >= 0) and (t < t_final):
                        t_final = t
                            
        intersect_point = p + normal * t_final
        boundary_points.append(intersect_point)
        
        plt.plot([vor.vertices[i,0], intersect_point[0]], [vor.vertices[i,1], intersect_point[1]], 'k-',linewidth=1)
            
    if boundaryFlag in ['on','On','Y','y','Yes','yes']:
        boundary_points = np.asarray(boundary_points)
        boundary_points_featured = [x[1][0] for x in boundaries]
        boundary_points_featured = np.reshape(np.asarray(boundary_points_featured),(-1,2))
        nboundary_pts_featured = boundary_points_featured.shape[0]
        boundary_points = np.vstack((boundary_points,boundary_points_featured))
        nvertices_in = voronoi_vertices_in.shape[0]
        nboundary_pts = boundary_points.shape[0]
        finite_ridges = np.asarray(finite_ridges)
        nfinite_ridge = finite_ridges.shape[0]
        infinite_ridges = np.asarray(infinite_ridges)
        ninfinite_ridge = infinite_ridges.shape[0]
           
    else:
        boundary_points = np.asarray(boundary_points)
        nboundary_pts_featured = 0
        nvertices_in = voronoi_vertices_in.shape[0]
        nboundary_pts = boundary_points.shape[0]
        finite_ridges = np.asarray(finite_ridges)
        nfinite_ridge = finite_ridges.shape[0]
        infinite_ridges = np.asarray(infinite_ridges)
        ninfinite_ridge = infinite_ridges.shape[0]
           
    # reconstruct the connectivity for ridges since the unique operation rearrange the order of vertices (need to find a more efficient way like vectorize)
    finite_ridges_new = np.copy(finite_ridges)
    for i in range(0,nfinite_ridge): # loop over voronoi ridges in original finite_ridge list 
        ii = np.where(np.all(voronoi_vertices_in==vor.vertices[finite_ridges[i,0]],axis=1))[0][0]  # index of first vertex in voronoi_vertices_in array
        jj = np.where(np.all(voronoi_vertices_in==vor.vertices[finite_ridges[i,1]],axis=1))[0][0]  # index of second vertex in voronoi_vertices_in array
        finite_ridges_new[i,:] = [ii,jj]
            
    infinite_ridges_new = np.copy(infinite_ridges)
    for i in range(0,ninfinite_ridge): # loop over voronoi ridges in original infinite_ridge list
        if (infinite_ridges[i,:] < 0).argmax(axis=0) == 0: # first vertex is of negative index
            ii = nvertices_in + i
            jj = np.where(np.all(voronoi_vertices_in==vor.vertices[infinite_ridges[i,1]],axis=1))[0][0]  # index of finite vertex in voronoi_vertices_in array
        else:
            ii = np.where(np.all(voronoi_vertices_in==vor.vertices[infinite_ridges[i,0]],axis=1))[0][0]  # index of finite vertex in voronoi_vertices_in array
            jj = nvertices_in + i
        infinite_ridges_new[i,:] = [ii,jj]
    
    if boundaryFlag in ['on','On','Y','y','Yes','yes']:
        # construct the connectivity for a line path consisting of the boundary points
        boundary_ridges_new = np.zeros(boundary_points.shape)
        boundary_points_new = np.copy(boundary_points) 
        next_point = np.copy(boundary_points_new[0])  # get first point
        boundary_points_new[0] = [np.inf,np.inf]
        
        for i in range(0,boundary_points_new.shape[0]-1):
            next_point_id = cdist([next_point], boundary_points_new).argmin()
            if check_iscollinear(next_point,boundary_points_new[next_point_id],boundaries) > 0: # check if the new line segment is collinear with any boundary line segment
                boundary_ridges_new[i,1] = next_point_id
                boundary_ridges_new[i+1,0] = next_point_id
                next_point = np.copy(boundary_points_new[next_point_id])
                boundary_points_new[next_point_id] = [np.inf,np.inf]
            else:
                boundary_points_new_check = np.copy(boundary_points_new)
                while check_iscollinear(next_point,boundary_points_new_check[next_point_id],boundaries) == 0:
                    boundary_points_new_check[next_point_id] = [np.inf,np.inf]
                    next_point_id = cdist([next_point], boundary_points_new_check).argmin()
                    
                boundary_ridges_new[i,1] = next_point_id
                boundary_ridges_new[i+1,0] = next_point_id
                next_point = np.copy(boundary_points_new[next_point_id])
                boundary_points_new[next_point_id] = [np.inf,np.inf]
                
        boundary_ridges_new = (boundary_ridges_new+nvertices_in).astype(int) # shift the indices with "nvertices_in"
        
        voronoi_vertices = np.vstack((voronoi_vertices_in,boundary_points)) # vertical stack "in" vetices and "cross" boundary points
        nvertex = voronoi_vertices.shape[0]
        boundary_ridges_new = np.vstack((infinite_ridges_new,boundary_ridges_new))
        nboundary_ridge = boundary_ridges_new.shape[0]
        voronoi_ridges = np.vstack((finite_ridges_new,boundary_ridges_new))
        nridge = voronoi_ridges.shape[0]
    else:
        voronoi_vertices = np.vstack((voronoi_vertices_in,boundary_points)) # vertical stack "in" vetices and "cross" boundary points
        nvertex = voronoi_vertices.shape[0]
        boundary_ridges_new = np.copy(infinite_ridges_new)
        nboundary_ridge = boundary_ridges_new.shape[0]
        voronoi_ridges = np.vstack((finite_ridges_new,boundary_ridges_new))
        nridge = voronoi_ridges.shape[0]
    
    return voronoi_vertices,finite_ridges,boundary_points,finite_ridges_new,\
        boundary_ridges_new,nvertex,nvertices_in,nfinite_ridge,nboundary_ridge,\
            nboundary_pts,nboundary_pts_featured,voronoi_ridges,nridge

def RebuildVoronoi_merge(vor,circles,boundaries,generation_center,x_min,x_max,y_min,y_max,box_center,box_shape,merge_tol,boundaryFlag):
    """Clip Voronoi mesh by the boundaries, merge short Voronoi ridges, and rebuild the new Voronoi mesh"""
    # Store indices of Voronoi vertices for each finite ridge
    finite_ridges = []
    finite_ridges_pointid = []
    infinite_ridges = []
    infinite_ridges_pointid = []
    boundary_points = []
    
    # Form two new Voronoi vertices lists with and without out-of-bound vertices
    voronoi_vertices_in = []
    voronoi_vertices_out = []
    
    if box_shape == 'notched_square': # notched
        
        x_min = boundaries[0][1][0][0]
        x_max = boundaries[1][1][0][0]
        y_min = boundaries[0][1][0][1]
        y_max = boundaries[2][1][0][1]
        x_indent = boundaries[5][1][0][0]
        y_indent_min = boundaries[6][1][0][1]
        y_indent_max = boundaries[4][1][0][1]
        
        for pointidx, simplex in zip(vor.ridge_points, vor.ridge_vertices):
            simplex = np.asarray(simplex)
            simplex = np.where(simplex == -1, -9999999999, simplex) # use a large negative number instead of -1 to represent the infinite ridge vertices 
            if np.all(simplex >= 0) and check_isinside_boundbox2Dindent(np.append(vor.vertices[simplex[0]],0),x_min,x_max,y_min,y_max,x_indent,y_indent_min,y_indent_max) and check_isinside_boundbox2Dindent(np.append(vor.vertices[simplex[1]],0),x_min,x_max,y_min,y_max,x_indent,y_indent_min,y_indent_max):
                voronoi_vertices_in.append(vor.vertices[simplex[0]])
                voronoi_vertices_in.append(vor.vertices[simplex[1]])
                finite_ridges.append(simplex)
                finite_ridges_pointid.append(pointidx)
                plt.plot(vor.vertices[simplex, 0], vor.vertices[simplex, 1], 'k-')
            else: # infinite Voronoi ridges and finite ridges with out-of-bound vetices
                if np.all(simplex >= 0) and (check_isinside_boundbox2Dindent(np.append(vor.vertices[simplex[0]],0),x_min,x_max,y_min,y_max,x_indent,y_indent_min,y_indent_max) or check_isinside_boundbox2Dindent(np.append(vor.vertices[simplex[1]],0),x_min,x_max,y_min,y_max,x_indent,y_indent_min,y_indent_max)): # finite Voronoi ridge but one vertex is far
                    if (check_isinside_boundbox2Dindent(np.append(vor.vertices[simplex[0]],0),x_min,x_max,y_min,y_max,x_indent,y_indent_min,y_indent_max) == 0):# if vertex is out-of-bound, inverse its index
                        voronoi_vertices_out.append(vor.vertices[simplex[0]])
                        if simplex[0] == 0:
                            simplex[0] = -9999999998 # since inverse of 0 is still zero, we use another large number to represent -0
                        else:
                            simplex[0] = -simplex[0]
                    elif (check_isinside_boundbox2Dindent(np.append(vor.vertices[simplex[1]],0),x_min,x_max,y_min,y_max,x_indent,y_indent_min,y_indent_max) == 0):
                        voronoi_vertices_out.append(vor.vertices[simplex[1]])
                        if simplex[1] == 0:
                            simplex[1] = -9999999998 # since inverse of 0 is still zero, we use another large number to represent -0
                        else:
                            simplex[1] = -simplex[1]
                    infinite_ridges.append(simplex)
                    infinite_ridges_pointid.append(pointidx) 
                elif check_isinside_boundbox2Dindent(np.append(vor.vertices[simplex[simplex >= 0][0]],0),x_min,x_max,y_min,y_max,x_indent,y_indent_min,y_indent_max): # index of finite end Voronoi vertex of the infinite ridge
                    infinite_ridges.append(simplex)
                    infinite_ridges_pointid.append(pointidx) 
    else:
        for pointidx, simplex in zip(vor.ridge_points, vor.ridge_vertices):
            simplex = np.asarray(simplex)
            simplex = np.where(simplex == -1, -9999999999, simplex) # use a large negative number instead of -1 to represent the infinite ridge vertices 
            if np.all(simplex >= 0) and check_isinside_boundbox2D(np.append(vor.vertices[simplex[0]],0),x_min,x_max,y_min,y_max) and check_isinside_boundbox2D(np.append(vor.vertices[simplex[1]],0),x_min,x_max,y_min,y_max):
                voronoi_vertices_in.append(vor.vertices[simplex[0]])
                voronoi_vertices_in.append(vor.vertices[simplex[1]])
                finite_ridges.append(simplex)
                finite_ridges_pointid.append(pointidx)
            else: # infinite Voronoi ridges and finite ridges with out-of-bound vetices
                if np.all(simplex >= 0) and (check_isinside_boundbox2D(np.append(vor.vertices[simplex[0]],0),x_min,x_max,y_min,y_max) or check_isinside_boundbox2D(np.append(vor.vertices[simplex[1]],0),x_min,x_max,y_min,y_max)): # finite Voronoi ridge but one vertex is far
        
                    if (check_isinside_boundbox2D(np.append(vor.vertices[simplex[0]],0),x_min,x_max,y_min,y_max) == 0):# if vertex is out-of-bound, inverse its index
                        voronoi_vertices_out.append(vor.vertices[simplex[0]])
                        voronoi_vertices_in.append(vor.vertices[simplex[1]])
                        if simplex[0] == 0:
                            simplex[0] = -9999999998 # since inverse of 0 is still zero, we use another large number to represent -0
                        else:
                            simplex[0] = -simplex[0]
                    elif (check_isinside_boundbox2D(np.append(vor.vertices[simplex[1]],0),x_min,x_max,y_min,y_max) == 0):
                        voronoi_vertices_out.append(vor.vertices[simplex[1]])
                        voronoi_vertices_in.append(vor.vertices[simplex[0]])
                        if simplex[1] == 0:
                            simplex[1] = -9999999998 # since inverse of 0 is still zero, we use another large number to represent -0
                        else:
                            simplex[1] = -simplex[1]
                    infinite_ridges.append(simplex)
                    infinite_ridges_pointid.append(pointidx) 
                elif check_isinside_boundbox2D(np.append(vor.vertices[simplex[simplex >= 0][0]],0),x_min,x_max,y_min,y_max): # index of finite end Voronoi vertex of the infinite ridge
                    infinite_ridges.append(simplex)
                    infinite_ridges_pointid.append(pointidx) 
     
    voronoi_vertices_in = np.unique(voronoi_vertices_in,axis=0)
    voronoi_vertices_out = np.unique(voronoi_vertices_out,axis=0)
            
    # Find the intersect point of infinite ridges (rays) with the boundary lines
    # Ref: https://stackoverflow.com/questions/14307158/how-do-you-check-for-intersection-between-a-line-segment-and-a-line-ray-emanatin#:~:text=Let%20r%20%3D%20(cos%20%CE%B8%2C,0%20%E2%89%A4%20u%20%E2%89%A4%201).&text=Then%20your%20line%20segment%20intersects,0%20%E2%89%A4%20u%20%E2%89%A4%201.
    for pointidx, simplex in zip(infinite_ridges_pointid, infinite_ridges):
        simplex = np.asarray(simplex)
        i = simplex[simplex >= 0][0] # index of finite end Voronoi vertex of the infinite ridge
        p = vor.vertices[i]
        s = circles[pointidx[1],0:2] - circles[pointidx[0],0:2]
        if simplex[simplex < 0][0] == -9999999999:
            tangent = s / np.linalg.norm(s)
            normal = np.array([-tangent[0,1], tangent[0,0]]) # normal
            midpoint = circles[:,0:2][pointidx].mean(axis=0) # midpoint between site points
            for boundary in boundaries: # loop over boundary lines
                q = np.asarray(boundary[1][0])
                s = np.asarray(boundary[1][1]) - np.asarray(boundary[1][0])
                t = np.cross((q-p),s)/np.cross(normal,s)
                u = np.cross((q-p),normal)/np.cross(normal,s)
                
                if (u>= 0) and (u<=1):
                    if np.sign(np.dot(midpoint - generation_center, normal)) == 1: # facing outwards
                        if (t >= 0) and math.isfinite(t):
                            t_final = t
                    elif np.sign(np.dot(midpoint - generation_center, normal)) == -1: # facing inwards
                        if (t < 0) and math.isfinite(t):
                            t_final = t
                            
        elif simplex[simplex < 0][0] == -9999999998:
            t_final = np.inf
            normal = (vor.vertices[0] - p)/np.linalg.norm(vor.vertices[0] - p)
            for boundary in boundaries: # loop over boundary lines
                q = np.asarray(boundary[1][0])
                s = np.asarray(boundary[1][1]) - np.asarray(boundary[1][0])
                t = np.cross((q-p),s)/np.cross(normal,s)
                u = np.cross((q-p),normal)/np.cross(normal,s)
                
                if (u >= 0) and (u <= 1):
                    if (t > 0) and (t < t_final):
                        t_final = t
        else:
            t_final = np.inf
            normal = (vor.vertices[-simplex[simplex < 0][0]] - p)/np.linalg.norm(vor.vertices[-simplex[simplex < 0][0]] - p)
            for boundary in boundaries: # loop over boundary lines
                q = np.asarray(boundary[1][0])
                s = np.asarray(boundary[1][1]) - np.asarray(boundary[1][0])
                t = np.cross((q-p),s)/np.cross(normal,s)
                u = np.cross((q-p),normal)/np.cross(normal,s)
                
                if (u >= 0) and (u <= 1):
                    if (t > 0) and (t < t_final):
                        t_final = t
                            
        intersect_point = p + normal * t_final
        boundary_points.append(intersect_point)
            
    if boundaryFlag in ['on','On','Y','y','Yes','yes']:            
        boundary_points = np.asarray(boundary_points)
        boundary_points_featured = [x[1][0] for x in boundaries]
        boundary_points_featured = np.reshape(np.asarray(boundary_points_featured),(-1,2))
        nboundary_pts_ini = boundary_points.shape[0]
        nboundary_pts_featured = boundary_points_featured.shape[0]
        boundary_points = np.vstack((boundary_points,boundary_points_featured))
        
        nvertices_in = voronoi_vertices_in.shape[0]
        nboundary_pts = boundary_points.shape[0]
        finite_ridges = np.asarray(finite_ridges)
        nfinite_ridge = finite_ridges.shape[0]
        infinite_ridges = np.asarray(infinite_ridges)
        ninfinite_ridge = infinite_ridges.shape[0]
    else:
        boundary_points = np.asarray(boundary_points)
        nboundary_pts_featured = 0
        nvertices_in = voronoi_vertices_in.shape[0]
        nboundary_pts = boundary_points.shape[0]
        finite_ridges = np.asarray(finite_ridges)
        nfinite_ridge = finite_ridges.shape[0]
        infinite_ridges = np.asarray(infinite_ridges)
        ninfinite_ridge = infinite_ridges.shape[0]
        
    # reconstruct the connectivity for ridges since the unique operation rearrange the order of vertices (need to find a more efficient way like vectorize)
    finite_ridges_new = np.copy(finite_ridges)
    for i in range(0,nfinite_ridge): # loop over voronoi ridges in original finite_ridge list 
        ii = np.where(np.all(voronoi_vertices_in==vor.vertices[finite_ridges[i,0]],axis=1))[0][0]  # index of first vertex in voronoi_vertices_in array
        jj = np.where(np.all(voronoi_vertices_in==vor.vertices[finite_ridges[i,1]],axis=1))[0][0]  # index of second vertex in voronoi_vertices_in array
        finite_ridges_new[i,:] = [ii,jj]
            
    infinite_ridges_new = np.copy(infinite_ridges)
    for i in range(0,ninfinite_ridge): # loop over voronoi ridges in original infinite_ridge list
        if (infinite_ridges[i,:] < 0).argmax(axis=0) == 0: # first vertex is of negative index
            ii = nvertices_in + i
            jj = np.where(np.all(voronoi_vertices_in==vor.vertices[infinite_ridges[i,1]],axis=1))[0][0]  # index of finite vertex in voronoi_vertices_in array
        else:
            ii = np.where(np.all(voronoi_vertices_in==vor.vertices[infinite_ridges[i,0]],axis=1))[0][0]  # index of finite vertex in voronoi_vertices_in array
            jj = nvertices_in + i
        infinite_ridges_new[i,:] = [ii,jj]
    
    
    # Merge operation
    # Case 1 merge
    voronoi_vertices_in_whole = np.copy(voronoi_vertices_in)
    tree_voronoi_vertices_in = cKDTree(voronoi_vertices_in)
    rows_to_fuse_voronoi_vertices_in = tree_voronoi_vertices_in.query_pairs(r=merge_tol,output_type='ndarray')
    rows_to_fuse_voronoi_vertices_in_whole = np.copy(rows_to_fuse_voronoi_vertices_in)
        
    finite_ridges_merged = np.copy(finite_ridges_new)
    infinite_ridges_merged = np.copy(infinite_ridges_new)
    while len(rows_to_fuse_voronoi_vertices_in) != 0:
        points_new = np.reshape(voronoi_vertices_in_whole[rows_to_fuse_voronoi_vertices_in_whole],(-1,2))
        voronoi_vertices_in_fused = 0.5*(points_new[0::2] + points_new[1::2])
        
        # voronoi_vertices_in_fused = voronoi_vertices_in_fused.round(decimals=0)
        voronoi_vertices_in_fused = np.round(voronoi_vertices_in_fused/merge_tol).astype(int)*merge_tol # here we have to round to the nearest merge_tol number to avoid infinite loop

        # check if there are identical merged points from different pairs of closeby points (usually happens when three points are clustered)
        if len(np.unique(voronoi_vertices_in_fused,axis=0)) != len(voronoi_vertices_in_fused):
            [unique_pts,unique_inverse] = np.unique(voronoi_vertices_in_fused,return_inverse=True,axis=0)
            
            # Offset the indices of boundary points as we inserted points in vertice_in list
            infinite_ridges_merged[infinite_ridges_merged >= nvertices_in] += len(unique_pts)
            # Replace the indices of two close points with the index of new merged point in both finite_ridges_new and infinite_ridges_new lists
            for i in range(0,len(rows_to_fuse_voronoi_vertices_in_whole)):
                finite_ridges_merged[finite_ridges_merged == rows_to_fuse_voronoi_vertices_in_whole[i,0]] = nvertices_in + unique_inverse[i]
                finite_ridges_merged[finite_ridges_merged == rows_to_fuse_voronoi_vertices_in_whole[i,1]] = nvertices_in + unique_inverse[i]
                infinite_ridges_merged[infinite_ridges_merged == rows_to_fuse_voronoi_vertices_in_whole[i,0]] = nvertices_in + unique_inverse[i]
                infinite_ridges_merged[infinite_ridges_merged == rows_to_fuse_voronoi_vertices_in_whole[i,1]] = nvertices_in + unique_inverse[i]
            voronoi_vertices_in_fused = unique_pts.astype(float)
        else:
            # Offset the indices of boundary points as we inserted points in vertice_in list
            infinite_ridges_merged[infinite_ridges_merged >= nvertices_in] += len(rows_to_fuse_voronoi_vertices_in_whole)
            # Replace the indices of two close points with the index of new merged point in both finite_ridges_new and infinite_ridges_new lists
            for i in range(0,len(rows_to_fuse_voronoi_vertices_in_whole)):
                finite_ridges_merged[finite_ridges_merged == rows_to_fuse_voronoi_vertices_in_whole[i,0]] = nvertices_in + i
                finite_ridges_merged[finite_ridges_merged == rows_to_fuse_voronoi_vertices_in_whole[i,1]] = nvertices_in + i
                infinite_ridges_merged[infinite_ridges_merged == rows_to_fuse_voronoi_vertices_in_whole[i,0]] = nvertices_in + i
                infinite_ridges_merged[infinite_ridges_merged == rows_to_fuse_voronoi_vertices_in_whole[i,1]] = nvertices_in + i
                
        voronoi_vertices_in_whole = np.vstack((voronoi_vertices_in_whole,voronoi_vertices_in_fused)) # create a whole vertex list
        voronoi_vertices_in_new = np.vstack((np.delete(voronoi_vertices_in,rows_to_fuse_voronoi_vertices_in,0),voronoi_vertices_in_fused)) # create a new vertex list
        # update the voronoi_vertices_in
        voronoi_vertices_in = np.copy(voronoi_vertices_in_new)
        nvertices_in = voronoi_vertices_in_whole.shape[0]
        # update the kd-tree
        tree_voronoi_vertices_in = cKDTree(voronoi_vertices_in)
        rows_to_fuse_voronoi_vertices_in = tree_voronoi_vertices_in.query_pairs(r=merge_tol,output_type='ndarray')
        # update the rows_to_fuse_voronoi_vertices_in with the rows-to-fuse in list voronoi_vertices_in_whole
        rows_to_fuse_voronoi_vertices_in_whole = np.copy(rows_to_fuse_voronoi_vertices_in)
        for i in range(0,len(rows_to_fuse_voronoi_vertices_in)):
            rows_to_fuse_voronoi_vertices_in_whole[i,0] = np.where(np.all(voronoi_vertices_in_whole==voronoi_vertices_in_new[rows_to_fuse_voronoi_vertices_in[i,0],:],axis=1))[0][0]
            rows_to_fuse_voronoi_vertices_in_whole[i,1] = np.where(np.all(voronoi_vertices_in_whole==voronoi_vertices_in_new[rows_to_fuse_voronoi_vertices_in[i,1],:],axis=1))[0][0]

    # delete rows with identical indices and reverse indices in finite_ridges_merged
    zero_rows = np.where((finite_ridges_merged[:,0] - finite_ridges_merged[:,1]) == 0)[0]
    finite_ridges_merged = np.delete(finite_ridges_merged,zero_rows,0)
    
    # delete rows with identical indices and reverse indices in infinite_ridges_merged
    zero_rows = np.where((infinite_ridges_merged[:,0] - infinite_ridges_merged[:,1]) == 0)[0]
    infinite_ridges_merged = np.delete(infinite_ridges_merged,zero_rows,0)
    
    if boundaryFlag in ['on','On','Y','y','Yes','yes']: 
        # construct the connectivity for a line path consisting of the boundary points
        ####################################################################################################################################################
        # The algorithm here needs improve: bug contains when searching for the closest point, the closest point may be on the other sides
        ####################################################################################################################################################
        boundary_points_ini = np.copy(boundary_points) 
        boundary_ridges_new = np.zeros(boundary_points.shape)
        boundary_points_new = np.copy(boundary_points) 
        next_point = np.copy(boundary_points_new[0])  # get first point
        boundary_points_new[0] = [np.inf,np.inf]
        
        for i in range(0,boundary_points_new.shape[0]-1):
            next_point_id = cdist([next_point], boundary_points_new).argmin()
            boundary_ridges_new[i,1] = next_point_id
            boundary_ridges_new[i+1,0] = next_point_id
            next_point = np.copy(boundary_points_new[next_point_id])
            boundary_points_new[next_point_id] = [np.inf,np.inf]
        boundary_ridges_new = (boundary_ridges_new+nvertices_in).astype(int) # shift the indices with "nvertices_in"
        
        # Case 2 merge
        tree_boundary_points = cKDTree(boundary_points)
        rows_to_fuse_boundary_points = tree_boundary_points.query_pairs(r=merge_tol,output_type='ndarray')
        rows_to_fuse_boundary_points_whole = np.copy(rows_to_fuse_boundary_points)
        boundary_ridges_merged = np.copy(boundary_ridges_new)
        nboundary_pts = boundary_points.shape[0]
        boundary_points_whole = np.copy(boundary_points)
        
        while len(rows_to_fuse_boundary_points) != 0:
            points_new = np.reshape(boundary_points_whole[rows_to_fuse_boundary_points_whole],(-1,2))
            boundary_points_fused = 0.5*(points_new[0::2] + points_new[1::2])
            
            # boundary_points_fused = boundary_points_fused.round(decimals=0)
            boundary_points_fused = np.round(boundary_points_fused/merge_tol).astype(int)*merge_tol # here we have to round to the nearest merge_tol number to avoid infinite loop
            
            # check if there are identical merged points from different pairs of closeby points (usually happens when three points are clustered)
            if len(np.unique(boundary_points_fused,axis=0)) != len(boundary_points_fused):
                [unique_pts,unique_inverse] = np.unique(boundary_points_fused,return_inverse=True,axis=0)
                # Replace the indices of closeby points with the index of new merged point in both finite_ridges_new and infinite_ridges_new lists
                for i in range(0,len(rows_to_fuse_boundary_points)):
                    boundary_ridges_merged[boundary_ridges_merged == (rows_to_fuse_boundary_points_whole[i,0]+nvertices_in)] = nvertices_in + nboundary_pts + unique_inverse[i]
                    boundary_ridges_merged[boundary_ridges_merged == (rows_to_fuse_boundary_points_whole[i,1]+nvertices_in)] = nvertices_in + nboundary_pts + unique_inverse[i]
                    infinite_ridges_merged[infinite_ridges_merged == (rows_to_fuse_boundary_points_whole[i,0]+nvertices_in)] = nvertices_in + nboundary_pts + unique_inverse[i]
                    infinite_ridges_merged[infinite_ridges_merged == (rows_to_fuse_boundary_points_whole[i,1]+nvertices_in)] = nvertices_in + nboundary_pts + unique_inverse[i]
                boundary_points_fused = unique_pts.astype(float)
            else:
                # Replace the indices of two close points with the index of new merged point in both boundary_ridges_new and infinite_ridges_new lists
                for i in range(0,len(rows_to_fuse_boundary_points)):
                    boundary_ridges_merged[boundary_ridges_merged == (rows_to_fuse_boundary_points_whole[i,0]+nvertices_in)] = nvertices_in + nboundary_pts + i
                    boundary_ridges_merged[boundary_ridges_merged == (rows_to_fuse_boundary_points_whole[i,1]+nvertices_in)] = nvertices_in + nboundary_pts + i
                    infinite_ridges_merged[infinite_ridges_merged == (rows_to_fuse_boundary_points_whole[i,0]+nvertices_in)] = nvertices_in + nboundary_pts + i
                    infinite_ridges_merged[infinite_ridges_merged == (rows_to_fuse_boundary_points_whole[i,1]+nvertices_in)] = nvertices_in + nboundary_pts + i
                    
            boundary_points_whole = np.vstack((boundary_points_whole,boundary_points_fused)) # create a whole vertex list        
            boundary_points_new = np.vstack((np.delete(boundary_points,rows_to_fuse_boundary_points,0),boundary_points_fused)) # create a new vertex list
            
            
            # update the boundary_points
            boundary_points = np.copy(boundary_points_new)
            nboundary_pts = boundary_points_whole.shape[0]
            # update the kd-tree
            tree_boundary_points = cKDTree(boundary_points)
            rows_to_fuse_boundary_points = tree_boundary_points.query_pairs(r=merge_tol,output_type='ndarray')
            # update the rows_to_fuse in list voronoi_vertices_in with the rows-to-fuse in list voronoi_vertices_in_whole
            rows_to_fuse_boundary_points_whole = np.copy(rows_to_fuse_boundary_points)
            for i in range(0,len(rows_to_fuse_boundary_points)):
                rows_to_fuse_boundary_points_whole[i,0] = np.where(np.all(boundary_points_whole==boundary_points_new[rows_to_fuse_boundary_points[i,0],:],axis=1))[0][0]
                rows_to_fuse_boundary_points_whole[i,1] = np.where(np.all(boundary_points_whole==boundary_points_new[rows_to_fuse_boundary_points[i,1],:],axis=1))[0][0]
        
                    
        # delete rows with identical indices and reverse indices in finite_ridges_merged
        zero_rows = np.where((boundary_ridges_merged[:,0] - boundary_ridges_merged[:,1]) == 0)[0]
        boundary_ridges_merged = np.delete(boundary_ridges_merged,zero_rows,0)
        
        # delete rows with identical indices and reverse indices in infinite_ridges_merged
        zero_rows = np.where((infinite_ridges_merged[:,0] - infinite_ridges_merged[:,1]) == 0)[0]
        infinite_ridges_merged = np.delete(infinite_ridges_merged,zero_rows,0)
          
        # replace too close boundary points with the featured
        tree_boundary_points_whole = cKDTree(boundary_points_whole)
        for i in range(0,boundary_points_featured.shape[0]):
            rows_to_replace_boundary_points_whole = tree_boundary_points_whole.query_ball_point(boundary_points_featured[i,:],merge_tol)
            if len(rows_to_replace_boundary_points_whole) != 0: 
                # Visualize nearby points
                nearby_points = boundary_points_whole[rows_to_replace_boundary_points_whole]
                # replace indices
                for j in range(0,len(rows_to_replace_boundary_points_whole)):
                    boundary_ridges_merged[boundary_ridges_merged == (rows_to_replace_boundary_points_whole[j]+nvertices_in)] = nvertices_in + nboundary_pts_ini + i
                    infinite_ridges_merged[infinite_ridges_merged == (rows_to_replace_boundary_points_whole[j]+nvertices_in)] = nvertices_in + nboundary_pts_ini + i
        
        # delete rows with identical indices and reverse indices in boundary_ridges_merged
        zero_rows = np.where((boundary_ridges_merged[:,0] - boundary_ridges_merged[:,1]) == 0)[0]
        boundary_ridges_merged = np.delete(boundary_ridges_merged,zero_rows,0)
        
        # delete rows with identical indices and reverse indices in infinite_ridges_merged
        zero_rows = np.where((infinite_ridges_merged[:,0] - infinite_ridges_merged[:,1]) == 0)[0]
        infinite_ridges_merged = np.delete(infinite_ridges_merged,zero_rows,0)
        
        # Case 3 merge (after loop of case 1 and 2)
        for i in range(0,len(boundary_points)):
            rows_to_replace_voronoi_vertices_in = tree_voronoi_vertices_in.query_ball_point(boundary_points[i,:],merge_tol)
            if len(rows_to_replace_voronoi_vertices_in) != 0:
                # update the rows_to_place_voronoi_vertices_in with the rows-to-replace in list voronoi_vertices_in_whole
                rows_to_replace_voronoi_vertices_in_whole = np.copy(rows_to_replace_voronoi_vertices_in)
                for j in range(0,len(rows_to_replace_voronoi_vertices_in_whole)):
                    rows_to_replace_voronoi_vertices_in_whole[j] = np.where(np.all(voronoi_vertices_in_whole==voronoi_vertices_in[rows_to_replace_voronoi_vertices_in[j],:],axis=1))[0][0]
        
                # Visualize nearby points
                nearby_points = voronoi_vertices_in_whole[rows_to_replace_voronoi_vertices_in_whole]
                # replace indices
                for rows in rows_to_replace_voronoi_vertices_in_whole:
                    finite_ridges_merged[finite_ridges_merged == rows] = np.where(np.all(boundary_points_whole == boundary_points[i,:],axis=1))[0][0] + nvertices_in
                    infinite_ridges_merged[infinite_ridges_merged == rows] = np.where(np.all(boundary_points_whole == boundary_points[i,:],axis=1))[0][0] + nvertices_in

        zero_rows = np.where((finite_ridges_merged[:,0] - finite_ridges_merged[:,1]) == 0)[0]
        finite_ridges_merged = np.delete(finite_ridges_merged,zero_rows,0)
        
        # delete rows with identical indices and reverse indices in infinite_ridges_merged
        zero_rows = np.where((infinite_ridges_merged[:,0] - infinite_ridges_merged[:,1]) == 0)[0]
        infinite_ridges_merged = np.delete(infinite_ridges_merged,zero_rows,0)

        voronoi_vertices = np.vstack((voronoi_vertices_in,boundary_points)) # vertical stack "in" vetices and "cross" boundary points
        voronoi_vertices_whole = np.vstack((voronoi_vertices_in_whole,boundary_points_whole)) # vertical stack "in" vetices and "cross" boundary points

        all_ridges_merged = np.vstack((finite_ridges_merged,infinite_ridges_merged,boundary_ridges_merged))
        
    else:
        boundary_points_whole = np.copy(boundary_points)
        voronoi_vertices_whole = np.vstack((voronoi_vertices_in_whole,boundary_points_whole)) # vertical stack "in" vetices and "cross" boundary points
        all_ridges_merged = np.vstack((finite_ridges_merged,infinite_ridges_merged))
        
    # Visualize merged mesh   
    for vpair in all_ridges_merged:
        plt.plot([voronoi_vertices_whole[vpair[0],0], voronoi_vertices_whole[vpair[1],0]],[voronoi_vertices_whole[vpair[0],1], voronoi_vertices_whole[vpair[1],1]], 'g--')
    
    
    voronoi_vertices = np.unique(np.reshape(voronoi_vertices_whole[all_ridges_merged],(-1,2)),axis=0)
    # update the indices from voronoi_vertices_whole to voronoi_vertices(which is the list with no useless points)
    for i in range(0,len(finite_ridges_merged)):
        finite_ridges_merged[i,0] = np.where(np.all(voronoi_vertices==voronoi_vertices_whole[finite_ridges_merged[i,0],:],axis=1))[0][0]
        finite_ridges_merged[i,1] = np.where(np.all(voronoi_vertices==voronoi_vertices_whole[finite_ridges_merged[i,1],:],axis=1))[0][0]
    for i in range(0,len(infinite_ridges_merged)):
        infinite_ridges_merged[i,0] = np.where(np.all(voronoi_vertices==voronoi_vertices_whole[infinite_ridges_merged[i,0],:],axis=1))[0][0]
        infinite_ridges_merged[i,1] = np.where(np.all(voronoi_vertices==voronoi_vertices_whole[infinite_ridges_merged[i,1],:],axis=1))[0][0]
    if boundaryFlag in ['on','On','Y','y','Yes','yes']:  
        for i in range(0,len(boundary_ridges_merged)):
            boundary_ridges_merged[i,0] = np.where(np.all(voronoi_vertices==voronoi_vertices_whole[boundary_ridges_merged[i,0],:],axis=1))[0][0]
            boundary_ridges_merged[i,1] = np.where(np.all(voronoi_vertices==voronoi_vertices_whole[boundary_ridges_merged[i,1],:],axis=1))[0][0]
    # delete rows with identical indices and reverse indices in finite_ridges_merged
    zero_rows = np.where((finite_ridges_merged[:,0] - finite_ridges_merged[:,1]) == 0)[0]
    finite_ridges_merged = np.delete(finite_ridges_merged,zero_rows,0)
    # delete rows with identical indices and reverse indices in infinite_ridges_merged
    zero_rows = np.where((infinite_ridges_merged[:,0] - infinite_ridges_merged[:,1]) == 0)[0]
    infinite_ridges_merged = np.delete(infinite_ridges_merged,zero_rows,0)
    
    if boundaryFlag in ['on','On','Y','y','Yes','yes']: 
        # delete rows with identical indices and reverse indices in boundary_ridges_merged
        zero_rows = np.where((boundary_ridges_merged[:,0] - boundary_ridges_merged[:,1]) == 0)[0]
        boundary_ridges_merged = np.delete(boundary_ridges_merged,zero_rows,0)
    
    # voronoi_vertices = voronoi_vertices_whole
    nvertex = voronoi_vertices.shape[0]
    finite_ridges_new = finite_ridges_merged
    nfinite_ridge = finite_ridges_new.shape[0]
    if boundaryFlag in ['on','On','Y','y','Yes','yes']: 
        boundary_ridges_new = np.vstack((infinite_ridges_merged,boundary_ridges_merged))
        nboundary_ridge = boundary_ridges_new.shape[0]
        voronoi_ridges = np.vstack((finite_ridges_new,boundary_ridges_new))
        nridge = voronoi_ridges.shape[0]
        
        intersect_temp = intersect2D(np.unique(np.reshape(voronoi_vertices_whole[finite_ridges_merged],(-1,2)),axis=0),np.unique(np.reshape(voronoi_vertices_whole[boundary_ridges_merged],(-1,2)),axis=0))
        nvertices_in = np.unique(np.reshape(voronoi_vertices_whole[finite_ridges_merged],(-1,2)),axis=0).shape[0] - len(intersect_temp)
        nboundary_pts = np.unique(np.reshape(voronoi_vertices_whole[boundary_ridges_merged],(-1,2)),axis=0).shape[0]
    else:
        boundary_ridges_new = np.copy(infinite_ridges_merged)  
        nboundary_ridge = boundary_ridges_new.shape[0]
        voronoi_ridges = np.vstack((finite_ridges_new,boundary_ridges_new))
        nridge = voronoi_ridges.shape[0]
        
        nvertices_in = np.unique(np.reshape(voronoi_vertices_whole[finite_ridges_merged],(-1,2)),axis=0).shape[0]

    return voronoi_vertices,finite_ridges,boundary_points,finite_ridges_new,\
           boundary_ridges_new,nvertex,nvertices_in,nfinite_ridge,nboundary_ridge,\
           nboundary_pts,nboundary_pts_featured,voronoi_ridges,nridge
           

def RidgeMidQuarterPts(voronoi_vertices,nvertex,nvertices_in,voronoi_ridges,\
                       finite_ridges_new,boundary_ridges_new,nfinite_ridge,\
                       nboundary_ridge,nboundary_pts,nboundary_pts_featured):
    ######################### For finite Voronoi ridges ###########################
    # Form a list of middle points of finite Voronoi edges
    finite_ridge_mid = []
    finite_midpt_indices = []
    count = nvertex

    for vpair in finite_ridges_new:
        midpoint = (voronoi_vertices[vpair[0]] + voronoi_vertices[vpair[1]])/2
        finite_ridge_mid.append(midpoint)
        finite_midpt_indices.append(count)
        count += 1
        
    finite_ridge_mid = np.tile(finite_ridge_mid, (2,1)) # duplicate the mid point list
    finite_second_midpt_indices = [x+nfinite_ridge for x in finite_midpt_indices]
    finite_midpt_indices = np.concatenate((finite_midpt_indices,finite_second_midpt_indices))
    count += nfinite_ridge
    nfinite_midpt = nfinite_ridge*2
    
    # Form a list of quarter points of Voronoi edges
    finite_ridge_quarter = []
    finite_quarterpt_indices = []
    for vpair in finite_ridges_new:
        quarterpoint = 3./4*voronoi_vertices[vpair[0]] + 1./4*voronoi_vertices[vpair[1]]
        finite_ridge_quarter.append(quarterpoint)
        finite_quarterpt_indices.append(count)
        count += 1
    
    for vpair in finite_ridges_new:
        quarterpoint = 1./4*voronoi_vertices[vpair[0]] + 3./4*voronoi_vertices[vpair[1]]
        finite_ridge_quarter.append(quarterpoint)
        finite_quarterpt_indices.append(count)
        count += 1
    finite_quarterpt_indices = np.asarray(finite_quarterpt_indices)
    
    nfinite_quarterpt = nfinite_ridge*2
    
    nfinitept_per_layer = nvertices_in + nfinite_midpt + nfinite_quarterpt
    
    ######################### For boundary Voronoi ridges #########################
    # Form a list of mid-points of boundary Voronoi edges
    boundary_ridge_mid = []
    boundary_midpt_indices = []
    for vpair in boundary_ridges_new:
        boundary_midpoint = (voronoi_vertices[vpair[0]] + voronoi_vertices[vpair[1]])/2
        boundary_ridge_mid.append(boundary_midpoint)
        boundary_midpt_indices.append(count)
        count += 1
    
    boundary_ridge_mid = np.tile(boundary_ridge_mid, (2,1))
    boundary_second_midpt_indices = [x+nboundary_ridge for x in boundary_midpt_indices]
    boundary_midpt_indices = np.concatenate((boundary_midpt_indices,boundary_second_midpt_indices))
    count += nboundary_ridge
    nboundary_midpt = nboundary_ridge*2
    
    # Form a list of quarter-points of boundary Voronoi edges
    boundary_ridge_quarter = []
    boundary_quarterpt_indices = []
    for vpair in boundary_ridges_new:
        boundary_quarterpoint = 3./4*voronoi_vertices[vpair[0]] + 1./4*voronoi_vertices[vpair[1]]
        boundary_ridge_quarter.append(boundary_quarterpoint)
        boundary_quarterpt_indices.append(count)
        count += 1
    
    for vpair in boundary_ridges_new:
        boundary_quarterpoint = 1./4*voronoi_vertices[vpair[0]] + 3./4*voronoi_vertices[vpair[1]]
        boundary_ridge_quarter.append(boundary_quarterpoint)
        boundary_quarterpt_indices.append(count)
        count += 1
    boundary_quarterpt_indices = np.asarray(boundary_quarterpt_indices)
        
        
    nboundary_quarterpt = nboundary_ridge*2
    
    nboundary_pt_per_layer = nboundary_pts + nboundary_midpt + nboundary_quarterpt
    
    
    npt_per_layer = nvertices_in + nboundary_pts
    npt_per_layer_normal = npt_per_layer - nboundary_pts_featured
    npt_per_layer_vtk = nfinitept_per_layer + nboundary_pt_per_layer
    
    all_midpt_indices = np.vstack((np.reshape(finite_midpt_indices,(2,-1)).T,np.reshape(boundary_midpt_indices,(2,-1)).T))
    all_quarterpt_indices = np.vstack((np.reshape(finite_quarterpt_indices,(2,-1)).T,np.reshape(boundary_quarterpt_indices,(2,-1)).T))
    
    all_pts_2D = np.vstack((voronoi_vertices,finite_ridge_mid,finite_ridge_quarter,boundary_ridge_mid,boundary_ridge_quarter))
    all_ridges = np.hstack((voronoi_ridges,all_midpt_indices,all_quarterpt_indices))
    
    return all_pts_2D,all_ridges,npt_per_layer,npt_per_layer_normal,npt_per_layer_vtk



def VertexandRidgeinfo(all_pts_2D,all_ridges,npt_per_layer,npt_per_layer_normal,\
                       npt_per_layer_vtk,nridge,geoName,radii,generation_center,\
                       cellwallthickness_sparse,cellwallthickness_dense):
    
    """Generate the vertex and ridge info 
       Vertex info includes: coordinates, ridge indices, indices of another vertex for each ridge, ridge lengths, ridge angles, ridge width
       Ridge info includes: 2 vertex indices, 2 mid point indices, 2 quarter point indices """

    
    # Calculate lengths of all Voronoi ridges
    vector = all_pts_2D[all_ridges[:,1],:] - all_pts_2D[all_ridges[:,0],:]
    all_ridge_lengths = np.linalg.norm(vector, axis=1)
    
    # Calculate angles of all Voronoi ridges (angles measured counter-clock wise, x-axis --> (1,0), y-axis --> (0,1))
    all_ridge_angles = np.arctan2(vector[:,1],vector[:,0]) # np.arctan2(y, x) * 180 / np.pi = the angle

#======== This part can be expanded to allow the smooth transition of cell wall thickness =========
    # Case 1: Abrupt transition of cell wall thickness between dense/sparse cells
    # thicknesses of all Voronoi vertices (here identical thickness for all wings belonging to the same vertex is assumed)
    vertex_cellwallthickness_2D = np.zeros(npt_per_layer)
    
    vertex_distance2logcenter = np.linalg.norm(all_pts_2D - generation_center, axis=1) # calc distance between vertex and logcenter
    for i in range(0,npt_per_layer):
        if (bisect.bisect(radii,vertex_distance2logcenter[i]) % 2) != 0: # if odd, sparse cells
            vertex_cellwallthickness_2D[i] = cellwallthickness_sparse
        else: # if even, dense cells
            vertex_cellwallthickness_2D[i] = cellwallthickness_dense
#==================================================================================================
     
    # Generate a list containing info for each vertex
    all_vertices_info_2D = []
    all_direction_info_2D = []
    for i in range(0,npt_per_layer): # loop over all vertices
        
        allrows = np.where(all_ridges == i)[0] # find all Voronoi ridges containing the ith vertex
        vertex_info = []
        vertex_info.append(len(allrows)) # append total number of wings for the ith vertex
        vertex_info.append(allrows) # append ridge indices
        vertex_info.append(all_ridges[allrows,0:2][np.nonzero(all_ridges[allrows,0:2]-i)]) # append the indices of vetices at the other end of each ridge
        direction_info = []
        direction = np.where(all_ridges[allrows,0:2] == i)[1] # the position of current vertex in each ridge, p1 --> 0, p2 --> 1 
        direction = -np.power(-1,direction+1) # convert the way to express the direction of each ridge to: 1/-1 --> pointing from/to the current vertex
        direction_info.append(direction) # append ridge direction
        vertex_info.append(all_ridge_lengths[allrows]/2) # append wing lengths
        vertex_info.append(np.ones(len(allrows))*vertex_cellwallthickness_2D[i]) # append wing widths
        vertex_info.append(all_ridge_angles[allrows]+math.pi*(-direction.clip(max=0))) # append wing angles
        
        all_vertices_info_2D.append(vertex_info)
        all_direction_info_2D.append(direction_info)
    
    # flattened 1D numpy array of lengths, widths and angles of wings
    # ref:https://stackoverflow.com/questions/952914/how-do-i-make-a-flat-list-out-of-a-list-of-lists
    ridge_inices = [vertex_info[1] for vertex_info in all_vertices_info_2D]
    ridge_inices = np.array([item for sublist in ridge_inices for item in sublist])
    another_vertex_inices = [vertex_info[2] for vertex_info in all_vertices_info_2D]
    another_vertex_inices = np.array([item for sublist in another_vertex_inices for item in sublist])
    ridge_lengths = [vertex_info[3] for vertex_info in all_vertices_info_2D]
    ridge_lengths = np.array([item for sublist in ridge_lengths for item in sublist])
    ridge_widths = [vertex_info[4] for vertex_info in all_vertices_info_2D]
    ridge_widths = np.array([item for sublist in ridge_widths for item in sublist])
    ridge_angles = [vertex_info[5] for vertex_info in all_vertices_info_2D]
    ridge_angles = np.array([item for sublist in ridge_angles for item in sublist])
    ridge_directions = [direction_info[0] for direction_info in all_direction_info_2D]
    ridge_directions = np.array([item for sublist in ridge_directions for item in sublist])
    
    flattened_all_vertices_2D = np.vstack((ridge_inices,another_vertex_inices,ridge_lengths,ridge_widths,ridge_angles,ridge_directions))
 
    # Convert ragged info list to a rectangular shape numpy array
    max_wings = max([vertex_info[0] for vertex_info in all_vertices_info_2D])
    
    all_vertices_info_2D_nparray = np.zeros((npt_per_layer,1+max_wings*5))
    all_vertices_info_2D_nparray[:,0] = [vertex_info[0] for vertex_info in all_vertices_info_2D]
    
    for i in range(0,npt_per_layer):
        nwings = all_vertices_info_2D[i][0]
        all_vertices_info_2D_nparray[i,1:nwings+1] = all_vertices_info_2D[i][1]
        all_vertices_info_2D_nparray[i,max_wings+1:max_wings+nwings+1] = all_vertices_info_2D[i][2]
        all_vertices_info_2D_nparray[i,2*max_wings+1:2*max_wings+nwings+1] = all_vertices_info_2D[i][3]
        all_vertices_info_2D_nparray[i,3*max_wings+1:3*max_wings+nwings+1] = all_vertices_info_2D[i][4]
        all_vertices_info_2D_nparray[i,4*max_wings+1:4*max_wings+nwings+1] = all_vertices_info_2D[i][5]
    
    # Save info to txt files
    all_vertices_2D = np.hstack((all_pts_2D[0:npt_per_layer,:],all_vertices_info_2D_nparray))
    np.savetxt(Path('meshes/' + geoName + '/' + geoName +'-vertex.mesh'), all_vertices_2D, fmt='%.16g', delimiter=' '\
        ,header='Vertex Data Generated with RingsPy Mesh Generation Tool\n\
Number of vertices\n'+ str(npt_per_layer) +
    '\n\
Max number of wings for one vertex\n'+ str(max_wings) +
    '\n\
[xcoord ycoord nwings ridge1 ... ridgen farvertex1 ... farvertexn length1 ... lengthn width1 ... widthn angle1 ... anglen]', comments='')
    
    np.savetxt(Path('meshes/' + geoName + '/' + geoName +'-ridge.mesh'), all_ridges, fmt='%d', delimiter=' '\
        ,header='Ridge Data Generated with RingsPy Mesh Generation Tool\n\
Number of ridges\n'+ str(nridge) +
    '\n\
[vertex1 vertex2 midpt1 midpt2 qrtrpt1 qrtrpt2]', comments='')
    
    return all_vertices_2D, max_wings, flattened_all_vertices_2D, all_ridges


def GenerateBeamElement(NURBS_degree,nctrlpt_per_beam,segment_length,theta_min,theta_max,\
                        z_min,z_max,long_connector_ratio,npt_per_layer,voronoi_vertices,\
                        nvertex,voronoi_ridges,nridge,generation_center,all_vertices_2D,max_wings,\
                        flattened_all_vertices_2D,all_ridges):
    

    nctrlpt_per_elem = NURBS_degree + 1
    nbeam_per_grain = int(round((z_max-z_min)/segment_length))
    
    nlayer = nctrlpt_per_beam*nbeam_per_grain
    
    nconnector_t_per_beam = int((nctrlpt_per_beam-1)/NURBS_degree+1)
    nconnector_t_per_grain = int(nconnector_t_per_beam*nbeam_per_grain)
    
    theta = np.linspace(theta_min,theta_max,nlayer-(nbeam_per_grain-1))
    z_coord = np.linspace(z_min,z_max,nlayer-(nbeam_per_grain-1))
    connector_l_length = segment_length*long_connector_ratio
    connector_l_angle = (theta_max-theta_min)/nbeam_per_grain*long_connector_ratio
    
    # Insert repeated layers for the longitudinal connectors
    for i in range(nlayer-(nbeam_per_grain-1)-(nctrlpt_per_beam-1),1,-(nctrlpt_per_beam-1)): 
        theta = np.insert(theta,i,theta[i-1])
        theta[i:] += connector_l_angle
    for i in range(nlayer-(nbeam_per_grain-1)-(nctrlpt_per_beam-1),1,-(nctrlpt_per_beam-1)):
        z_coord = np.insert(z_coord,i,z_coord[i-1])
        z_coord[i:] += connector_l_length
    
    # Data for IGA use
    vertices_new = np.zeros((nlayer,npt_per_layer,3))
    for i in range(0,nlayer):
        for j in range(0,npt_per_layer):
            vertices_new[i,j,:2] = rotate_around_point_highperf(voronoi_vertices[j,:], theta[i], generation_center)
            vertices_new[i,j,2] = z_coord[i]
    
    # Vertex Data for IGA
    IGAvertices = np.reshape(vertices_new,(-1,3))
    
    # Connectivity for IGA
    npt_total = nlayer*npt_per_layer
    vertex_connectivity = np.linspace(0,npt_total-1,npt_total)
    
    # Beams
    ngrain = nvertex
    nbeam_total = ngrain*nbeam_per_grain
    beam_connectivity_original = np.zeros((nbeam_total,nctrlpt_per_beam))
    for i in range(0,nbeam_per_grain):
        for j in range(0, ngrain):
            irow = i*ngrain + j
            ivertex = i*npt_per_layer*nctrlpt_per_beam + j
            for icol in range(0,nctrlpt_per_beam):
                beam_connectivity_original[irow,icol] = ivertex+icol*npt_per_layer
    
    # Rearange beam connectivity such that each row is corresponding to a Bezier beam element
    beam_connectivity = np.copy(beam_connectivity_original)
    for ictrlpt in range(nctrlpt_per_beam-NURBS_degree,1,-NURBS_degree):
        beam_connectivity = np.insert(beam_connectivity,ictrlpt,beam_connectivity[:,ictrlpt-1],axis=1)
    
    beam_connectivity_original = (beam_connectivity_original+1).astype(int) # +1 because of in abaqus index starts from 1
    beam_connectivity = (np.reshape(beam_connectivity,(-1,nctrlpt_per_elem))+1).astype(int) # +1 because of in abaqus index starts from 1
    
    
    nbeamElem = beam_connectivity.shape[0]
    
    
    # Transverse connectors 
    nconnector_t = nridge*nconnector_t_per_grain
    connector_t_connectivity = np.zeros((nconnector_t,2))
    for i in range(0,nbeam_per_grain):
        for j in range(0,nconnector_t_per_beam):
            for k in range(0,nridge):
                irow = i*nconnector_t_per_beam*nridge + j*nridge + k
                connector_t_connectivity[irow,:] = voronoi_ridges[k]+i*npt_per_layer*nctrlpt_per_beam+j*NURBS_degree*npt_per_layer
                
    connector_t_connectivity = (connector_t_connectivity+1).astype(int)
    # Slicing block indices, https://stackoverflow.com/questions/39692769/efficient-numpy-indexing-take-first-n-rows-of-every-block-of-m-rows
    connector_t_index = np.linspace(0,nconnector_t-1,nconnector_t).astype(int)
    connector_t_bot_index = connector_t_index[np.mod(np.arange(connector_t_index.size),nridge*nconnector_t_per_beam)<nridge]
    connector_t_top_index = connector_t_index[np.mod(np.arange(connector_t_index.size),nridge*nconnector_t_per_beam)>=(nconnector_t_per_beam-1)*nridge]
    connector_t_reg_index = np.setdiff1d(np.setdiff1d(connector_t_index, connector_t_bot_index),connector_t_top_index)
    connector_t_bot_connectivity = np.copy(connector_t_connectivity)[connector_t_bot_index,:]
    connector_t_top_connectivity = np.copy(connector_t_connectivity)[connector_t_top_index,:]
    connector_t_reg_connectivity = np.copy(connector_t_connectivity)[connector_t_reg_index,:]

    
    # Longitudinal connectors
    nwings = all_vertices_2D[:,2].astype(int)
    nconnector_l_per_layer = sum(nwings)
    nconnector_l = nconnector_l_per_layer * (nbeam_per_grain-1)
    connector_l_connectivity = np.zeros((nconnector_l,2))
    connector_l_vertex_dict = np.zeros(nconnector_l)
    irow_conn_l = 0
    for i in range(0,nbeam_per_grain-1): # loop over layers of longitudinal connectors
        for j in range(0,ngrain): # loop over grains in each layer
            n = nwings[j]
            irow = i*ngrain + j
            for k in range(0,n): # loop over wings of each grain
                connector_l_connectivity[irow_conn_l,:] = (beam_connectivity_original[irow,-1],beam_connectivity_original[irow,-1]+npt_per_layer)
                connector_l_vertex_dict[irow_conn_l] = irow
                irow_conn_l += 1
    
    connector_l_vertex_dict = connector_l_vertex_dict.astype(int)
    connector_l_connectivity = connector_l_connectivity.astype(int)
    
    nconnector_total = nconnector_t + nconnector_l
    
    return IGAvertices,vertex_connectivity,beam_connectivity_original,nbeam_total,\
        beam_connectivity,nbeamElem,nlayer,connector_t_connectivity,\
        connector_t_bot_connectivity,connector_t_top_connectivity,\
        connector_t_reg_connectivity,connector_l_connectivity,nconnector_t_per_beam,\
        nconnector_t_per_grain,nconnector_t,nconnector_l,nconnector_total,\
        theta,z_coord,nbeam_per_grain,connector_l_vertex_dict


def ConnectorMeshFile(geoName,IGAvertices,connector_t_bot_connectivity,\
                   connector_t_reg_connectivity,connector_t_top_connectivity,\
                   height_connector_t,connector_l_connectivity,all_vertices_2D,\
                   max_wings,flattened_all_vertices_2D,nbeam_per_grain,nridge,\
                   connector_l_vertex_dict):
######### txt File stores the connector data for Abaqus analyses ##############
    nelem_connector_t_bot = connector_t_bot_connectivity.shape[0]
    nelem_connector_t_reg = connector_t_reg_connectivity.shape[0]
    nelem_connector_t_top = connector_t_top_connectivity.shape[0]
    nelem_connector_l = connector_l_connectivity.shape[0]
    nelem_total = nelem_connector_t_bot + nelem_connector_t_reg + nelem_connector_t_top + nelem_connector_l
    
    conn_l_ridge_index = flattened_all_vertices_2D[0,:].astype(int)
    conn_l_lengths = flattened_all_vertices_2D[2,:]
    # conn_l_widths = flattened_all_vertices_2D[3,:]  
    conn_l_angles = flattened_all_vertices_2D[4,:]
    rot = np.array([[np.cos(conn_l_angles), -np.sin(conn_l_angles)], [np.sin(conn_l_angles), np.cos(conn_l_angles)]]).T
    rot = np.hstack((rot[:,0,:],np.zeros(len(conn_l_angles))[:, np.newaxis])) # convert to 3D rotation matrix, assumes rotation remains still in-plane
    conn_l_tangents = rot
    
# Meshdata = [nodex1 nodey1 nodez1 nodex2 nodey2 nodez2 centerx centery centerz dx1 dy1 dz1 dx2 dy2 dz2 width height]    
    Meshdata = np.zeros((nelem_total,17))
    
    for i in range(0,nelem_connector_t_bot):
        Meshdata[i,0:3] = np.copy(IGAvertices)[connector_t_bot_connectivity[i,0]-1,:]
        Meshdata[i,3:6] = np.copy(IGAvertices)[connector_t_bot_connectivity[i,1]-1,:]
        Meshdata[i,16] = height_connector_t
    for i in range(0,nelem_connector_t_reg):
        Meshdata[i+nelem_connector_t_bot,0:3] = np.copy(IGAvertices)[connector_t_reg_connectivity[i,0]-1,:]
        Meshdata[i+nelem_connector_t_bot,3:6] = np.copy(IGAvertices)[connector_t_reg_connectivity[i,1]-1,:]
        Meshdata[i+nelem_connector_t_bot,16] = height_connector_t*2
    for i in range(0,nelem_connector_t_top):
        Meshdata[i+nelem_connector_t_bot+nelem_connector_t_reg,0:3] = np.copy(IGAvertices)[connector_t_top_connectivity[i,0]-1,:]
        Meshdata[i+nelem_connector_t_bot+nelem_connector_t_reg,3:6] = np.copy(IGAvertices)[connector_t_top_connectivity[i,1]-1,:]
        Meshdata[i+nelem_connector_t_bot+nelem_connector_t_reg,16] = height_connector_t
        
    for i in range(0,nelem_connector_l):
        Meshdata[i+nelem_connector_t_bot+nelem_connector_t_reg+nelem_connector_t_top,0:3] = np.copy(IGAvertices)[connector_l_connectivity[i,0]-1,:]
        Meshdata[i+nelem_connector_t_bot+nelem_connector_t_reg+nelem_connector_t_top,3:6] = np.copy(IGAvertices)[connector_l_connectivity[i,1]-1,:]
        Meshdata[i+nelem_connector_t_bot+nelem_connector_t_reg+nelem_connector_t_top,16] = conn_l_lengths[i]
        
    Meshdata[:,6:9] = (Meshdata[:,0:3] + Meshdata[:,3:6])/2
    Meshdata[0:nelem_connector_t_bot,8] = Meshdata[0:nelem_connector_t_bot,8]
    Meshdata[nelem_total-nelem_connector_t_top:nelem_total,8] = Meshdata[nelem_total-nelem_connector_t_top:nelem_total,8]
    Meshdata[:,9:12] = Meshdata[:,6:9] - Meshdata[:,0:3]
    Meshdata[:,12:15] = Meshdata[:,6:9] - Meshdata[:,3:6]

    # Add z-variation/random field for cellwall thicknesses/connector widths
    for i in range(0,nridge):
        Meshdata[i,15] = all_vertices_2D[connector_t_bot_connectivity[i,1]-1,3+max_wings*3] # assign widths to all bot transverse connectors
    Meshdata[:nelem_total-nelem_connector_l,15] = np.tile(Meshdata[0:nridge,15],3*nbeam_per_grain) # use the same widths for reg and top transverse connectors
    
    conn_l_widths = Meshdata[conn_l_ridge_index,15]
    Meshdata[nelem_total-nelem_connector_l:nelem_total,15] = conn_l_widths

    # add the eccentricity to the centers for longitudinal connectors (new center is correct)
    Meshdata[nelem_total-nelem_connector_l:nelem_total,6:9] += conn_l_tangents*Meshdata[nelem_total-nelem_connector_l:nelem_total,16][:, np.newaxis]/2

    # Replace nodal coordinates with nodal indices
    Meshdata[:,0:2] = np.concatenate((connector_t_bot_connectivity,connector_t_reg_connectivity,connector_t_top_connectivity,connector_l_connectivity))
    Meshdata = np.delete(Meshdata,[2,3,4,5],1)
    
    np.savetxt(Path('meshes/' + geoName + '/' + geoName+'-mesh.txt'), Meshdata, fmt='%.16g', delimiter=' '\
    ,header='# Connector Data Generated with RingsPy Mesh Generation Tool\n\
Number of bot connectors\n'+ str(nelem_connector_t_bot) +
'\n\
Number of reg connectors\n'+ str(nelem_connector_t_reg) +
'\n\
Number of top connectors\n'+ str(nelem_connector_t_top) +
'\n\
Number of long connectors\n'+ str(nelem_connector_l) +
'\n\
[inode jnode centerx centery centerz dx1 dy1 dz1 dx2 dy2 dz2 width height]', comments='')  

    return Meshdata


def insert_precracks(all_pts_2D,all_ridges,nridge,npt_per_layer,\
                     npt_per_layer_normal,npt_per_layer_vtk,\
                     nlayer,precrack_nodes,precrack_widths,\
                     cellsize_sparse,):
    
    precrack_midpts = (precrack_nodes[:,0:2]+precrack_nodes[:,2:4])/2.0
    ridge_midpts = all_pts_2D[all_ridges[:,2]]
    ridge_midpts_tree = cKDTree(ridge_midpts)
    near_ridges = []
    for midpt in precrack_midpts:
        near_ridges.append(ridge_midpts_tree.query_ball_point(midpt,3*cellsize_sparse))
        
    # Find the intersect point of neighboring ridges (lines) with the precrack line
    precracked_elem = []
    for i in range(0,len(precrack_nodes)):
        p = precrack_nodes[i,0:2]
        normal = (precrack_nodes[i,2:4] - p)/np.linalg.norm(precrack_nodes[i,2:4] - p)
        for ridge in near_ridges[i]:
            q = all_pts_2D[all_ridges[ridge,0],:]
            s = all_pts_2D[all_ridges[ridge,1],:] - q
            u = np.cross((q-p),normal)/np.cross(normal,s)
            if (u >= 0) and (u <= 1):
                precracked_elem.append(ridge)
                
    nconnector_t_precrack = len(precracked_elem)
    nconnector_l_precrack = 0

    # Visualize the precracks in the preview plot
    plt.plot(precrack_nodes[0,0::2],precrack_nodes[0,1::2],'r-',linewidth=2)
    
    return precracked_elem, nconnector_t_precrack, nconnector_l_precrack
        
        
def VisualizationFiles(geoName,NURBS_degree,nlayer,npt_per_layer_vtk,all_pts_2D,\
                       segment_length,theta,z_coord,nbeam_per_grain,nridge,\
                       voronoi_ridges,generation_center,all_ridges,nvertex,\
                       nconnector_t,nconnector_l,nctrlpt_per_beam,ConnMeshData,\
                       all_vertices_2D,max_wings,flattened_all_vertices_2D):
    ngrain = nvertex
    ninterval_per_beam_vtk = int((nctrlpt_per_beam-1)/2) # 2 layers ----> 1 interval

    nconnector_t_per_beam = int((nctrlpt_per_beam-1)/NURBS_degree+1)
    
    # Data for VTK use
    vertices_new = np.zeros((nlayer,npt_per_layer_vtk,3))
    for i in range(0,nlayer):
        for j in range(0,npt_per_layer_vtk):
            vertices_new[i,j,:2] = rotate_around_point_highperf(all_pts_2D[j,:], theta[i], generation_center)
            vertices_new[i,j,2] = z_coord[i]
    
    # Vertex Data for VTK use
    VTKvertices = np.reshape(vertices_new,(-1,3))
    
    # Connectivity for VTK use
    # Vertex
    npt_total_vtk = nlayer*npt_per_layer_vtk
    vertex_connectivity_vtk = np.linspace(0,npt_total_vtk-1,npt_total_vtk)
    
    # Beam
    nbeam_total_vtk = ngrain*nbeam_per_grain*ninterval_per_beam_vtk
    beam_connectivity_vtk = np.zeros((nbeam_total_vtk,3)) # 3 is the number of points per beam in Paraview
    for i in range(0,nbeam_per_grain):
        for j in range(0,ninterval_per_beam_vtk):
            for k in range(0,ngrain):
                irow = i*ninterval_per_beam_vtk*ngrain + j*ngrain + k
                ivertex = (i + (i*ninterval_per_beam_vtk+j)*2)*npt_per_layer_vtk + k
                beam_connectivity_vtk[irow,:] = (ivertex,ivertex+2*npt_per_layer_vtk,ivertex+npt_per_layer_vtk)
    
    # Transverse connectors 
    connector_t_connectivity_vtk = np.zeros((nconnector_t,2))
    for i in range(0,nbeam_per_grain):
        for j in range(0,nconnector_t_per_beam):
            for k in range(0,nridge):
                irow = i*nconnector_t_per_beam*nridge + j*nridge + k
                connector_t_connectivity_vtk[irow,:] = (voronoi_ridges[k][0]+i*npt_per_layer_vtk*nctrlpt_per_beam+j*NURBS_degree*npt_per_layer_vtk,voronoi_ridges[k][1]+i*npt_per_layer_vtk*nctrlpt_per_beam+j*NURBS_degree*npt_per_layer_vtk)
    
    # Longitudinal connectors
    connector_l_connectivity_vtk = np.zeros((nconnector_l,2))
    nwings = all_vertices_2D[:,2].astype(int)
    irow_conn_l = 0
    for i in range(0,nbeam_per_grain-1): # loop over layers of longitudinal connectors
        for j in range(0,ngrain): # loop over grains in each layer
            nw = nwings[j]
            irow = i*ngrain + j
            for k in range(0,nw): # loop over wings of each grain
                connector_l_connectivity_vtk[irow_conn_l,:] = (beam_connectivity_vtk[irow+(ninterval_per_beam_vtk-1)*ngrain,1],beam_connectivity_vtk[irow+(ninterval_per_beam_vtk-1)*ngrain,1]+npt_per_layer_vtk)
                irow_conn_l += 1
    
    # Quad - Beam Wings
    nquad = nridge*2
    nquad_total = nquad*nbeam_per_grain*ninterval_per_beam_vtk
    quad_connectivity = np.zeros((nquad_total,8))
    for i in range(0,nbeam_per_grain):
        for j in range(0,ninterval_per_beam_vtk):
            for k in range(0,nridge):
                for l in range(0,2):
                    irow = i*ninterval_per_beam_vtk*nquad + j*nquad + l*nridge + k
                    ivertex = i + (i*ninterval_per_beam_vtk+j)*2
                    quad_connectivity[irow,:] = (all_ridges[k][l]+ivertex*npt_per_layer_vtk,all_ridges[k][l]+(ivertex+2)*npt_per_layer_vtk, \
                                                 all_ridges[k][l+2]+(ivertex+2)*npt_per_layer_vtk,all_ridges[k][l+2]+ivertex*npt_per_layer_vtk, \
                                                 all_ridges[k][l]+(ivertex+1)*npt_per_layer_vtk,all_ridges[k][l+4]+(ivertex+2)*npt_per_layer_vtk, \
                                                 all_ridges[k][l+2]+(ivertex+1)*npt_per_layer_vtk,all_ridges[k][l+4]+ivertex*npt_per_layer_vtk)

    # Quad - Connector Cross-sections
    nquad_conn = nconnector_t + nconnector_l
    nquad_conn_total = nquad_conn
    
    conn_l_angles = flattened_all_vertices_2D[4,:]
    
    Quad_center = np.copy(ConnMeshData[:,2:5])
    Quad_normal1 = np.copy(ConnMeshData[:,5:8])
    Quad_normal2 = np.copy(ConnMeshData[:,8:11])
    Quad_length1 = np.linalg.norm(Quad_normal1, axis=1)
    Quad_length2 = np.linalg.norm(Quad_normal2, axis=1)
    Quad_width = np.copy(ConnMeshData[:,11])
    Quad_height = np.copy(ConnMeshData[:,12])
    Quad_normal = Quad_normal1/Quad_length1[:, np.newaxis]
    Quad_tangent = np.zeros((nquad_conn_total,3))
    Quad_tangent[0:nconnector_t] = np.tile(np.array([0.0,0.0,1.0]),(nconnector_t,1)) # assume a tangent of (0,0,1) for all conn_t
    
    rot = np.array([[np.cos(conn_l_angles), -np.sin(conn_l_angles)], [np.sin(conn_l_angles), np.cos(conn_l_angles)]]).T
    rot = np.hstack((rot[:,0,:],np.zeros(len(conn_l_angles))[:, np.newaxis])) # convert to 3D rotation matrix, assumes rotation remains still in-plane
    Quad_tangent[nconnector_t:nquad_conn] = rot
    Quad_bitangent = np.cross(Quad_normal,Quad_tangent)
    
    nconnector_t_bot = int(nconnector_t/3)
    nconnector_t_top = int(nconnector_t/3)
    # Add eccentricities to bot/top transverse connectors
    Quad_center[0:nconnector_t_bot,:] = Quad_center[0:nconnector_t_bot,:] + Quad_tangent[0:nconnector_t_bot,:]*Quad_height[0:nconnector_t_bot,np.newaxis]/2 # bot connectors
    Quad_center[(nconnector_t-nconnector_t_top):nconnector_t,:] = Quad_center[(nconnector_t-nconnector_t_top):nconnector_t,:] - Quad_tangent[(nconnector_t-nconnector_t_top):nconnector_t,:]*Quad_height[(nconnector_t-nconnector_t_top):nconnector_t,np.newaxis]/2
    
    quad_conn_vertices = np.zeros((nquad_conn_total*4,3))
    quad_conn_vertices[0:nquad_conn_total,:] = Quad_center - Quad_tangent*Quad_height[:, np.newaxis]/2 - Quad_bitangent*Quad_width[:, np.newaxis]/2
    quad_conn_vertices[nquad_conn_total:2*nquad_conn_total,:] = Quad_center - Quad_tangent*Quad_height[:, np.newaxis]/2 + Quad_bitangent*Quad_width[:, np.newaxis]/2
    quad_conn_vertices[2*nquad_conn_total:3*nquad_conn_total,:] = Quad_center + Quad_tangent*Quad_height[:, np.newaxis]/2 + Quad_bitangent*Quad_width[:, np.newaxis]/2
    quad_conn_vertices[3*nquad_conn_total:4*nquad_conn_total,:] = Quad_center + Quad_tangent*Quad_height[:, np.newaxis]/2 - Quad_bitangent*Quad_width[:, np.newaxis]/2
    
    quad_conn_connectivity = np.linspace(0,nquad_conn_total*4-1,nquad_conn_total*4)
    quad_conn_connectivity = np.reshape(quad_conn_connectivity, [4,-1]).T
    quad_conn_connectivity = quad_conn_connectivity + npt_total_vtk

    # Modify Vertex Data by adding quad_conn_vertices
    VTKvertices_quad_conn = np.vstack((VTKvertices,quad_conn_vertices))
    npt_total_vtk_quad_conn = VTKvertices_quad_conn.shape[0]
    vertex_connectivity_vtk_quad_conn = np.linspace(0,npt_total_vtk_quad_conn-1,npt_total_vtk_quad_conn)

    # Hex - Connector Volumes
    nhex_conn_total = nquad_conn_total*2
    
    hex_conn_vertices = np.zeros((nquad_conn_total*4*3,3))
    hex_conn_vertices[0:nquad_conn_total,:] = Quad_center - Quad_tangent*Quad_height[:, np.newaxis]/2 - Quad_bitangent*Quad_width[:, np.newaxis]/2
    hex_conn_vertices[nquad_conn_total:2*nquad_conn_total,:] = Quad_center - Quad_tangent*Quad_height[:, np.newaxis]/2 + Quad_bitangent*Quad_width[:, np.newaxis]/2
    hex_conn_vertices[2*nquad_conn_total:3*nquad_conn_total,:] = Quad_center + Quad_tangent*Quad_height[:, np.newaxis]/2 + Quad_bitangent*Quad_width[:, np.newaxis]/2
    hex_conn_vertices[3*nquad_conn_total:4*nquad_conn_total,:] = Quad_center + Quad_tangent*Quad_height[:, np.newaxis]/2 - Quad_bitangent*Quad_width[:, np.newaxis]/2
    
    hex_conn_vertices[4*nquad_conn_total:5*nquad_conn_total,:] = Quad_center - Quad_tangent*Quad_height[:, np.newaxis]/2 - Quad_bitangent*Quad_width[:, np.newaxis]/2 - Quad_normal*Quad_length1[:, np.newaxis]
    hex_conn_vertices[5*nquad_conn_total:6*nquad_conn_total,:] = Quad_center - Quad_tangent*Quad_height[:, np.newaxis]/2 + Quad_bitangent*Quad_width[:, np.newaxis]/2 - Quad_normal*Quad_length1[:, np.newaxis]
    hex_conn_vertices[6*nquad_conn_total:7*nquad_conn_total,:] = Quad_center + Quad_tangent*Quad_height[:, np.newaxis]/2 + Quad_bitangent*Quad_width[:, np.newaxis]/2 - Quad_normal*Quad_length1[:, np.newaxis]
    hex_conn_vertices[7*nquad_conn_total:8*nquad_conn_total,:] = Quad_center + Quad_tangent*Quad_height[:, np.newaxis]/2 - Quad_bitangent*Quad_width[:, np.newaxis]/2 - Quad_normal*Quad_length1[:, np.newaxis]

    hex_conn_vertices[8*nquad_conn_total:9*nquad_conn_total,:] = Quad_center - Quad_tangent*Quad_height[:, np.newaxis]/2 - Quad_bitangent*Quad_width[:, np.newaxis]/2 + Quad_normal*Quad_length2[:, np.newaxis]
    hex_conn_vertices[9*nquad_conn_total:10*nquad_conn_total,:] = Quad_center - Quad_tangent*Quad_height[:, np.newaxis]/2 + Quad_bitangent*Quad_width[:, np.newaxis]/2 + Quad_normal*Quad_length2[:, np.newaxis]
    hex_conn_vertices[10*nquad_conn_total:11*nquad_conn_total,:] = Quad_center + Quad_tangent*Quad_height[:, np.newaxis]/2 + Quad_bitangent*Quad_width[:, np.newaxis]/2 + Quad_normal*Quad_length2[:, np.newaxis]
    hex_conn_vertices[11*nquad_conn_total:12*nquad_conn_total,:] = Quad_center + Quad_tangent*Quad_height[:, np.newaxis]/2 - Quad_bitangent*Quad_width[:, np.newaxis]/2 + Quad_normal*Quad_length2[:, np.newaxis]
    
    hex_conn1_connectivity = np.linspace(0,nquad_conn_total*4*2-1,nquad_conn_total*4*2)
    hex_conn1_connectivity = np.reshape(hex_conn1_connectivity, [8,-1]).T
    hex_conn2_connectivity = np.linspace(0,nquad_conn_total*4*2-1,nquad_conn_total*4*2)
    hex_conn2_connectivity = np.reshape(hex_conn2_connectivity, [8,-1]).T
    hex_conn2_connectivity[:,4:8] = hex_conn2_connectivity[:,4:8] + nquad_conn_total*4
    hex_conn_connectivity = np.vstack((hex_conn1_connectivity,hex_conn2_connectivity))
    hex_conn_connectivity = hex_conn_connectivity + npt_total_vtk

    # Modify Vertex Data by adding hex_conn_vertices
    VTKvertices_hex_conn = np.vstack((VTKvertices,hex_conn_vertices))
    npt_total_vtk_hex_conn = VTKvertices_hex_conn.shape[0]
    vertex_connectivity_vtk_hex_conn = np.linspace(0,npt_total_vtk_hex_conn-1,npt_total_vtk_hex_conn)
    
# =============================================================================
    # Paraview Visualization File
    collocation_flag_vtk = np.concatenate((np.ones(ngrain),np.zeros(npt_per_layer_vtk*NURBS_degree-ngrain)))
    collocation_flag_vtk = np.concatenate((np.tile(collocation_flag_vtk, nconnector_t_per_beam-1),np.concatenate((np.ones(ngrain),np.zeros(npt_per_layer_vtk-ngrain)))))
    collocation_flag_vtk = np.tile(collocation_flag_vtk, nbeam_per_grain)

# =============================================================================
    # Paraview Vertices File
    VTKcell_types_vertices = (np.ones(npt_total_vtk)).astype(int)

    ncell_vertices = VTKcell_types_vertices.shape[0]

    vtkfile_vertices = open (Path('meshes/' + geoName + '/' + geoName + '_vertices'+'.vtu'),'w')
    
    vtkfile_vertices.write('<VTKFile type="UnstructuredGrid" version="2.0" byte_order="LittleEndian">'+'\n')
    vtkfile_vertices.write('<UnstructuredGrid>'+'\n')
    vtkfile_vertices.write('<Piece NumberOfPoints="'+str(npt_total_vtk)+'"'+' '+'NumberOfCells="'+str(ncell_vertices)+'">'+'\n')
    
    # <Points>
    vtkfile_vertices.write('<Points>'+'\n')
    vtkfile_vertices.write('<DataArray type="Float64" NumberOfComponents="3" format="ascii">'+'\n')
    for i in range(0,npt_total_vtk):
        X,Y,Z = VTKvertices[i]
        vtkfile_vertices.write(' '+'%11.8e'%X+'  '+'%11.8e'%Y+'  '+'%11.8e'%Z+'\n')
    vtkfile_vertices.write('</DataArray>'+'\n')
    vtkfile_vertices.write('</Points>'+'\n')
    # </Points>
    
    # <PointData> 
    vtkfile_vertices.write("<"+"PointData"\
        +" "+"Tensors="+'"'+""+'"'\
        +" "+"Vevtors="+'"'+""+'"'\
        +" "+"Scalars="+'"'+"IGAcollocation_flag"+'"'+">"+'\n')
    
    # Point Data
    vtkfile_vertices.write("<"+"DataArray"+" "+"type="+'"'+"Float32"+'"'+" "+"Name="+'"IGAcollocation_flag" format='+'"'+"ascii"+'"'+">"+'\n')
    for i in range(0,npt_total_vtk):
        X = collocation_flag_vtk[i]
        vtkfile_vertices.write('%11.8e'%X+'\n')
    vtkfile_vertices.write("</DataArray>"+'\n')
    
    vtkfile_vertices.write('</PointData>'+'\n')
    # </PointData> 
    
    
    # <Cells>
    vtkfile_vertices.write('<Cells>'+'\n')
    
    # Cell connectivity
    vtkfile_vertices.write('<DataArray type="Int32" Name="connectivity" format="ascii">'+'\n')
    for x in vertex_connectivity_vtk.astype(int):
        vtkfile_vertices.write("%s\n" % x)
    vtkfile_vertices.write('\n')  
    vtkfile_vertices.write('</DataArray>'+'\n')
    
    # Cell offsets
    vtkfile_vertices.write('<DataArray type="Int32" Name="offsets" format="ascii">'+'\n')
    current_offset = 0
    for element in vertex_connectivity_vtk:
        element_offset = 1
        current_offset += element_offset
        vtkfile_vertices.write(str(current_offset)+'\n')
    vtkfile_vertices.write('</DataArray>'+'\n')
    
    # Cell type
    vtkfile_vertices.write('<DataArray type="UInt8" Name="types" format="ascii">'+'\n')
    for i in range(0,ncell_vertices):
        element = VTKcell_types_vertices[i]
        vtkfile_vertices.write(str(element)+'\n')

    vtkfile_vertices.write('</DataArray>'+'\n')
    
    vtkfile_vertices.write('</Cells>'+'\n')
    # </Cells>
    
    vtkfile_vertices.write('</Piece>'+'\n')
    #</Piece>
    
    vtkfile_vertices.write('</UnstructuredGrid>'+'\n')
    #</UnstructuredGrid>
    
    vtkfile_vertices.write('</VTKFile>'+'\n')
    #</VTKFile>
    
    vtkfile_vertices.close()

# =============================================================================
    # Paraview Beam File
    VTKcell_types_beams = np.concatenate((21*np.ones(nbeam_total_vtk),23*np.ones(nquad_total))).astype(int)
  
    ncell_beams = VTKcell_types_beams.shape[0]

    vtkfile_beams = open (Path('meshes/' + geoName + '/' + geoName + '_beams'+'.vtu'),'w')
    
    vtkfile_beams.write('<VTKFile type="UnstructuredGrid" version="2.0" byte_order="LittleEndian">'+'\n')
    vtkfile_beams.write('<UnstructuredGrid>'+'\n')
    vtkfile_beams.write('<Piece NumberOfPoints="'+str(npt_total_vtk)+'"'+' '+'NumberOfCells="'+str(ncell_beams)+'">'+'\n')
    
    # <Points>
    vtkfile_beams.write('<Points>'+'\n')
    vtkfile_beams.write('<DataArray type="Float64" NumberOfComponents="3" format="ascii">'+'\n')
    for i in range(0,npt_total_vtk):
        X,Y,Z = VTKvertices[i]
        vtkfile_beams.write(' '+'%11.8e'%X+'  '+'%11.8e'%Y+'  '+'%11.8e'%Z+'\n')
    vtkfile_beams.write('</DataArray>'+'\n')
    vtkfile_beams.write('</Points>'+'\n')
    # </Points>
    
    # <PointData> 
    # </PointData> 
      
    # <Cells>
    vtkfile_beams.write('<Cells>'+'\n')
    
    # Cell connectivity
    vtkfile_beams.write('<DataArray type="Int32" Name="connectivity" format="ascii">'+'\n')
    vtkfile_beams.write("\n".join(" ".join(map(str, x)) for x in beam_connectivity_vtk.astype(int)))
    vtkfile_beams.write('\n')
    vtkfile_beams.write("\n".join(" ".join(map(str, x)) for x in quad_connectivity.astype(int)))
    vtkfile_beams.write('\n') 
    vtkfile_beams.write('</DataArray>'+'\n')
    
    # Cell offsets
    vtkfile_beams.write('<DataArray type="Int32" Name="offsets" format="ascii">'+'\n')
    current_offset = 0
    for element in beam_connectivity_vtk:
        element_offset = len(element)
        current_offset += element_offset
        vtkfile_beams.write(str(current_offset)+'\n')
    for element in quad_connectivity:
        element_offset = len(element)
        current_offset += element_offset
        vtkfile_beams.write(str(current_offset)+'\n')
    vtkfile_beams.write('</DataArray>'+'\n')
    
    # Cell type
    vtkfile_beams.write('<DataArray type="UInt8" Name="types" format="ascii">'+'\n')
    for i in range(0,ncell_beams):
        element = VTKcell_types_beams[i]
        vtkfile_beams.write(str(element)+'\n')

    vtkfile_beams.write('</DataArray>'+'\n')
    
    vtkfile_beams.write('</Cells>'+'\n')
    # </Cells>
    
    vtkfile_beams.write('</Piece>'+'\n')
    #</Piece>
    
    vtkfile_beams.write('</UnstructuredGrid>'+'\n')
    #</UnstructuredGrid>
    
    vtkfile_beams.write('</VTKFile>'+'\n')
    #</VTKFile>
    
    vtkfile_beams.close()
    
# =============================================================================
    # Paraview Connector (Axis + Center section) File
    VTKcell_types_conns = np.concatenate((3*np.ones(nconnector_t),3*np.ones(nconnector_l),9*np.ones(nconnector_t),9*np.ones(nconnector_l))).astype(int)
    ncell_conns = VTKcell_types_conns.shape[0]
    
    Quad_width_vtk = np.tile(Quad_width,(2))

    vtkfile_conns = open (Path('meshes/' + geoName + '/' + geoName + '_conns'+'.vtu'),'w')
    
    vtkfile_conns.write('<VTKFile type="UnstructuredGrid" version="2.0" byte_order="LittleEndian">'+'\n')
    vtkfile_conns.write('<UnstructuredGrid>'+'\n')
    vtkfile_conns.write('<Piece NumberOfPoints="'+str(npt_total_vtk_quad_conn)+'"'+' '+'NumberOfCells="'+str(ncell_conns)+'">'+'\n')
    
    # <Points>
    vtkfile_conns.write('<Points>'+'\n')
    vtkfile_conns.write('<DataArray type="Float64" NumberOfComponents="3" format="ascii">'+'\n')
    for i in range(0,npt_total_vtk_quad_conn):
        X,Y,Z = VTKvertices_quad_conn[i]
        vtkfile_conns.write(' '+'%11.8e'%X+'  '+'%11.8e'%Y+'  '+'%11.8e'%Z+'\n')
    vtkfile_conns.write('</DataArray>'+'\n')
    vtkfile_conns.write('</Points>'+'\n')
    # </Points>
    
    # <PointData> 
    # </PointData> 
      
    # <Cells>
    vtkfile_conns.write('<Cells>'+'\n')
    
    # Cell connectivity
    vtkfile_conns.write('<DataArray type="Int32" Name="connectivity" format="ascii">'+'\n')
    vtkfile_conns.write("\n".join(" ".join(map(str, x)) for x in connector_t_connectivity_vtk.astype(int)))
    vtkfile_conns.write('\n')
    vtkfile_conns.write("\n".join(" ".join(map(str, x)) for x in connector_l_connectivity_vtk.astype(int)))
    vtkfile_conns.write('\n')
    vtkfile_conns.write("\n".join(" ".join(map(str, x)) for x in quad_conn_connectivity.astype(int)))
    vtkfile_conns.write('\n')   
    vtkfile_conns.write('</DataArray>'+'\n')
    
    # Cell offsets
    vtkfile_conns.write('<DataArray type="Int32" Name="offsets" format="ascii">'+'\n')
    current_offset = 0
    for element in connector_t_connectivity_vtk:
        element_offset = len(element)
        current_offset += element_offset
        vtkfile_conns.write(str(current_offset)+'\n')
    for element in connector_l_connectivity_vtk:
        element_offset = len(element)
        current_offset += element_offset
        vtkfile_conns.write(str(current_offset)+'\n')
    for element in quad_conn_connectivity:
        element_offset = len(element)
        current_offset += element_offset
        vtkfile_conns.write(str(current_offset)+'\n')
    vtkfile_conns.write('</DataArray>'+'\n')
    
    # Cell type
    vtkfile_conns.write('<DataArray type="UInt8" Name="types" format="ascii">'+'\n')
    for i in range(0,ncell_conns):
        element = VTKcell_types_conns[i]
        vtkfile_conns.write(str(element)+'\n')

    vtkfile_conns.write('</DataArray>'+'\n')
    
    vtkfile_conns.write('</Cells>'+'\n')
    # </Cells>

    # <CellData>
    vtkfile_conns.write("<"+"CellData"\
            +" "+"Tensors="+'"'+""+'"'\
            +" "+"Vevtors="+'"'+""+'"'\
            +" "+"Scalars="+'"'+"Connector_width"+'"'+">"+'\n')
    
    vtkfile_conns.write("<"+"DataArray"+" "+"type="+'"'+"Float32"+'"'+" "+"Name="+'"Connector_width"'+" "+"format="+'"'+"ascii"+'"'+">"+'\n')
    for i in range(0,ncell_conns):
        X = Quad_width_vtk[i]
        vtkfile_conns.write('%11.8e'%X+'\n')
    vtkfile_conns.write('</DataArray>'+'\n')

    vtkfile_conns.write('</CellData>'+'\n')
        
    vtkfile_conns.write('</Piece>'+'\n')
    #</Piece>
    
    vtkfile_conns.write('</UnstructuredGrid>'+'\n')
    #</UnstructuredGrid>
    
    vtkfile_conns.write('</VTKFile>'+'\n')
    #</VTKFile>
    
    vtkfile_conns.close()
    
# =============================================================================
    # Paraview Connector (Volume) File
    VTKcell_types_conns_vol = np.concatenate((12*np.ones(nconnector_t),12*np.ones(nconnector_l),12*np.ones(nconnector_t),12*np.ones(nconnector_l))).astype(int)
    ncell_conns_vol = VTKcell_types_conns_vol.shape[0]
    
    Quad_width_vtk = np.tile(Quad_width,(2))

    vtkfile_conns_vol = open (Path('meshes/' + geoName + '/' + geoName + '_conns_vol'+'.vtu'),'w')
    
    vtkfile_conns_vol.write('<VTKFile type="UnstructuredGrid" version="2.0" byte_order="LittleEndian">'+'\n')
    vtkfile_conns_vol.write('<UnstructuredGrid>'+'\n')
    vtkfile_conns_vol.write('<Piece NumberOfPoints="'+str(npt_total_vtk_hex_conn)+'"'+' '+'NumberOfCells="'+str(ncell_conns_vol)+'">'+'\n')
    
    # <Points>
    vtkfile_conns_vol.write('<Points>'+'\n')
    vtkfile_conns_vol.write('<DataArray type="Float64" NumberOfComponents="3" format="ascii">'+'\n')
    for i in range(0,npt_total_vtk_hex_conn):
        X,Y,Z = VTKvertices_hex_conn[i]
        vtkfile_conns_vol.write(' '+'%11.8e'%X+'  '+'%11.8e'%Y+'  '+'%11.8e'%Z+'\n')
    vtkfile_conns_vol.write('</DataArray>'+'\n')
    vtkfile_conns_vol.write('</Points>'+'\n')
    # </Points>
    
    # <PointData> 
    # </PointData> 
      
    # <Cells>
    vtkfile_conns_vol.write('<Cells>'+'\n')
    
    # Cell connectivity
    vtkfile_conns_vol.write('<DataArray type="Int32" Name="connectivity" format="ascii">'+'\n')
    vtkfile_conns_vol.write("\n".join(" ".join(map(str, x)) for x in hex_conn_connectivity.astype(int)))
    vtkfile_conns_vol.write('\n')   
    vtkfile_conns_vol.write('</DataArray>'+'\n')
    
    # Cell offsets
    vtkfile_conns_vol.write('<DataArray type="Int32" Name="offsets" format="ascii">'+'\n')
    current_offset = 0
    for element in hex_conn_connectivity:
        element_offset = len(element)
        current_offset += element_offset
        vtkfile_conns_vol.write(str(current_offset)+'\n')
    vtkfile_conns_vol.write('</DataArray>'+'\n')
    
    # Cell type
    vtkfile_conns_vol.write('<DataArray type="UInt8" Name="types" format="ascii">'+'\n')
    for i in range(0,ncell_conns_vol):
        element = VTKcell_types_conns_vol[i]
        vtkfile_conns_vol.write(str(element)+'\n')

    vtkfile_conns_vol.write('</DataArray>'+'\n')
    
    vtkfile_conns_vol.write('</Cells>'+'\n')
    # </Cells>
    
    # <CellData>
    vtkfile_conns_vol.write("<"+"CellData"\
            +" "+"Tensors="+'"'+""+'"'\
            +" "+"Vectors="+'"'+""+'"'\
            +" "+"Scalars="+'"'+"Connector_width"+'"'+">"+'\n')
    
    vtkfile_conns_vol.write("<"+"DataArray"+" "+"type="+'"'+"Float32"+'"'+" "+"Name="+'"Connector_width"'+" "+"format="+'"'+"ascii"+'"'+">"+'\n')
    for i in range(0,ncell_conns_vol):
        X = Quad_width_vtk[i]
        vtkfile_conns_vol.write('%11.8e'%X+'\n')
    vtkfile_conns_vol.write('</DataArray>'+'\n')

    vtkfile_conns_vol.write('</CellData>'+'\n')
    
    vtkfile_conns_vol.write('</Piece>'+'\n')
    #</Piece>
    
    vtkfile_conns_vol.write('</UnstructuredGrid>'+'\n')
    #</UnstructuredGrid>
    
    vtkfile_conns_vol.write('</VTKFile>'+'\n')
    #</VTKFile>
    
    vtkfile_conns_vol.close()
    
    
def BezierExtraction(NURBS_degree,nctrlpt_per_beam,nbeam_total):
    knotVec_per_beam = np.concatenate((np.zeros(NURBS_degree),(np.linspace(0,1,int((nctrlpt_per_beam-1)/NURBS_degree+1))),np.ones(NURBS_degree)))
    knotVec = np.tile(knotVec_per_beam, (nbeam_total,1))
    return knotVec
        

def BezierBeamFile(geoName,NURBS_degree,nctrlpt_per_beam,\
                   nconnector_t_per_beam,npatch,knotVec):
############################# Abaqus txt File #################################
    txtfile = open (Path('meshes/' + geoName + '/' + geoName + 'IGA.txt'),'w')
    txtfile.write('# Dimension of beam elements \n')
    txtfile.write('1 \n')
    txtfile.write('# Order of basis function \n')
    txtfile.write('{:d} \n'.format(NURBS_degree))
    txtfile.write('# Number of control points per patch \n')
    txtfile.write('{:d} \n'.format(nctrlpt_per_beam))
    txtfile.write('# Number of elements per patch \n') 
    txtfile.write('{:d} \n'.format(nconnector_t_per_beam-1))
    txtfile.write('# Number of Patches \n') 
    txtfile.write('{:d} \n'.format(npatch))
    # Loop over patches
    for ipatch in range(0,npatch):
        txtfile.write('{:s} \n'.format('PATCH-'+str(ipatch+1)))
        txtfile.write('Size of knot vectors \n') 
        txtfile.write('{:d} \n'.format(knotVec.shape[1])) 
        txtfile.write('knot vectors \n')
        for j in range(0,knotVec.shape[1]):
            txtfile.write('{:f} '.format(knotVec[ipatch,j])) 
        txtfile.write('\n')
    txtfile.close()
    
            
def AbaqusFile(geoName,NURBS_degree,npatch,nbeam_per_grain,IGAvertices,beam_connectivity,\
               connector_t_bot_connectivity,connector_t_reg_connectivity,\
               connector_t_top_connectivity,connector_l_connectivity,nbeamElem,\
               nconnector_t,nconnector_l,nconnector_t_precrack,nconnector_l_precrack,\
               segment_length,height_connector_t,long_connector_ratio,\
               x_max,x_min,y_max,y_min,z_coord,box_shape,box_size,\
               cellwallthickness_sparse,cellwallthickness_dense,\
               merge_operation,merge_tol,\
               precrackFlag='off',precrack_elem=[]):
    
    # IGA beam parameters
    ninstance = 1
    nsvars_beam = 27 # number of svars for beam
    nsecgp = 4 # number of gauss points for beam sectional integration
    nsvars_secgp = 16 # number of svars at each gauss point
    iprops_beam = np.zeros(6)
    iprops_beam[0]    = 1 # section type index
    iprops_beam[1]    = 1 # number of instances
    iprops_beam[2]    = nconnector_t # number of transverse connectors
    iprops_beam[3]    = nconnector_t_precrack # number of precracked transverse connectors
    iprops_beam[4]    = nconnector_l # number of longitudinal connectors
    iprops_beam[5]    = nconnector_l_precrack # number of precracked longitudinal connectors
    iprops_beam = [int(x) for x in iprops_beam] 
    
    props_beam = np.zeros(16)
    props_beam[0]     = 1.5E-9 # Cell substance density [tonne/mm^3]
    props_beam[1]     = 1.75E+6 # Mesoscale elastic modulus [MPa]
    props_beam[2]     = 0.3 # Beam Poisson's ratio
    props_beam[3]     = cellwallthickness_sparse # Cross-sectional height [mm]
    props_beam[4]     = cellwallthickness_sparse # Cross-sectional width [mm]
    props_beam[5]     = 100 # Tensile Strength [MPa]
    props_beam[6]     = 200 # Tensile fracture energy [mJ/mm^2]
    props_beam[7]     = 4.1 # Shear Strength Ratio
    props_beam[8]     = 0.2 # Softening Exponent
    props_beam[9]     = 0.2 # Initial Friction
    props_beam[10]    = 0.0 # Asymptotic Friction
    props_beam[11]    = 600 # Transitional Stress [MPa]
    props_beam[12]    = 0.0 # Tensile Unloading
    props_beam[13]    = 0.0 # Shear Unloading
    props_beam[14]    = 0.0 # Shear Softening
    props_beam[15]    = 1.0 # Elastic Analysis Flag
    
    # Transverse connector parameters
    nsvars_conn_t = 32  # number of svars for transverse connector
    iprops_connector_t_bot = np.zeros(7)
    props_connector_t_bot = np.zeros(26)
    
    props_connector_t_bot[0]     = 1.5E-9 # Cell substance density [tonne/mm^3]
    props_connector_t_bot[1]     = 6.6E+8 # Mesoscale elastic modulus [MPa]
    props_connector_t_bot[2]     = 0.25 # Shear-Normal coupling coefficient
    props_connector_t_bot[3]     = height_connector_t # Connector height [mm]
    props_connector_t_bot[4]     = 0 # M-Distance [mm]
    props_connector_t_bot[5]     = -height_connector_t/2 # L-Distance [mm]
    props_connector_t_bot[6]     = 10.0 # Tensile Strength [MPa]
    props_connector_t_bot[7]     = 5.0 # Tensile characteristic length [mm] will be updated to # Tensile fracture energy [mJ/mm^2]
    props_connector_t_bot[8]     = 2.6 # Shear Strength Ratio
    props_connector_t_bot[9]     = 0.2 # Softening Exponent
    props_connector_t_bot[10]    = 0.2 # Initial Friction
    props_connector_t_bot[11]    = 0.0 # Asymptotic Friction
    props_connector_t_bot[12]    = 600 # Transitional Stress [MPa]
    props_connector_t_bot[12]    = 0.0 # Tensile Unloading
    props_connector_t_bot[14]    = 0.0 # Shear Unloading
    props_connector_t_bot[15]    = 0.0 # Shear Softening
    props_connector_t_bot[16]    = 0.0 # Elastic Analysis Flag
    props_connector_t_bot[17]    = 0.2 # Compressive Yielding Strength [MPa]
    props_connector_t_bot[18]    = 600 # Initial Hardening Modulus Ratio
    props_connector_t_bot[19]    = 0.0 # Transitional Strain Ratio
    props_connector_t_bot[20]    = 0.0 # Deviatoric Strain Threshold Ratio
    props_connector_t_bot[21]    = 0.0 # Deviatoric Damage Parameter
    props_connector_t_bot[22]    = 0.0 # Final Hardening Modulus Ratio
    props_connector_t_bot[23]    = 0.0 # Densification Ratio
    props_connector_t_bot[24]    = 0.0 # Volumetric Deviatoric Coupling
    props_connector_t_bot[25]    = 0.0 # Compressive Unloading
    
    iprops_connector_t_bot = [3,ninstance,nbeamElem,nconnector_t,nconnector_t_precrack,nconnector_l,nconnector_l_precrack]
    iprops_connector_t_bot = [int(x) for x in iprops_connector_t_bot] 
    iprops_connector_t_top = [2,ninstance,nbeamElem,nconnector_t,nconnector_t_precrack,nconnector_l,nconnector_l_precrack]
    iprops_connector_t_top = [int(x) for x in iprops_connector_t_top] 
    iprops_connector_t_reg = [1,ninstance,nbeamElem,nconnector_t,nconnector_t_precrack,nconnector_l,nconnector_l_precrack]
    iprops_connector_t_reg = [int(x) for x in iprops_connector_t_reg] 
    
    props_connector_t_reg = np.copy(props_connector_t_bot)
    props_connector_t_reg[3] = height_connector_t*2
    props_connector_t_reg[5] = 0
    props_connector_t_top = np.copy(props_connector_t_bot)
    props_connector_t_top[5] = height_connector_t/2
    
    # Longitudinal connector parameters
    nsvars_conn_l = 32  # number of svars for transverse connector
    iprops_connector_l = np.zeros(7)
    props_connector_l = np.zeros(24)
    # height_connector_l = segment_length/((nconnector_l_per_grain-1)*2)
    
    props_connector_l[0]     = 1.5E-9 # Cell substance density [tonne/mm^3]
    props_connector_l[1]     = 6.0E+06 # Mesoscale elastic modulus [MPa]
    props_connector_l[2]     = 0.25 # Shear-Normal coupling coefficient
    props_connector_l[3]     = cellwallthickness_sparse*cellwallthickness_dense # Connector sectional area [mm^2]
    props_connector_l[4]     = 20 # Tensile Strength [MPa]
    props_connector_l[5]     = long_connector_ratio*segment_length*1.05# 0.0105 # Tensile characteristic length [mm] will be updated to # Tensile fracture energy [mJ/mm^2]
    props_connector_l[6]     = 2.6 # Shear Strength Ratio
    props_connector_l[7]     = 0.2 # Softening Exponent
    props_connector_l[8]     = 0.2 # Initial Friction
    props_connector_l[9]     = 0.0 # Asymptotic Friction
    props_connector_l[10]    = 600 # Transitional Stress [MPa]
    props_connector_l[11]    = 0.0 # Tensile Unloading
    props_connector_l[12]    = 0.0 # Shear Unloading
    props_connector_l[13]    = 0.0 # Shear Softening
    props_connector_l[14]    = 0.0 # Elastic Analysis Flag
    props_connector_l[15]    = 0.2 # Compressive Yielding Strength [MPa]
    props_connector_l[16]    = 600 # Initial Hardening Modulus Ratio
    props_connector_l[17]    = 0.0 # Transitional Strain Ratio
    props_connector_l[18]    = 0.0 # Deviatoric Strain Threshold Ratio
    props_connector_l[19]    = 0.0 # Deviatoric Damage Parameter
    props_connector_l[20]    = 0.0 # Final Hardening Modulus Ratio
    props_connector_l[21]    = 0.0 # Densification Ratio
    props_connector_l[22]    = 0.0 # Volumetric Deviatoric Coupling
    props_connector_l[23]    = 0.0 # Compressive Unloading
    
    iprops_connector_l = [1,ninstance,nbeamElem,nconnector_t,nconnector_t_precrack,nconnector_l,nconnector_l_precrack]
    iprops_connector_l = [int(x) for x in iprops_connector_l] 
    
    # Strain rate effect parameters
    props_strainrate = np.zeros(4)
    props_strainrate[0] = 1.0 # Strain rate effect flag 
    props_strainrate[1] = 5.0 # Physical time scaling factor
    props_strainrate[2] = 1.0E-5 # Strain rate effect constant c0
    props_strainrate[3] = 5.0E-2 # Strain rate effect constant c1
    
    
    timestep     = 5.0E-8
    totaltime    = 5.0E-2
    z_min = np.min(z_coord)
    z_max = np.max(z_coord)
    # boundary_conditions = ['Hydrostatic']
    boundary_conditions = ['Bottom','Top','Left','Right','Front','Back']
    BC_velo_dof = 1 # 1-x, 2-y, 3-z
    BC_velo_value = box_size*0.03 # mm/s

    # Beam
    nelnode = NURBS_degree + 1
    noGPs = NURBS_degree + 1
    if NURBS_degree == 2:
        eltype = 201
    elif NURBS_degree == 3:
        eltype = 202
    else:
        print('Current NURBS degree is not supported. Check variable NURBS_degree.')
        exit()
    nelnode_connector = 2
    eltype_connector_t = 501
    eltype_connector_l = 601
    
    numnode = IGAvertices.shape[0]
    nelem = beam_connectivity.shape[0]
    nnode = beam_connectivity.shape[1]
    nelem_connector_t_bot = connector_t_bot_connectivity.shape[0]
    nelem_connector_t_reg = connector_t_reg_connectivity.shape[0]
    nelem_connector_t_top = connector_t_top_connectivity.shape[0]
    nelem_connector_l = connector_l_connectivity.shape[0]
    nconnector_t_precrack = iprops_connector_l[4]
    nconnector_l_precrack = iprops_connector_l[6]
    nelem_total = nelem + nelem_connector_t_bot + nelem_connector_t_reg + nelem_connector_t_top\
                + nelem_connector_l + nconnector_t_precrack + nconnector_l_precrack
    
    xaxis = 1
    yaxis = 2
    zaxis = 3

    # Generate a .inp file which can be directly imported and played in Abaqus
    meshinpfile = open(Path('meshes/' + geoName + '/' + geoName + '.inp'),'w')
    meshinpfile.write('*HEADING'+'\n')
    meshinpfile.write('** Job name: {0} Model name: Model-{1}\n'.format(geoName,geoName))
    meshinpfile.write('** Generated with RingsPy Mesh Generation Tool V{}\n'.format(pkg_resources.get_distribution("RingsPy").version))
    # PART
    meshinpfile.write('** \n')
    meshinpfile.write('** PARTS\n') 
    meshinpfile.write('** \n') 
    meshinpfile.write('*Part, name=Part-1\n')
    meshinpfile.write('*End Part\n')
    meshinpfile.write('** \n')
    # ASSEMBLY
    meshinpfile.write('** \n')
    meshinpfile.write('** ASSEMBLY\n') 
    meshinpfile.write('** \n') 
    meshinpfile.write('*Assembly, name=Assembly\n')
    meshinpfile.write('** \n') 
    # INSTANCE
    meshinpfile.write('*Instance, name=Part-1-1, part=Part-1\n')
    # nodes
    meshinpfile.write('*Node,nset=AllNodes\n')
    for i in range(0,numnode):
        meshinpfile.write('{:d}, {:#.9e}, {:#.9e},  {:#.9e}\n'.format(i+1,IGAvertices[i,0],IGAvertices[i,1],IGAvertices[i,2]))
    # beam element connectivity 
    count = 1
    meshinpfile.write('*Element, type=B32,elset=AllBeams\n') 
    for i in range(0,nelem):
        meshinpfile.write('{:d}'.format(count))
        count += 1
        for j in range(0,nnode):
            meshinpfile.write(', {:d}'.format(beam_connectivity[i,j])) 
        meshinpfile.write('\n')
        
    if precrackFlag in ['on','On','Y','y','Yes','yes']:
        connector_t_precracked_index = []
        connector_t_precracked_connectivity = []
        connector_l_precracked_index = []
        connector_l_precracked_connectivity = []
        # bottom transverse connector element connectivity 
        meshinpfile.write('*Element, type=T3D2,elset=AllBotConns\n') 
        for i in range(0,nelem_connector_t_bot):
            if i in precrack_elem:
                connector_t_precracked_index.append(count)
                count += 1
                connector_t_precracked_connectivity.append(connector_t_bot_connectivity[i,:])
            else:
                meshinpfile.write('{:d}'.format(count))
                count += 1
                for j in range(0,nelnode_connector):
                    meshinpfile.write(', {:d}'.format(connector_t_bot_connectivity[i,j])) 
                meshinpfile.write('\n')
        # regular transverse connector element connectivity 
        meshinpfile.write('*Element, type=T3D2,elset=AllRegConns\n') 
        for i in range(0,nelem_connector_t_reg):
            if i in precrack_elem:
                connector_t_precracked_index.append(count)
                count += 1
                connector_t_precracked_connectivity.append(connector_t_reg_connectivity[i,:])
            else:
                meshinpfile.write('{:d}'.format(count))
                count += 1
                for j in range(0,nelnode_connector):
                    meshinpfile.write(', {:d}'.format(connector_t_reg_connectivity[i,j])) 
                meshinpfile.write('\n')
        # top transverse connector element connectivity 
        meshinpfile.write('*Element, type=T3D2,elset=AllTopConns\n') 
        for i in range(0,nelem_connector_t_top):
            if i in precrack_elem:
                connector_t_precracked_index.append(count)
                count += 1
                connector_t_precracked_connectivity.append(connector_t_top_connectivity[i,:])
            else:
                meshinpfile.write('{:d}'.format(count))
                count += 1
                for j in range(0,nelnode_connector):
                    meshinpfile.write(', {:d}'.format(connector_t_top_connectivity[i,j])) 
                meshinpfile.write('\n')
        # precracked transverse connector element connectivity 
        meshinpfile.write('*Element, type=T3D2,elset=AllPrecrackTConns\n') 
        for i in range(0,len(connector_t_precracked_connectivity)):
            meshinpfile.write('{:d}'.format(connector_t_precracked_index[i]))
            for j in range(0,nelnode_connector):
                meshinpfile.write(', {:d}'.format(connector_t_precracked_connectivity[i][j]))
            meshinpfile.write('\n')
        # longitudinal connector element connectivity 
        meshinpfile.write('*Element, type=T3D2,elset=AllLongConns\n') 
        for i in range(0,nelem_connector_l):
            if i in precrack_elem:
                connector_l_precracked_index.append(count)
                count += 1
                connector_l_precracked_connectivity.append(connector_l_connectivity[i,:])
            else:
                meshinpfile.write('{:d}'.format(count))
                count += 1
                for j in range(0,nelnode_connector):
                    meshinpfile.write(', {:d}'.format(connector_l_connectivity[i,j])) 
                meshinpfile.write('\n')
        # precracked longitudinal connector element connectivity 
        meshinpfile.write('*Element, type=T3D2,elset=AllPrecrackLConns\n') 
        for i in range(0,len(connector_l_precracked_connectivity)):
            meshinpfile.write('{:d}'.format(connector_l_precracked_index[i]))
            for j in range(0,nelnode_connector):
                meshinpfile.write(', {:d}'.format(connector_l_precracked_connectivity[i][j]))
            meshinpfile.write('\n')
    else:
        # bottom transverse connector element connectivity 
        meshinpfile.write('*Element, type=T3D2,elset=AllBotConns\n') 
        for i in range(0,nelem_connector_t_bot):
            meshinpfile.write('{:d}'.format(count))
            count += 1
            for j in range(0,nelnode_connector):
                meshinpfile.write(', {:d}'.format(connector_t_bot_connectivity[i,j])) 
            meshinpfile.write('\n')
        # regular transverse connector element connectivity 
        meshinpfile.write('*Element, type=T3D2,elset=AllRegConns\n') 
        for i in range(0,nelem_connector_t_reg):
            meshinpfile.write('{:d}'.format(count))
            count += 1
            for j in range(0,nelnode_connector):
                meshinpfile.write(', {:d}'.format(connector_t_reg_connectivity[i,j])) 
            meshinpfile.write('\n')
        # top transverse connector element connectivity 
        meshinpfile.write('*Element, type=T3D2,elset=AllTopConns\n') 
        for i in range(0,nelem_connector_t_top):
            meshinpfile.write('{:d}'.format(count))
            count += 1
            for j in range(0,nelnode_connector):
                meshinpfile.write(', {:d}'.format(connector_t_top_connectivity[i,j])) 
            meshinpfile.write('\n')
        # longitudinal connector element connectivity 
        meshinpfile.write('*Element, type=T3D2,elset=AllLongConns\n') 
        for i in range(0,nelem_connector_l):
            meshinpfile.write('{:d}'.format(count))
            count += 1
            for j in range(0,nelnode_connector):
                meshinpfile.write(', {:d}'.format(connector_l_connectivity[i,j])) 
            meshinpfile.write('\n')

    # Ghost mesh for easy visualization
    count_offset = 10**(int(math.log10(count))+1) # set an offset with the order of magnitude of the max number + 1
    count = count_offset + 1
    
    if NURBS_degree == 2:
        meshinpfile.write('*ELEMENT, TYPE=B32, ELSET=VisualBeams\n')
        for i in range(0,nelem):
            meshinpfile.write('{:d}'.format(count))
            count += 1
            for j in range(0,nnode):
                meshinpfile.write(', {:d}'.format(beam_connectivity[i,j])) 
            meshinpfile.write('\n')
            
    if precrackFlag in ['on','On','Y','y','Yes','yes']:
        visual_connector_precracked_index = []
        meshinpfile.write('*ELEMENT, TYPE=T3D2, ELSET=VisualBotConns\n')
        for i in range(0,nelem_connector_t_bot):
            if i in precrack_elem:
                visual_connector_precracked_index.append(count)
                count += 1
            else:
                meshinpfile.write('{:d}'.format(count))
                count += 1
                for j in range(0,nelnode_connector):
                    meshinpfile.write(', {:d}'.format(connector_t_bot_connectivity[i,j])) 
                meshinpfile.write('\n')
        meshinpfile.write('*ELEMENT, TYPE=T3D2, ELSET=VisualRegConns\n')
        for i in range(0,nelem_connector_t_reg):
            if i in precrack_elem:
                visual_connector_precracked_index.append(count)
                count += 1
            else:
                meshinpfile.write('{:d}'.format(count))
                count += 1
                for j in range(0,nelnode_connector):
                    meshinpfile.write(', {:d}'.format(connector_t_reg_connectivity[i,j])) 
                meshinpfile.write('\n')
        meshinpfile.write('*ELEMENT, TYPE=T3D2, ELSET=VisualTopConns\n')
        for i in range(0,nelem_connector_t_top):
            if i in precrack_elem:
                visual_connector_precracked_index.append(count)
                count += 1
            else:
                meshinpfile.write('{:d}'.format(count))
                count += 1
                for j in range(0,nelnode_connector):
                    meshinpfile.write(', {:d}'.format(connector_t_top_connectivity[i,j])) 
                meshinpfile.write('\n')
        # precracked transverse connector element connectivity 
        meshinpfile.write('*ELEMENT, TYPE=T3D2, ELSET=VisualPrecrackTConns\n')
        for i in range(0,len(connector_t_precracked_connectivity)):
            meshinpfile.write('{:d}'.format(visual_connector_precracked_index[i]))
            for j in range(0,nelnode_connector):
                meshinpfile.write(', {:d}'.format(connector_t_precracked_connectivity[i][j]))
            meshinpfile.write('\n')
        # longitudinal connector element connectivity 
        meshinpfile.write('*Element, type=T3D2,elset=VisualLConns\n') 
        for i in range(0,nelem_connector_l):
            if i in precrack_elem:
                connector_l_precracked_index.append(count)
                count += 1
                connector_l_precracked_connectivity.append(connector_l_connectivity[i,:])
            else:
                meshinpfile.write('{:d}'.format(count))
                count += 1
                for j in range(0,nelnode_connector):
                    meshinpfile.write(', {:d}'.format(connector_l_connectivity[i,j])) 
                meshinpfile.write('\n')
        # precracked longitudinal connector element connectivity 
        meshinpfile.write('*Element, type=T3D2,elset=VisualPrecrackLConns\n')
        meshinpfile.write('\n')
        # for i in range(0,len(connector_l_precracked_connectivity)):
        #     meshinpfile.write('{:d}'.format(connector_l_precracked_index[i]))
        #     for j in range(0,nelnode_connector):
        #         meshinpfile.write(', {:d}'.format(connector_l_precracked_connectivity[i][j]))
        #     meshinpfile.write('\n')      
    else:
        meshinpfile.write('*ELEMENT, TYPE=T3D2, ELSET=VisualBotConns\n')
        for i in range(0,nelem_connector_t_bot):
            meshinpfile.write('{:d}'.format(count))
            count += 1
            for j in range(0,nelnode_connector):
                meshinpfile.write(', {:d}'.format(connector_t_bot_connectivity[i,j])) 
            meshinpfile.write('\n')
        meshinpfile.write('*ELEMENT, TYPE=T3D2, ELSET=VisualRegConns\n')
        for i in range(0,nelem_connector_t_reg):
            meshinpfile.write('{:d}'.format(count))
            count += 1
            for j in range(0,nelnode_connector):
                meshinpfile.write(', {:d}'.format(connector_t_reg_connectivity[i,j])) 
            meshinpfile.write('\n')
        meshinpfile.write('*ELEMENT, TYPE=T3D2, ELSET=VisualTopConns\n')
        for i in range(0,nelem_connector_t_top):
            meshinpfile.write('{:d}'.format(count))
            count += 1
            for j in range(0,nelnode_connector):
                meshinpfile.write(', {:d}'.format(connector_t_top_connectivity[i,j])) 
            meshinpfile.write('\n')
        meshinpfile.write('*ELEMENT, TYPE=T3D2, ELSET=VisualLongConns\n')
        for i in range(0,nelem_connector_l):
            meshinpfile.write('{:d}'.format(count))
            count += 1
            for j in range(0,nelnode_connector):
                meshinpfile.write(', {:d}'.format(connector_l_connectivity[i,j])) 
            meshinpfile.write('\n')
            
    meshinpfile.write('** Section: Section-1\n')
    meshinpfile.write('*Solid Section, elset=VisualBotConns, material=VisualTConns\n')
    meshinpfile.write('{:8.4E}\n'.format(props_connector_t_bot[3]*cellwallthickness_sparse))
    meshinpfile.write('** Section: Section-2\n')
    meshinpfile.write('*Solid Section, elset=VisualRegConns, material=VisualTConns\n')
    meshinpfile.write('{:8.4E}\n'.format(props_connector_t_bot[3]*cellwallthickness_sparse))
    meshinpfile.write('** Section: Section-3\n')
    meshinpfile.write('*Solid Section, elset=VisualTopConns, material=VisualTConns\n')
    meshinpfile.write('{:8.4E}\n'.format(props_connector_t_bot[3]*cellwallthickness_sparse))
    meshinpfile.write('** Section: Section-11\n')
    meshinpfile.write('*Solid Section, elset=AllBotConns, material=VisualTConns\n')
    meshinpfile.write('{:8.4E}\n'.format(props_connector_t_bot[3]*cellwallthickness_sparse))
    meshinpfile.write('** Section: Section-12\n')
    meshinpfile.write('*Solid Section, elset=AllRegConns, material=VisualTConns\n')
    meshinpfile.write('{:8.4E}\n'.format(props_connector_t_bot[3]*cellwallthickness_sparse))
    meshinpfile.write('** Section: Section-13\n')
    meshinpfile.write('*Solid Section, elset=AllTopConns, material=VisualTConns\n')
    meshinpfile.write('{:8.4E}\n'.format(props_connector_t_bot[3]*cellwallthickness_sparse))
    if precrackFlag in ['on','On','Y','y','Yes','yes']:
        meshinpfile.write('** Section: Section-4\n')
        meshinpfile.write('*Solid Section, elset=VisualPrecrackTConns, material=VisualTConns\n')
        meshinpfile.write('{:8.4E}\n'.format(props_connector_t_bot[3]*cellwallthickness_sparse))
        meshinpfile.write('** Section: Section-5\n')
        meshinpfile.write('*Solid Section, elset=VisualLongConns, material=VisualLConns\n')
        meshinpfile.write('{:8.4E}\n'.format(props_connector_l[3]))
        meshinpfile.write('** Section: Section-6\n')
        meshinpfile.write('*Solid Section, elset=VisualPrecrackLConns, material=VisualLConns\n')
        meshinpfile.write('{:8.4E}\n'.format(props_connector_l[3]))
        meshinpfile.write('** Section: Section-14\n')
        meshinpfile.write('*Solid Section, elset=AllPrecrackTConns, material=VisualTConns\n')
        meshinpfile.write('{:8.4E}\n'.format(props_connector_t_bot[3]*cellwallthickness_sparse))
        meshinpfile.write('** Section: Section-15\n')
        meshinpfile.write('*Solid Section, elset=AllLongConns, material=VisualLConns\n')
        meshinpfile.write('{:8.4E}\n'.format(props_connector_l[3]))
        meshinpfile.write('** Section: Section-16\n')
        meshinpfile.write('*Solid Section, elset=AllPrecrackLConns, material=VisualLConns\n')
        meshinpfile.write('{:8.4E}\n'.format(props_connector_l[3]))
        if NURBS_degree == 2:
            meshinpfile.write('** Section: Section-7  Profile: Profile-1\n')
            meshinpfile.write('*Beam Section, elset=VisualBeams, material=VisualBeams, temperature=GRADIENTS, section=RECT\n')
            meshinpfile.write('{:8.4E}, {:8.4E}\n'.format(props_beam[3],props_beam[4]))
            meshinpfile.write('1.,0.,0.\n')
            meshinpfile.write('** Section: Section-17  Profile: Profile-1\n')
            meshinpfile.write('*Beam Section, elset=AllBeams, material=VisualBeams, temperature=GRADIENTS, section=RECT\n')
            meshinpfile.write('{:8.4E}, {:8.4E}\n'.format(props_beam[3],props_beam[4]))
            meshinpfile.write('1.,0.,0.\n')  
    else:
        meshinpfile.write('** Section: Section-4\n')
        meshinpfile.write('*Solid Section, elset=VisualLongConns, material=VisualLConns\n')
        meshinpfile.write('{:8.4E}\n'.format(props_connector_l[3]))
        meshinpfile.write('** Section: Section-14\n')
        meshinpfile.write('*Solid Section, elset=AllLongConns, material=VisualLConns\n')
        meshinpfile.write('{:8.4E}\n'.format(props_connector_l[3]))
        if NURBS_degree == 2:
            meshinpfile.write('** Section: Section-5  Profile: Profile-1\n')
            meshinpfile.write('*Beam Section, elset=VisualBeams, material=VisualBeams, temperature=GRADIENTS, section=RECT\n')
            meshinpfile.write('{:8.4E}, {:8.4E}\n'.format(props_beam[3],props_beam[4]))
            meshinpfile.write('1.,0.,0.\n')
            meshinpfile.write('** Section: Section-15  Profile: Profile-1\n')
            meshinpfile.write('*Beam Section, elset=AllBeams, material=VisualBeams, temperature=GRADIENTS, section=RECT\n')
            meshinpfile.write('{:8.4E}, {:8.4E}\n'.format(props_beam[3],props_beam[4]))
            meshinpfile.write('1.,0.,0.\n')
            
    meshinpfile.write('*End Instance\n') 
    meshinpfile.write('** \n')
    # NODE SETS
    meshinpfile.write('*Nset, nset=AllNodes, instance=Part-1-1, generate\n') 
    meshinpfile.write('{:d}, {:d}, {:d} \n'.format(1,numnode,1))
    
    if any(x in geoName for x in ['hydrostatic_', 'Hydrostatic_']):
        # boundary nodes
        if box_shape == 'square':
            if any(item in boundary_conditions for item in ['left','Left','L']):
                # Nodes on the left
                LeftNodes = (np.where(IGAvertices[:,xaxis-1] <= x_min)[0]+1).reshape(1,-1)
                meshinpfile.write('*Nset, nset=LeftNodes, instance=Part-1-1\n')
                for i in range(0,math.ceil(len(LeftNodes[0])/15)): # Abaqus only accepts maximum 15 items per row
                    meshinpfile.write(''.join(','.join(map(str, LeftNodes[0][15*i:15*(i+1)])))+',\n')
            if any(item in boundary_conditions for item in ['right','Right','R']):
                # Nodes on the right
                RightNodes = (np.where(IGAvertices[:,xaxis-1] >= x_max)[0]+1).reshape(1,-1)
                meshinpfile.write('*Nset, Nset=RightNodes, instance=Part-1-1\n')
                for i in range(0,math.ceil(len(RightNodes[0])/15)): # Abaqus only accepts maximum 15 items per row
                    meshinpfile.write(''.join(','.join(map(str, RightNodes[0][15*i:15*(i+1)])))+',\n')
            if any(item in boundary_conditions for item in ['bottom','Bottom','Bot']):
                # Nodes on the bottom
                BottomNodes = (np.where(IGAvertices[:,zaxis-1] <= z_min)[0]+1).reshape(1,-1)
                meshinpfile.write('*Nset, nset=BottomNodes, instance=Part-1-1\n')
                for i in range(0,math.ceil(len(BottomNodes[0])/15)): # Abaqus only accepts maximum 15 items per row
                    meshinpfile.write(''.join(','.join(map(str, BottomNodes[0][15*i:15*(i+1)])))+',\n')
            if any(item in boundary_conditions for item in ['top','Top','T']):
                # Nodes on the top
                TopNodes = (np.where(IGAvertices[:,zaxis-1] >= z_max)[0]+1).reshape(1,-1)
                meshinpfile.write('*Nset, nset=TopNodes, instance=Part-1-1\n')
                for i in range(0,math.ceil(len(TopNodes[0])/15)): # Abaqus only accepts maximum 15 items per row
                    meshinpfile.write(''.join(','.join(map(str, TopNodes[0][15*i:15*(i+1)])))+',\n')
            if any(item in boundary_conditions for item in ['back','Back']):
                # Nodes on the back
                BackNodes = (np.where(IGAvertices[:,yaxis-1] <= y_min)[0]+1).reshape(1,-1)
                meshinpfile.write('*Nset, nset=BackNodes, instance=Part-1-1\n')
                for i in range(0,math.ceil(len(BackNodes[0])/15)): # Abaqus only accepts maximum 15 items per row
                    meshinpfile.write(''.join(','.join(map(str, BackNodes[0][15*i:15*(i+1)])))+',\n')
            if any(item in boundary_conditions for item in ['front','Front','F']):
                # Nodes on the front
                FrontNodes = (np.where(IGAvertices[:,yaxis-1] >= y_max)[0]+1).reshape(1,-1)
                meshinpfile.write('*Nset, nset=FrontNodes, instance=Part-1-1\n')
                for i in range(0,math.ceil(len(FrontNodes[0])/15)): # Abaqus only accepts maximum 15 items per row
                    meshinpfile.write(''.join(','.join(map(str, FrontNodes[0][15*i:15*(i+1)])))+',\n')
        elif box_shape == 'hexagon':      
            if any(item in boundary_conditions for item in ['Hydrostatic','hydrostatic']):
                offset = x_max*0.05
                # Nodes on the bottom-left
                BottomLeftNodes = (np.where(IGAvertices[:,yaxis-1] <= (-np.sqrt(3)*IGAvertices[:,xaxis-1] - np.sqrt(3)*x_max + offset))[0]+1).reshape(1,-1)
                meshinpfile.write('*Nset, nset=BottomLeftNodes, instance=Part-1-1\n')
                for i in range(0,math.ceil(len(BottomLeftNodes[0])/15)): # Abaqus only accepts maximum 15 items per row
                    meshinpfile.write(''.join(','.join(map(str, BottomLeftNodes[0][15*i:15*(i+1)])))+',\n')
                # Nodes on the bottom-right
                BottomRightNodes = (np.where(IGAvertices[:,yaxis-1] <= (np.sqrt(3)*IGAvertices[:,xaxis-1] - np.sqrt(3)*x_max + offset))[0]+1).reshape(1,-1)
                meshinpfile.write('*Nset, Nset=BottomRightNodes, instance=Part-1-1\n')
                for i in range(0,math.ceil(len(BottomRightNodes[0])/15)): # Abaqus only accepts maximum 15 items per row
                    meshinpfile.write(''.join(','.join(map(str, BottomRightNodes[0][15*i:15*(i+1)])))+',\n')
                # Nodes on the bottom
                BottomNodes = (np.where(IGAvertices[:,yaxis-1] <= y_min)[0]+1).reshape(1,-1)
                meshinpfile.write('*Nset, nset=BottomNodes, instance=Part-1-1\n')
                for i in range(0,math.ceil(len(BottomNodes[0])/15)): # Abaqus only accepts maximum 15 items per row
                    meshinpfile.write(''.join(','.join(map(str, BottomNodes[0][15*i:15*(i+1)])))+',\n')
                # Nodes on the top-left
                TopLeftNodes = (np.where(IGAvertices[:,yaxis-1] >= (np.sqrt(3)*IGAvertices[:,xaxis-1] + np.sqrt(3)*x_max - offset))[0]+1).reshape(1,-1)
                meshinpfile.write('*Nset, nset=TopLeftNodes, instance=Part-1-1\n')
                for i in range(0,math.ceil(len(TopLeftNodes[0])/15)): # Abaqus only accepts maximum 15 items per row
                    meshinpfile.write(''.join(','.join(map(str, TopLeftNodes[0][15*i:15*(i+1)])))+',\n')
                # Nodes on the top-right
                TopRightNodes = (np.where(IGAvertices[:,yaxis-1] >= (-np.sqrt(3)*IGAvertices[:,xaxis-1] + np.sqrt(3)*x_max - offset))[0]+1).reshape(1,-1)
                meshinpfile.write('*Nset, nset=TopRightNodes, instance=Part-1-1\n')
                for i in range(0,math.ceil(len(TopRightNodes[0])/15)): # Abaqus only accepts maximum 15 items per row
                    meshinpfile.write(''.join(','.join(map(str, TopRightNodes[0][15*i:15*(i+1)])))+',\n')
                # Nodes on the top
                TopNodes = (np.where(IGAvertices[:,yaxis-1] >= y_max)[0]+1).reshape(1,-1)
                meshinpfile.write('*Nset, nset=TopNodes, instance=Part-1-1\n')
                for i in range(0,math.ceil(len(TopNodes[0])/15)): # Abaqus only accepts maximum 15 items per row
                    meshinpfile.write(''.join(','.join(map(str, TopNodes[0][15*i:15*(i+1)])))+',\n')
                    
                # Nodes on the front
                FrontNodes = (np.where(IGAvertices[:,zaxis-1] >= z_max)[0]+1).reshape(1,-1)
                meshinpfile.write('*Nset, nset=FrontNodes, instance=Part-1-1\n')
                for i in range(0,math.ceil(len(FrontNodes[0])/15)): # Abaqus only accepts maximum 15 items per row
                    meshinpfile.write(''.join(','.join(map(str, FrontNodes[0][15*i:15*(i+1)])))+',\n')
                # Nodes on the back
                BackNodes = (np.where(IGAvertices[:,zaxis-1] <= z_min)[0]+1).reshape(1,-1)
                meshinpfile.write('*Nset, nset=BackNodes, instance=Part-1-1\n')
                for i in range(0,math.ceil(len(BackNodes[0])/15)): # Abaqus only accepts maximum 15 items per row
                    meshinpfile.write(''.join(','.join(map(str, BackNodes[0][15*i:15*(i+1)])))+',\n')
            
    elif any(x in geoName for x in ['uniaxial_', 'Uniaxial_']):
        # boundary nodes
        if box_shape == 'square':
            if any(item in boundary_conditions for item in ['left','Left','L']):
                # Nodes on the left
                LeftNodes = (np.where(IGAvertices[:,xaxis-1] <= x_min)[0]+1).reshape(1,-1)
                meshinpfile.write('*Nset, nset=LeftNodes, instance=Part-1-1\n')
                for i in range(0,math.ceil(len(LeftNodes[0])/15)): # Abaqus only accepts maximum 15 items per row
                    meshinpfile.write(''.join(','.join(map(str, LeftNodes[0][15*i:15*(i+1)])))+',\n')
            if any(item in boundary_conditions for item in ['right','Right','R']):
                # Nodes on the right
                RightNodes = (np.where(IGAvertices[:,xaxis-1] >= x_max)[0]+1).reshape(1,-1)
                meshinpfile.write('*Nset, Nset=RightNodes, instance=Part-1-1\n')
                for i in range(0,math.ceil(len(RightNodes[0])/15)): # Abaqus only accepts maximum 15 items per row
                    meshinpfile.write(''.join(','.join(map(str, RightNodes[0][15*i:15*(i+1)])))+',\n')
            if any(item in boundary_conditions for item in ['bottom','Bottom','Bot']):
                # Nodes on the bottom
                BottomNodes = (np.where(IGAvertices[:,zaxis-1] <= z_min)[0]+1).reshape(1,-1)
                meshinpfile.write('*Nset, nset=BottomNodes, instance=Part-1-1\n')
                for i in range(0,math.ceil(len(BottomNodes[0])/15)): # Abaqus only accepts maximum 15 items per row
                    meshinpfile.write(''.join(','.join(map(str, BottomNodes[0][15*i:15*(i+1)])))+',\n')
            if any(item in boundary_conditions for item in ['top','Top','T']):
                # Nodes on the top
                TopNodes = (np.where(IGAvertices[:,zaxis-1] >= z_max)[0]+1).reshape(1,-1)
                meshinpfile.write('*Nset, nset=TopNodes, instance=Part-1-1\n')
                for i in range(0,math.ceil(len(TopNodes[0])/15)): # Abaqus only accepts maximum 15 items per row
                    meshinpfile.write(''.join(','.join(map(str, TopNodes[0][15*i:15*(i+1)])))+',\n')
            if any(item in boundary_conditions for item in ['back','Back']):
                # Nodes on the back
                BackNodes = (np.where(IGAvertices[:,yaxis-1] <= y_min)[0]+1).reshape(1,-1)
                meshinpfile.write('*Nset, nset=BackNodes, instance=Part-1-1\n')
                for i in range(0,math.ceil(len(BackNodes[0])/15)): # Abaqus only accepts maximum 15 items per row
                    meshinpfile.write(''.join(','.join(map(str, BackNodes[0][15*i:15*(i+1)])))+',\n')
            if any(item in boundary_conditions for item in ['front','Front','F']):
                # Nodes on the front
                FrontNodes = (np.where(IGAvertices[:,yaxis-1] >= y_max)[0]+1).reshape(1,-1)
                meshinpfile.write('*Nset, nset=FrontNodes, instance=Part-1-1\n')
                for i in range(0,math.ceil(len(FrontNodes[0])/15)): # Abaqus only accepts maximum 15 items per row
                    meshinpfile.write(''.join(','.join(map(str, FrontNodes[0][15*i:15*(i+1)])))+',\n')
        elif box_shape == 'hexagon':      
            if any(item in boundary_conditions for item in ['Hydrostatic','hydrostatic']):
                offset = x_max*0.05
                # Nodes on the bottom-left
                BottomLeftNodes = (np.where(IGAvertices[:,yaxis-1] <= (-np.sqrt(3)*IGAvertices[:,xaxis-1] - np.sqrt(3)*x_max + offset))[0]+1).reshape(1,-1)
                meshinpfile.write('*Nset, nset=BottomLeftNodes, instance=Part-1-1\n')
                for i in range(0,math.ceil(len(BottomLeftNodes[0])/15)): # Abaqus only accepts maximum 15 items per row
                    meshinpfile.write(''.join(','.join(map(str, BottomLeftNodes[0][15*i:15*(i+1)])))+',\n')
                # Nodes on the bottom-right
                BottomRightNodes = (np.where(IGAvertices[:,yaxis-1] <= (np.sqrt(3)*IGAvertices[:,xaxis-1] - np.sqrt(3)*x_max + offset))[0]+1).reshape(1,-1)
                meshinpfile.write('*Nset, Nset=BottomRightNodes, instance=Part-1-1\n')
                for i in range(0,math.ceil(len(BottomRightNodes[0])/15)): # Abaqus only accepts maximum 15 items per row
                    meshinpfile.write(''.join(','.join(map(str, BottomRightNodes[0][15*i:15*(i+1)])))+',\n')
                # Nodes on the bottom
                BottomNodes = (np.where(IGAvertices[:,yaxis-1] <= y_min)[0]+1).reshape(1,-1)
                meshinpfile.write('*Nset, nset=BottomNodes, instance=Part-1-1\n')
                for i in range(0,math.ceil(len(BottomNodes[0])/15)): # Abaqus only accepts maximum 15 items per row
                    meshinpfile.write(''.join(','.join(map(str, BottomNodes[0][15*i:15*(i+1)])))+',\n')
                # Nodes on the top-left
                TopLeftNodes = (np.where(IGAvertices[:,yaxis-1] >= (np.sqrt(3)*IGAvertices[:,xaxis-1] + np.sqrt(3)*x_max - offset))[0]+1).reshape(1,-1)
                meshinpfile.write('*Nset, nset=TopLeftNodes, instance=Part-1-1\n')
                for i in range(0,math.ceil(len(TopLeftNodes[0])/15)): # Abaqus only accepts maximum 15 items per row
                    meshinpfile.write(''.join(','.join(map(str, TopLeftNodes[0][15*i:15*(i+1)])))+',\n')
                # Nodes on the top-right
                TopRightNodes = (np.where(IGAvertices[:,yaxis-1] >= (-np.sqrt(3)*IGAvertices[:,xaxis-1] + np.sqrt(3)*x_max - offset))[0]+1).reshape(1,-1)
                meshinpfile.write('*Nset, nset=TopRightNodes, instance=Part-1-1\n')
                for i in range(0,math.ceil(len(TopRightNodes[0])/15)): # Abaqus only accepts maximum 15 items per row
                    meshinpfile.write(''.join(','.join(map(str, TopRightNodes[0][15*i:15*(i+1)])))+',\n')
                # Nodes on the top
                TopNodes = (np.where(IGAvertices[:,yaxis-1] >= y_max)[0]+1).reshape(1,-1)
                meshinpfile.write('*Nset, nset=TopNodes, instance=Part-1-1\n')
                for i in range(0,math.ceil(len(TopNodes[0])/15)): # Abaqus only accepts maximum 15 items per row
                    meshinpfile.write(''.join(','.join(map(str, TopNodes[0][15*i:15*(i+1)])))+',\n')
                    
                # Nodes on the front
                FrontNodes = (np.where(IGAvertices[:,zaxis-1] >= z_max)[0]+1).reshape(1,-1)
                meshinpfile.write('*Nset, nset=FrontNodes, instance=Part-1-1\n')
                for i in range(0,math.ceil(len(FrontNodes[0])/15)): # Abaqus only accepts maximum 15 items per row
                    meshinpfile.write(''.join(','.join(map(str, FrontNodes[0][15*i:15*(i+1)])))+',\n')
                # Nodes on the back
                BackNodes = (np.where(IGAvertices[:,zaxis-1] <= z_min)[0]+1).reshape(1,-1)
                meshinpfile.write('*Nset, nset=BackNodes, instance=Part-1-1\n')
                for i in range(0,math.ceil(len(BackNodes[0])/15)): # Abaqus only accepts maximum 15 items per row
                    meshinpfile.write(''.join(','.join(map(str, BackNodes[0][15*i:15*(i+1)])))+',\n')
    else:
        if merge_operation in ['on','On','Y','y','Yes','yes']:
            # boundary nodes
            if box_shape == 'square':
                if any(item in boundary_conditions for item in ['left','Left','L']):
                    # Nodes on the left
                    LeftNodes = (np.where(IGAvertices[:,xaxis-1] <= x_min+merge_tol)[0]+1).reshape(1,-1)
                    meshinpfile.write('*Nset, nset=LeftNodes, instance=Part-1-1\n')
                    for i in range(0,math.ceil(len(LeftNodes[0])/15)): # Abaqus only accepts maximum 15 items per row
                        meshinpfile.write(''.join(','.join(map(str, LeftNodes[0][15*i:15*(i+1)])))+',\n')
                if any(item in boundary_conditions for item in ['right','Right','R']):
                    # Nodes on the right
                    RightNodes = (np.where(IGAvertices[:,xaxis-1] >= x_max-merge_tol)[0]+1).reshape(1,-1)
                    meshinpfile.write('*Nset, Nset=RightNodes, instance=Part-1-1\n')
                    for i in range(0,math.ceil(len(RightNodes[0])/15)): # Abaqus only accepts maximum 15 items per row
                        meshinpfile.write(''.join(','.join(map(str, RightNodes[0][15*i:15*(i+1)])))+',\n')
                if any(item in boundary_conditions for item in ['bottom','Bottom','Bot']):
                    # Nodes on the bottom
                    BottomNodes = (np.where(IGAvertices[:,zaxis-1] <= z_min+merge_tol)[0]+1).reshape(1,-1)
                    meshinpfile.write('*Nset, nset=BottomNodes, instance=Part-1-1\n')
                    for i in range(0,math.ceil(len(BottomNodes[0])/15)): # Abaqus only accepts maximum 15 items per row
                        meshinpfile.write(''.join(','.join(map(str, BottomNodes[0][15*i:15*(i+1)])))+',\n')
                if any(item in boundary_conditions for item in ['top','Top','T']):
                    # Nodes on the top
                    TopNodes = (np.where(IGAvertices[:,zaxis-1] >= z_max-merge_tol)[0]+1).reshape(1,-1)
                    meshinpfile.write('*Nset, nset=TopNodes, instance=Part-1-1\n')
                    for i in range(0,math.ceil(len(TopNodes[0])/15)): # Abaqus only accepts maximum 15 items per row
                        meshinpfile.write(''.join(','.join(map(str, TopNodes[0][15*i:15*(i+1)])))+',\n')
                if any(item in boundary_conditions for item in ['back','Back']):
                    # Nodes on the back
                    BackNodes = (np.where(IGAvertices[:,yaxis-1] <= y_min+merge_tol)[0]+1).reshape(1,-1)
                    meshinpfile.write('*Nset, nset=BackNodes, instance=Part-1-1\n')
                    for i in range(0,math.ceil(len(BackNodes[0])/15)): # Abaqus only accepts maximum 15 items per row
                        meshinpfile.write(''.join(','.join(map(str, BackNodes[0][15*i:15*(i+1)])))+',\n')
                if any(item in boundary_conditions for item in ['front','Front','F']):
                    # Nodes on the front
                    FrontNodes = (np.where(IGAvertices[:,yaxis-1] >= y_max-merge_tol)[0]+1).reshape(1,-1)
                    meshinpfile.write('*Nset, nset=FrontNodes, instance=Part-1-1\n')
                    for i in range(0,math.ceil(len(FrontNodes[0])/15)): # Abaqus only accepts maximum 15 items per row
                        meshinpfile.write(''.join(','.join(map(str, FrontNodes[0][15*i:15*(i+1)])))+',\n')
    
            # boundary nodes
            elif box_shape == 'notched_square':
                # if any(item in boundary_conditions for item in ['left','Left','L']):
                    # Nodes on the left bottom
                    LeftBottomNodes = (np.where((IGAvertices[:,xaxis-1] <= x_min+merge_tol) & (IGAvertices[:,yaxis-1] <= 0))[0]+1).reshape(1,-1)
                    meshinpfile.write('*Nset, nset=LeftBottomNodes, instance=Part-1-1\n')
                    for i in range(0,math.ceil(len(LeftBottomNodes[0])/15)): # Abaqus only accepts maximum 15 items per row
                        meshinpfile.write(''.join(','.join(map(str, LeftBottomNodes[0][15*i:15*(i+1)])))+',\n')
                    # Nodes on the left top
                    LeftTopNodes = (np.where((IGAvertices[:,xaxis-1] <= x_min+merge_tol) & (IGAvertices[:,yaxis-1] >= 0))[0]+1).reshape(1,-1)
                    meshinpfile.write('*Nset, nset=LeftTopNodes, instance=Part-1-1\n')
                    for i in range(0,math.ceil(len(LeftTopNodes[0])/15)): # Abaqus only accepts maximum 15 items per row
                        meshinpfile.write(''.join(','.join(map(str, LeftTopNodes[0][15*i:15*(i+1)])))+',\n')
                # if any(item in boundary_conditions for item in ['right','Right','R']):
                    # Nodes on the right
                    RightNodes = (np.where(IGAvertices[:,xaxis-1] >= x_max-merge_tol)[0]+1).reshape(1,-1)
                    meshinpfile.write('*Nset, Nset=RightNodes, instance=Part-1-1\n')
                    for i in range(0,math.ceil(len(RightNodes[0])/15)): # Abaqus only accepts maximum 15 items per row
                        meshinpfile.write(''.join(','.join(map(str, RightNodes[0][15*i:15*(i+1)])))+',\n')
                # if any(item in boundary_conditions for item in ['bottom','Bottom','Bot']):
                    # Nodes on the bottom
                    BottomNodes = (np.where(IGAvertices[:,zaxis-1] <= z_min+merge_tol)[0]+1).reshape(1,-1)
                    meshinpfile.write('*Nset, nset=BottomNodes, instance=Part-1-1\n')
                    for i in range(0,math.ceil(len(BottomNodes[0])/15)): # Abaqus only accepts maximum 15 items per row
                        meshinpfile.write(''.join(','.join(map(str, BottomNodes[0][15*i:15*(i+1)])))+',\n')
                # if any(item in boundary_conditions for item in ['top','Top','T']):
                    # Nodes on the top
                    TopNodes = (np.where(IGAvertices[:,zaxis-1] >= z_max-merge_tol)[0]+1).reshape(1,-1)
                    meshinpfile.write('*Nset, nset=TopNodes, instance=Part-1-1\n')
                    for i in range(0,math.ceil(len(TopNodes[0])/15)): # Abaqus only accepts maximum 15 items per row
                        meshinpfile.write(''.join(','.join(map(str, TopNodes[0][15*i:15*(i+1)])))+',\n')
                # if any(item in boundary_conditions for item in ['back','Back']):
                    # Nodes on the back
                    BackNodes = (np.where(IGAvertices[:,yaxis-1] <= y_min+merge_tol)[0]+1).reshape(1,-1)
                    meshinpfile.write('*Nset, nset=BackNodes, instance=Part-1-1\n')
                    for i in range(0,math.ceil(len(BackNodes[0])/15)): # Abaqus only accepts maximum 15 items per row
                        meshinpfile.write(''.join(','.join(map(str, BackNodes[0][15*i:15*(i+1)])))+',\n')
                # if any(item in boundary_conditions for item in ['front','Front','F']):
                    # Nodes on the front
                    FrontNodes = (np.where(IGAvertices[:,yaxis-1] >= y_max-merge_tol)[0]+1).reshape(1,-1)
                    meshinpfile.write('*Nset, nset=FrontNodes, instance=Part-1-1\n')
                    for i in range(0,math.ceil(len(FrontNodes[0])/15)): # Abaqus only accepts maximum 15 items per row
                        meshinpfile.write(''.join(','.join(map(str, FrontNodes[0][15*i:15*(i+1)])))+',\n')
        else: # not merged
            # boundary nodes
            if box_shape == 'square':
                if any(item in boundary_conditions for item in ['left','Left','L']):
                    # Nodes on the left
                    LeftNodes = (np.where(IGAvertices[:,xaxis-1] <= x_min)[0]+1).reshape(1,-1)
                    meshinpfile.write('*Nset, nset=LeftNodes, instance=Part-1-1\n')
                    for i in range(0,math.ceil(len(LeftNodes[0])/15)): # Abaqus only accepts maximum 15 items per row
                        meshinpfile.write(''.join(','.join(map(str, LeftNodes[0][15*i:15*(i+1)])))+',\n')
                if any(item in boundary_conditions for item in ['right','Right','R']):
                    # Nodes on the right
                    RightNodes = (np.where(IGAvertices[:,xaxis-1] >= x_max)[0]+1).reshape(1,-1)
                    meshinpfile.write('*Nset, Nset=RightNodes, instance=Part-1-1\n')
                    for i in range(0,math.ceil(len(RightNodes[0])/15)): # Abaqus only accepts maximum 15 items per row
                        meshinpfile.write(''.join(','.join(map(str, RightNodes[0][15*i:15*(i+1)])))+',\n')
                if any(item in boundary_conditions for item in ['bottom','Bottom','Bot']):
                    # Nodes on the bottom
                    BottomNodes = (np.where(IGAvertices[:,zaxis-1] <= z_min)[0]+1).reshape(1,-1)
                    meshinpfile.write('*Nset, nset=BottomNodes, instance=Part-1-1\n')
                    for i in range(0,math.ceil(len(BottomNodes[0])/15)): # Abaqus only accepts maximum 15 items per row
                        meshinpfile.write(''.join(','.join(map(str, BottomNodes[0][15*i:15*(i+1)])))+',\n')
                if any(item in boundary_conditions for item in ['top','Top','T']):
                    # Nodes on the top
                    TopNodes = (np.where(IGAvertices[:,zaxis-1] >= z_max)[0]+1).reshape(1,-1)
                    meshinpfile.write('*Nset, nset=TopNodes, instance=Part-1-1\n')
                    for i in range(0,math.ceil(len(TopNodes[0])/15)): # Abaqus only accepts maximum 15 items per row
                        meshinpfile.write(''.join(','.join(map(str, TopNodes[0][15*i:15*(i+1)])))+',\n')
                if any(item in boundary_conditions for item in ['back','Back']):
                    # Nodes on the back
                    BackNodes = (np.where(IGAvertices[:,yaxis-1] <= y_min)[0]+1).reshape(1,-1)
                    meshinpfile.write('*Nset, nset=BackNodes, instance=Part-1-1\n')
                    for i in range(0,math.ceil(len(BackNodes[0])/15)): # Abaqus only accepts maximum 15 items per row
                        meshinpfile.write(''.join(','.join(map(str, BackNodes[0][15*i:15*(i+1)])))+',\n')
                if any(item in boundary_conditions for item in ['front','Front','F']):
                    # Nodes on the front
                    FrontNodes = (np.where(IGAvertices[:,yaxis-1] >= y_max)[0]+1).reshape(1,-1)
                    meshinpfile.write('*Nset, nset=FrontNodes, instance=Part-1-1\n')
                    for i in range(0,math.ceil(len(FrontNodes[0])/15)): # Abaqus only accepts maximum 15 items per row
                        meshinpfile.write(''.join(','.join(map(str, FrontNodes[0][15*i:15*(i+1)])))+',\n')
    
            # boundary nodes
            elif box_shape == 'notched_square':
                if any(item in boundary_conditions for item in ['left','Left','L']):
                    # Nodes on the left bottom
                    LeftBottomNodes = (np.where((IGAvertices[:,xaxis-1] <= x_min) & (IGAvertices[:,yaxis-1] <= 0))[0]+1).reshape(1,-1)
                    meshinpfile.write('*Nset, nset=LeftBottomNodes, instance=Part-1-1\n')
                    for i in range(0,math.ceil(len(LeftBottomNodes[0])/15)): # Abaqus only accepts maximum 15 items per row
                        meshinpfile.write(''.join(','.join(map(str, LeftBottomNodes[0][15*i:15*(i+1)])))+',\n')
                    # Nodes on the left top
                    LeftTopNodes = (np.where((IGAvertices[:,xaxis-1] <= x_min) & (IGAvertices[:,yaxis-1] >= 0))[0]+1).reshape(1,-1)
                    meshinpfile.write('*Nset, nset=LeftTopNodes, instance=Part-1-1\n')
                    for i in range(0,math.ceil(len(LeftTopNodes[0])/15)): # Abaqus only accepts maximum 15 items per row
                        meshinpfile.write(''.join(','.join(map(str, LeftTopNodes[0][15*i:15*(i+1)])))+',\n')
                if any(item in boundary_conditions for item in ['right','Right','R']):
                    # Nodes on the right
                    RightNodes = (np.where(IGAvertices[:,xaxis-1] >= x_max)[0]+1).reshape(1,-1)
                    meshinpfile.write('*Nset, Nset=RightNodes, instance=Part-1-1\n')
                    for i in range(0,math.ceil(len(RightNodes[0])/15)): # Abaqus only accepts maximum 15 items per row
                        meshinpfile.write(''.join(','.join(map(str, RightNodes[0][15*i:15*(i+1)])))+',\n')
                if any(item in boundary_conditions for item in ['bottom','Bottom','Bot']):
                    # Nodes on the bottom
                    BottomNodes = (np.where(IGAvertices[:,zaxis-1] <= z_min)[0]+1).reshape(1,-1)
                    meshinpfile.write('*Nset, nset=BottomNodes, instance=Part-1-1\n')
                    for i in range(0,math.ceil(len(BottomNodes[0])/15)): # Abaqus only accepts maximum 15 items per row
                        meshinpfile.write(''.join(','.join(map(str, BottomNodes[0][15*i:15*(i+1)])))+',\n')
                if any(item in boundary_conditions for item in ['top','Top','T']):
                    # Nodes on the top
                    TopNodes = (np.where(IGAvertices[:,zaxis-1] >= z_max)[0]+1).reshape(1,-1)
                    meshinpfile.write('*Nset, nset=TopNodes, instance=Part-1-1\n')
                    for i in range(0,math.ceil(len(TopNodes[0])/15)): # Abaqus only accepts maximum 15 items per row
                        meshinpfile.write(''.join(','.join(map(str, TopNodes[0][15*i:15*(i+1)])))+',\n')
                if any(item in boundary_conditions for item in ['back','Back']):
                    # Nodes on the back
                    BackNodes = (np.where(IGAvertices[:,yaxis-1] <= y_min)[0]+1).reshape(1,-1)
                    meshinpfile.write('*Nset, nset=BackNodes, instance=Part-1-1\n')
                    for i in range(0,math.ceil(len(BackNodes[0])/15)): # Abaqus only accepts maximum 15 items per row
                        meshinpfile.write(''.join(','.join(map(str, BackNodes[0][15*i:15*(i+1)])))+',\n')
                if any(item in boundary_conditions for item in ['front','Front','F']):
                    # Nodes on the front
                    FrontNodes = (np.where(IGAvertices[:,yaxis-1] >= y_max)[0]+1).reshape(1,-1)
                    meshinpfile.write('*Nset, nset=FrontNodes, instance=Part-1-1\n')
                    for i in range(0,math.ceil(len(FrontNodes[0])/15)): # Abaqus only accepts maximum 15 items per row
                        meshinpfile.write(''.join(','.join(map(str, FrontNodes[0][15*i:15*(i+1)])))+',\n')      
    # ELEMENT SETS
    meshinpfile.write('*Elset, elset=AllElles, instance=Part-1-1, generate\n') 
    meshinpfile.write('{:d}, {:d}, {:d} \n'.format(1,nelem_total,1))
    meshinpfile.write('*Elset, elset=AllVisualElle, instance=Part-1-1, generate\n') 
    meshinpfile.write('{:d}, {:d}, {:d} \n'.format(count_offset+1,count-1,1))
    meshinpfile.write('*End Assembly\n')


def ReadSavedSites(radial_growth_rule):
    """Search in the current directory and all directories above it 
    for a file of a particular name.

    Arguments:
    ---------
    filename :: string, the filename to look for.

    Returns
    -------
    pathlib.Path, the location of the first file found or
    None, if none was found
    """
    
    sitesfile = radial_growth_rule
    radiifile = radial_growth_rule.strip().replace("_sites", '_radii')
    
    d = Path(__file__).parent #Path.cwd()
    root = Path(d.root)

    while d != root:
        site_attempt = d / sitesfile
        radii_attempt = d / radiifile
        if site_attempt.exists() and site_attempt.exists():
            print('Sites info from file: {:s} has been loaded.'.format(str(site_attempt)))
            return np.load(site_attempt), np.load(radii_attempt)
        d = d.parent
    
    print('Could not find file: {:s}, please check if the existing site file is under the same directory with the input script.'.format(radial_growth_rule))
    print('Now exitting...')
    exit()
    
    return None


def StlModelFile(geoName):
    
    """Generate the 3D model file (geoName.stl file) for the generated geometry.

    Arguments:
    ---------
    geoName :: string, geometry name.

    Returns
    -------
    None or error message, if generation failed (TBD)
    """
    
    import vtk
    import os
    
    vtufilename = geoName + '_conns_vol'+'.vtu'
    cwd = os.getcwd()
    vtupath = os.path.join(cwd,'meshes',geoName)
    stlfilename = os.path.join(vtupath, geoName+".stl")
    file_name = os.path.join(vtupath,vtufilename)
    
    # create a new 'XML Unstructured Grid Reader'
    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(file_name)
    reader.Update()
    output = reader.GetOutputPort()

    # Extract the outer (polygonal) surface.
    surface = vtk.vtkDataSetSurfaceFilter()
    surface.SetInputConnection(0, output)
    surface.Update()
    sufaceoutput = surface.GetOutputPort()

    # Write the stl file to disk
    stlWriter = vtk.vtkSTLWriter()
    stlWriter.SetFileName(stlfilename)
    stlWriter.SetInputConnection(sufaceoutput)
    stlWriter.Write()


def ModelInfo(box_shape,boundary_points,z_min,z_max,skeleton_density,MeshData):

    """Calculate the material properties of generated geometry

    Arguments:
    ---------

    Returns
    -------
    total mass, bulk volume, bulk density, and porosity of the model
    """
    
    mass = 0
    
    for i in range(0,MeshData.shape[0]):
        mass += skeleton_density*MeshData[i,11]*MeshData[i,12]*np.linalg.norm(MeshData[i,5:8] - MeshData[i,8:11])

    
    hull = ConvexHull(boundary_points)
    
    volume = (z_max-z_min)*hull.volume

    bulk_density = mass/volume
    
    mass_if_all_solid = skeleton_density*volume # if all skeleton phase
    porosity = 1 - mass/mass_if_all_solid
    
    return mass,volume,bulk_density,porosity


def LogFile(geoName,iter_max,r_min,r_max,nrings,width_heart,width_sparse,width_dense,\
        generation_center,box_shape,box_center,box_size,x_min,x_max,y_min,y_max,
        cellsize_sparse,cellsize_dense,cellwallthickness_sparse,cellwallthickness_dense,\
        merge_operation,merge_tol,precrackFlag,precrack_widths,boundaryFlag,\
        segment_length,theta_min,theta_max,z_min,z_max,long_connector_ratio,\
        NURBS_degree,nctrlpt_per_beam,nconnector_t_precrack,nconnector_l_precrack,\
        nParticles,nbeamElem,skeleton_density,mass,volume,density,porosity,\
        stlFlag,inpFlag,inpType,radial_growth_rule,\
        startTime,placementTime,voronoiTime,RebuildvorTime,BeamTime,FileTime):
    
    """Generate the log file (geoName.log file) for the generation procedure.

    Arguments:
    ---------
    geoName :: string, geometry name.
    -

    Returns
    -------
    None or error message, if generation failed (TBD)
    """
    
    # get current local time
    current_time = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    
    # Generate log file
    logfile = open(Path('meshes/' + geoName + '/' + geoName + '-report.log'),'w')                                      
    logfile.write('################################################################################\n')
    logfile.write('##             Log Created with RingsPy Mesh Generation Tool V' + pkg_resources.get_distribution("ringspy").version + '           ##\n')
    logfile.write('##             Files generated at local time: ' + current_time + '             ##\n')
    logfile.write('################################################################################\n')
    logfile.write('\n')
    logfile.write('GEOMETRY:\n') 
    logfile.write('geoName:                                 ' + geoName + '\n')  
    logfile.write('\n')
    logfile.write('GENERATION PARAMETERS:\n')
    logfile.write('iter_max:                                ' + str(iter_max) + '\n')
    logfile.write('radial_growth_rule:                      ' + str(radial_growth_rule) + '\n')
    logfile.write('merge_operation:                         ' + merge_operation + '\n')
    logfile.write('boundaryFlag:                            ' + boundaryFlag + '\n')
    logfile.write('\n')
    logfile.write('GEOMETRY PARAMETERS:\n')                   
    logfile.write('r_min:                                   ' + str(r_min) + '\n')
    logfile.write('r_max:                                   ' + str(r_max) + '\n')
    logfile.write('nrings:                                  ' + str(nrings) + '\n')
    logfile.write('width_heart:                             ' + str(width_heart) + '\n')
    logfile.write('width_sparse:                            ' + str(width_sparse) + '\n')
    logfile.write('width_dense:                             ' + str(width_dense) + '\n')
    logfile.write('generation_center:                       ' + str(generation_center) + '\n')
    logfile.write('cellsize_sparse:                         ' + str(cellsize_sparse) + '\n')
    logfile.write('cellsize_dense:                          ' + str(cellsize_dense) + '\n')
    logfile.write('cellwallthickness_sparse:                ' + str(cellwallthickness_sparse) + '\n')
    logfile.write('cellwallthickness_dense:                 ' + str(cellwallthickness_dense) + '\n')

    if merge_operation in ['on','On','Y','y','Yes','yes']: 
        logfile.write('merge_tol:                               ' + str(merge_tol) + '\n')

    logfile.write('\n')
    logfile.write('CLIPPING PARAMETERS:\n')
    logfile.write('box_shape:                               ' + str(box_shape) + '\n')
    logfile.write('box_center:                              ' + str(box_center) + '\n')
    logfile.write('box_size:                                ' + str(box_size) + '\n')
    logfile.write('precrackFlag:                            ' + precrackFlag + '\n')
    if precrackFlag in ['on','On','Y','y','Yes','yes']:
        logfile.write('\n')
        logfile.write('PRECRACK INSERTION PARAMETERS:\n')
        logfile.write('precrack_widths:                         ' + str(precrack_widths) + '\n')
        logfile.write('nconnector_t_precrack:                   ' + str(nconnector_t_precrack) + '\n')
        logfile.write('nconnector_l_precrack:                   ' + str(nconnector_l_precrack) + '\n')
    logfile.write('\n')
    logfile.write('GRAIN EXTRUSION PARAMETERS:\n')
    logfile.write('segment_length:                          ' + str(segment_length) + '\n')
    logfile.write('theta_min:                               ' + str(theta_min) + '\n')
    logfile.write('theta_max:                               ' + str(theta_max) + '\n')
    logfile.write('long_connector_ratio:                    ' + str(long_connector_ratio) + '\n')
    logfile.write('NURBS_degree:                            ' + str(NURBS_degree) + '\n')
    logfile.write('nctrlpt_per_beam:                        ' + str(nctrlpt_per_beam) + '\n')
    logfile.write('\n')
    logfile.write('GENERATION:\n')  
    logfile.write('nParticles:                              ' + str(nParticles) + '\n')
    logfile.write('nbeamElem:                               ' + str(nbeamElem) + '\n')
    logfile.write('MODEL PROPERTIES:\n')
    logfile.write('skeleton density:                        ' + str(skeleton_density) + '\n')
    logfile.write('total mass:                              ' + str(mass) + '\n')
    logfile.write('bulk volume:                             ' + str(volume) + '\n')
    logfile.write('bulk density:                            ' + str(density) + '\n')
    logfile.write('porosity:                                ' + str(porosity) + '\n')
    logfile.write('\n')
    logfile.write('PERFORMANCE:\n')  
    logfile.write('Placement Time:                          ' + str(placementTime - startTime) + '\n')
    logfile.write('Tessellation Time:                       ' + str(voronoiTime - placementTime) + '\n')
    logfile.write('Data Reconstruction Time:                ' + str(RebuildvorTime - voronoiTime) + '\n')
    logfile.write('Grain Extrusion Time:                    ' + str(BeamTime - RebuildvorTime) + '\n')
    logfile.write('File Writing Time:                       ' + str(FileTime - BeamTime) + '\n')
    logfile.write('Total Time:                              ' + str(FileTime - startTime) + '\n')
    logfile.write('\n')
    logfile.write('FILE CREATION:\n')  
    logfile.write('Cross-sectional Image File:              ./meshes/' + geoName + '/' + geoName \
        + '.png\n')
    logfile.write('NURBS Beam File:                         ./meshes/' + geoName + '/' + geoName \
        + 'IGA.txt\n')
    logfile.write('Connector Data File:                     ./meshes/' + geoName + '/' + geoName \
        + '-mesh.txt\n')
    logfile.write('Grain-ridge Data File:                   ./meshes/' + geoName + '/' + geoName \
        + '-vertex.mesh\n')
    logfile.write('Ridge Data File:                         ./meshes/' + geoName + '/' + geoName \
        + '-ridge.mesh\n')
    logfile.write('Paraview Vertices File:                  ./meshes/' + geoName + '/' + geoName
                  
        + '_vertices'+'.vtu\n')
    logfile.write('Paraview Beams File:                     ./meshes/' + geoName + '/' \
        + geoName + '_beams'+'.vtu\n')
    logfile.write('Paraview Connectors File:                ./meshes/' + geoName + '/' + geoName \
        + '_conns'+'.vtu\n')
    logfile.write('Paraview Connectors (Volume) File:       ./meshes/' + geoName + '/' + geoName \
        + '_conns_vol'+'.vtu\n')
    if stlFlag in ['on','On','Y','y','Yes','yes']:
        logfile.write('3D Model File:                           ./meshes/' + geoName + '/' + geoName \
            + '.stl\n')
    if inpFlag in ['on','On','Y','y','Yes','yes']:
        if inpType in ['abaqus','Abaqus','ABQ','abq','ABAQUS','Abq']:
            logfile.write('Abaqus Input File:                       ./meshes/' + geoName + '/' + geoName \
            + '.inp\n') 
    logfile.write('Log File:                                ./meshes/' + geoName + '/' + geoName \
        + '-report.log\n')
    logfile.close()