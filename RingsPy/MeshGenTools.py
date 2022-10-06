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
from scipy.spatial import cKDTree
from pathlib import Path
import datetime


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

def check_isinside_boundcircle(circC,log_center,rmin,rmax):
    ''' Make sure this circle does not protrude outside the ring's bounding circle '''
    x = circC[0]
    y = circC[1]
    w = circC[2]
    within = ( math.sqrt((x-log_center[0])**2 + (y-log_center[1])**2) > rmin+w ) and ( math.sqrt((x-log_center[0])**2 + (y-log_center[1])**2)<rmax-w )
    
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
    ''' Make sure this circle does not protrude outside the bounding box '''
    x = circC[0]
    y = circC[1]
    w = circC[2]
    within = ( x > x_min+w ) and ( x < x_max-w ) and ( y > y_min+w ) and ( y < y_max-w )
    
    return within

def check_isinside_boundbox2Dindent(circC,x_min,x_max,y_min,y_max,x_indent,y_indent_min,y_indent_max):
    ''' Make sure this circle does not protrude outside the indent bounding box '''
    x = circC[0]
    y = circC[1]
    w = circC[2]
    within = ( x > x_min+w ) and ( x < x_max-w ) and ( y > y_min+w ) and ( y < y_max-w )
    within_indent = ( x > x_min+w ) and ( x < x_indent-w ) and ( y > y_indent_min+w ) and ( y < y_indent_max-w )
    within = np.subtract(within,within_indent,dtype=np.float32)
    return within

def check_isinside_boundbox2Dprecrack(circC,x_min,x_max,y_min,y_max,x_indent,y_indent,x_precrack,y_precrack):
    ''' Make sure this circle does not protrude outside the ring's bounding box '''
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
    ''' Make sure this circle does not protrude outside the bounding hexagon 
    credit to: http://www.playchilla.com/how-to-check-if-a-point-is-inside-a-hexagon
    '''
    x = np.abs(circC[0] - box_center[0])
    y = np.abs(circC[1] - box_center[1])
    v = L/2
    h = L/np.sqrt(3)/2
    within = ( x < 2*h ) and ( y < v ) and ((2*v*h - x*v - h*y) > 0)
    return within

def check_overlap(circles,circC):
    ''' Make sure the distance between the current circle's center and all
        other circle centers is greater than or equal to the circle's perimeter (2r)
    '''
    x = circC[0]
    y = circC[1]
    w = circC[2]
    nooverlap = all( (x - c[0])**2 + (y - c[1])**2 >= (w + c[2])**2 for c in circles )
    return nooverlap

def check_iscollinear(p1,p2,boundaries):
    ''' Make sure the line segment x1-x2 is collinear with 
        any line segment of the boundary polygon
    '''
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


def CellPlacement_Wood(log_center,r_max,r_min,nrings,width_heart,
                       width_early,width_late,cellsize_early,cellsize_late,\
                       iter_max,print_interval):
    """
    packing cells in annual rings
    """
    # Annual ring width distribution
    width = np.concatenate(([width_heart],np.tile([width_early,width_late],nrings)))
    noise = np.random.normal(1,0.25,len(width))
    width = np.multiply(width,noise)
    radius = np.concatenate(([0],np.cumsum(width)))
    
    # Generate circle seeds in each annual ring (can be parallelized in the future)
    # circles: list = list()
    circles = []
    for ibin in range(0,nrings*2+1):
        iter = 0 
        while iter < iter_max:
            r = (radius[ibin+1]-radius[ibin])*math.sqrt(np.random.random()) + radius[ibin]  # randomly generate a point in the domain
            t = np.random.random()*2*math.pi 
            
            x = r*math.cos(t) + log_center[0]
            y = r*math.sin(t) + log_center[1]
    
            if (ibin % 2) == 0: # if even, latewood
                w = cellsize_late #scipy.stats.truncnorm.rvs((we_lower - we_mu) / we_sigma, (we_upper - we_mu) / we_sigma, loc=we_mu, scale=we_sigma)
            else:
                w = cellsize_early #scipy.stats.truncnorm.rvs((wl_lower - wl_mu) / wl_sigma, (wl_upper - wl_mu) / wl_sigma, loc=wl_mu, scale=wl_sigma)
    
            circC = [x, y, w]
            
            if check_overlap(circles, circC) and check_isinside_boundcircle(circC,log_center,radius[ibin],radius[ibin+1]):
                circles.append(circC)
                
                # regularly print 
                if (len(circles) % print_interval == 0):
                    if (nrings*2+1-ibin == 1):
                        print('{:d} ring remaining; {:d} cells/particles placed.'.format(nrings*2+1-ibin,len(circles)))
                    else:
                        print('{:d} rings remaining; {:d} cells/particles placed.'.format(nrings*2+1-ibin,len(circles)))
            iter += 1 
        
    sites = np.array(circles)

    # take off out of bound Voronoi vertices
    outofbound = []
    for i in range(0,sites.shape[0]):
        if( ((sites[i,0]-log_center[0])**2+(sites[i,1]-log_center[1])**2) > (1.2*r_max)**2 ):
            outofbound.append(i)
    sites = np.delete(sites,outofbound,0)

    return sites, radius


def CellPlacement_Honeycomb(log_center,r_max,r_min,nrings,box_center,box_size,\
                            width_heart,width_early,width_late,\
                            cellsize_early,cellsize_late,\
                            cellwallthickness_early,cellwallthickness_late,\
                            cellangle,iter_max,print_interval):
    from hexalattice.hexalattice import create_hex_grid
    """
    packing cells in annual rings
    """
    # Annual ring width distribution
    width = np.concatenate(([width_heart],np.tile([width_early,width_late],nrings)))
    noise = np.random.normal(1,0.25,len(width))
    width = np.multiply(width,noise)
    radius = np.concatenate(([0],np.cumsum(width)))
    
    cellsize = (cellsize_early+cellsize_late)/2
    nx = int(2*r_max/cellsize)
    ny = int(2*r_max/cellsize)
    # The hexagonal lattice sites are generated via hexalattice: https://pypi.org/project/hexalattice/
    hex_centers, _ = create_hex_grid(nx=nx,
                                     ny=ny,
                                     min_diam=cellsize,
                                     rotate_deg=cellangle,
                                     do_plot=False)
    sites = np.hstack((hex_centers+log_center,np.ones([hex_centers.shape[0],1])))
    sites = np.asarray(sites)

    # take off out of bound Voronoi vertices
    outofbound = []
    for i in range(0,sites.shape[0]):
        if( ((sites[i,0]-log_center[0])**2+(sites[i,1]-log_center[1])**2) > (1.2*r_max)**2 ):
            outofbound.append(i)
    sites = np.delete(sites,outofbound,0)

    return sites, radius


def RebuildVoronoi(vor,circles,boundaries,log_center,x_min,x_max,y_min,y_max,box_center,boundaryFlag):
    '''Clip Voronoi mesh by the boundaries, rebuild the new Voronoi mesh'''
    # Store indices of Voronoi vertices for each finite ridge
    finite_ridges = []
    finite_ridges_pointid = []
    infinite_ridges = []
    infinite_ridges_pointid = []
    boundary_points = []
    
    # Form two new Voronoi vertices lists with and without out-of-bound vertices
    voronoi_vertices_in = []
    voronoi_vertices_out = []
    
    if len(boundaries) == 3:
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
    elif len(boundaries) == 4:
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
    elif len(boundaries) == 6:
        L = y_max - y_min
        box_center = log_center
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
    elif len(boundaries) == 8: # grooved
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
                    if np.sign(np.dot(midpoint - log_center, normal)) == 1: # facing outwards
                        if (t >= 0) and math.isfinite(t):
                            t_final = t
                    elif np.sign(np.dot(midpoint - log_center, normal)) == -1: # facing inwards
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

def RebuildVoronoi_new(vor,circles,boundaries,log_center,x_min,x_max,y_min,y_max,box_center,merge_tol,boundaryFlag):
    '''Clip Voronoi mesh by the boundaries, rebuild the new Voronoi mesh'''
    # Store indices of Voronoi vertices for each finite ridge
    finite_ridges = []
    finite_ridges_pointid = []
    infinite_ridges = []
    infinite_ridges_pointid = []
    boundary_points = []
    
    # Form two new Voronoi vertices lists with and without out-of-bound vertices
    voronoi_vertices_in = []
    voronoi_vertices_out = []
    
    if len(boundaries) == 8: # grooved
        
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
                    if np.sign(np.dot(midpoint - log_center, normal)) == 1: # facing outwards
                        if (t >= 0) and math.isfinite(t):
                            t_final = t
                    elif np.sign(np.dot(midpoint - log_center, normal)) == -1: # facing inwards
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
                       npt_per_layer_vtk,nridge,geoName,radius,log_center,\
                       cellwallthickness_early,cellwallthickness_late):
    
    '''Generate the vertex and ridge info 
       Vertex info includes: coordinates, ridge indices, indices of another vertex for each ridge, ridge lengths, ridge angles, ridge width
       Ridge info includes: 2 vertex indices, 2 mid point indices, 2 quarter point indices '''

    
    # Calculate lengths of all Voronoi ridges
    vector = all_pts_2D[all_ridges[:,1],:] - all_pts_2D[all_ridges[:,0],:]
    all_ridge_lengths = np.linalg.norm(vector, axis=1)
    
    # Calculate angles of all Voronoi ridges (angles measured counter-clock wise, x-axis --> (1,0), y-axis --> (0,1))
    all_ridge_angles = np.arctan2(vector[:,1],vector[:,0]) # np.arctan2(y, x) * 180 / np.pi = the angle

#======== This part can be expanded to allow the smooth transition of cell wall thickness =========
    # Case 1: Abrupt transition of cell wall thickness between earlywood/latewood
    # thicknesses of all Voronoi vertices (here identical thickness for all wings belonging to the same vertex is assumed)
    vertex_cellwallthickness_2D = np.zeros(npt_per_layer)
    
    vertex_distance2logcenter = np.linalg.norm(all_pts_2D - log_center, axis=1) # calc distance between vertex and logcenter
    for i in range(0,npt_per_layer):
        if (bisect.bisect(radius,vertex_distance2logcenter[i]) % 2) != 0: # if odd, earlywood
            vertex_cellwallthickness_2D[i] = cellwallthickness_early
        else: # if even, latewood
            vertex_cellwallthickness_2D[i] = cellwallthickness_late
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
        ,header='Vertex Data Generated with Wood Mesh Generation Tool\n\
Number of vertices\n'+ str(npt_per_layer) +
    '\n\
Max number of wings for one vertex\n'+ str(max_wings) +
    '\n\
[xcoord ycoord nwings ridge1 ... ridgen farvertex1 ... farvertexn length1 ... lengthn width1 ... widthn angle1 ... anglen]', comments='')
    
    np.savetxt(Path('meshes/' + geoName + '/' + geoName +'-ridge.mesh'), all_ridges, fmt='%d', delimiter=' '\
        ,header='Ridge Data Generated with Wood Mesh Generation Tool\n\
Number of ridges\n'+ str(nridge) +
    '\n\
[vertex1 vertex2 midpt1 midpt2 qrtrpt1 qrtrpt2]', comments='')
    
    return all_vertices_2D, max_wings, flattened_all_vertices_2D, all_ridges


def GenerateBeamElement(NURBS_degree,nctrlpt_per_beam,fiberlength,theta_min,theta_max,\
                        z_min,z_max,long_connector_ratio,npt_per_layer,voronoi_vertices,\
                        nvertex,voronoi_ridges,nridge,log_center,all_vertices_2D,max_wings,\
                        flattened_all_vertices_2D,all_ridges):
    

    nctrlpt_per_elem = NURBS_degree + 1
    nbeam_per_grain = int(round((z_max-z_min)/fiberlength))
    
    nlayer = nctrlpt_per_beam*nbeam_per_grain
    
    nconnector_t_per_beam = int((nctrlpt_per_beam-1)/NURBS_degree+1)
    nconnector_t_per_grain = int(nconnector_t_per_beam*nbeam_per_grain)
    
    theta = np.linspace(theta_min,theta_max,nlayer-(nbeam_per_grain-1))
    z_coord = np.linspace(z_min,z_max,nlayer-(nbeam_per_grain-1))
    connector_l_length = fiberlength*long_connector_ratio
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
            vertices_new[i,j,:2] = rotate_around_point_highperf(voronoi_vertices[j,:], theta[i], log_center)
            vertices_new[i,j,2] = z_coord[i]
    
    # Vertex Data for IGA
    woodIGAvertices = np.reshape(vertices_new,(-1,3))
    
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
    
    return woodIGAvertices,vertex_connectivity,beam_connectivity_original,nbeam_total,\
        beam_connectivity,nbeamElem,nlayer,connector_t_connectivity,\
        connector_t_bot_connectivity,connector_t_top_connectivity,\
        connector_t_reg_connectivity,connector_l_connectivity,nconnector_t_per_beam,\
        nconnector_t_per_grain,nconnector_t,nconnector_l,nconnector_total,\
        theta,z_coord,nbeam_per_grain,connector_l_vertex_dict


def ConnectorMeshFile(geoName,woodIGAvertices,connector_t_bot_connectivity,\
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
        Meshdata[i,0:3] = np.copy(woodIGAvertices)[connector_t_bot_connectivity[i,0]-1,:]
        Meshdata[i,3:6] = np.copy(woodIGAvertices)[connector_t_bot_connectivity[i,1]-1,:]
        Meshdata[i,16] = height_connector_t
    for i in range(0,nelem_connector_t_reg):
        Meshdata[i+nelem_connector_t_bot,0:3] = np.copy(woodIGAvertices)[connector_t_reg_connectivity[i,0]-1,:]
        Meshdata[i+nelem_connector_t_bot,3:6] = np.copy(woodIGAvertices)[connector_t_reg_connectivity[i,1]-1,:]
        Meshdata[i+nelem_connector_t_bot,16] = height_connector_t*2
    for i in range(0,nelem_connector_t_top):
        Meshdata[i+nelem_connector_t_bot+nelem_connector_t_reg,0:3] = np.copy(woodIGAvertices)[connector_t_top_connectivity[i,0]-1,:]
        Meshdata[i+nelem_connector_t_bot+nelem_connector_t_reg,3:6] = np.copy(woodIGAvertices)[connector_t_top_connectivity[i,1]-1,:]
        Meshdata[i+nelem_connector_t_bot+nelem_connector_t_reg,16] = height_connector_t
        
    for i in range(0,nelem_connector_l):
        Meshdata[i+nelem_connector_t_bot+nelem_connector_t_reg+nelem_connector_t_top,0:3] = np.copy(woodIGAvertices)[connector_l_connectivity[i,0]-1,:]
        Meshdata[i+nelem_connector_t_bot+nelem_connector_t_reg+nelem_connector_t_top,3:6] = np.copy(woodIGAvertices)[connector_l_connectivity[i,1]-1,:]
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
    ,header='# Connector Data Generated with Wood Mesh Generation Tool\n\
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
                     cellsize_early,):
    
    precrack_midpts = (precrack_nodes[:,0:2]+precrack_nodes[:,2:4])/2.0
    ridge_midpts = all_pts_2D[all_ridges[:,2]]
    ridge_midpts_tree = cKDTree(ridge_midpts)
    near_ridges = []
    for midpt in precrack_midpts:
        near_ridges.append(ridge_midpts_tree.query_ball_point(midpt,3*cellsize_early))
        
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
                       fiberlength,theta,z_coord,nbeam_per_grain,nridge,\
                       voronoi_ridges,log_center,all_ridges,nvertex,\
                       nconnector_t,nconnector_l,nctrlpt_per_beam,ConnMeshData,\
                       all_vertices_2D,max_wings,flattened_all_vertices_2D):
    ngrain = nvertex
    ninterval_per_beam_vtk = int((nctrlpt_per_beam-1)/2) # 2 layers ----> 1 interval

    nconnector_t_per_beam = int((nctrlpt_per_beam-1)/NURBS_degree+1)
    
    # Data for VTK use
    vertices_new = np.zeros((nlayer,npt_per_layer_vtk,3))
    for i in range(0,nlayer):
        for j in range(0,npt_per_layer_vtk):
            vertices_new[i,j,:2] = rotate_around_point_highperf(all_pts_2D[j,:], theta[i], log_center)
            vertices_new[i,j,2] = z_coord[i]
    
    # Vertex Data for VTK use
    woodVTKvertices = np.reshape(vertices_new,(-1,3))
    
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
    woodVTKvertices_quad_conn = np.vstack((woodVTKvertices,quad_conn_vertices))
    npt_total_vtk_quad_conn = woodVTKvertices_quad_conn.shape[0]
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
    woodVTKvertices_hex_conn = np.vstack((woodVTKvertices,hex_conn_vertices))
    npt_total_vtk_hex_conn = woodVTKvertices_hex_conn.shape[0]
    vertex_connectivity_vtk_hex_conn = np.linspace(0,npt_total_vtk_hex_conn-1,npt_total_vtk_hex_conn)
    
# =============================================================================
    # Paraview Visualization File
    collocation_flag_vtk = np.concatenate((np.ones(ngrain),np.zeros(npt_per_layer_vtk*NURBS_degree-ngrain)))
    collocation_flag_vtk = np.concatenate((np.tile(collocation_flag_vtk, nconnector_t_per_beam-1),np.concatenate((np.ones(ngrain),np.zeros(npt_per_layer_vtk-ngrain)))))
    collocation_flag_vtk = np.tile(collocation_flag_vtk, nbeam_per_grain)

# =============================================================================
    # Paraview Vertices File
    woodVTKcell_types_vertices = (np.ones(npt_total_vtk)).astype(int)

    ncell_vertices = woodVTKcell_types_vertices.shape[0]

    vtkfile_vertices = open (Path('meshes/' + geoName + '/' + geoName + '_vertices'+'.vtu'),'w')
    
    vtkfile_vertices.write('<VTKFile type="UnstructuredGrid" version="2.0" byte_order="LittleEndian">'+'\n')
    vtkfile_vertices.write('<UnstructuredGrid>'+'\n')
    vtkfile_vertices.write('<Piece NumberOfPoints="'+str(npt_total_vtk)+'"'+' '+'NumberOfCells="'+str(ncell_vertices)+'">'+'\n')
    
    # <Points>
    vtkfile_vertices.write('<Points>'+'\n')
    vtkfile_vertices.write('<DataArray type="Float64" NumberOfComponents="3" format="ascii">'+'\n')
    for i in range(0,npt_total_vtk):
        X,Y,Z = woodVTKvertices[i]
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
        element = woodVTKcell_types_vertices[i]
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
    woodVTKcell_types_beams = np.concatenate((21*np.ones(nbeam_total_vtk),23*np.ones(nquad_total))).astype(int)
  
    ncell_beams = woodVTKcell_types_beams.shape[0]

    vtkfile_beams = open (Path('meshes/' + geoName + '/' + geoName + '_beams'+'.vtu'),'w')
    
    vtkfile_beams.write('<VTKFile type="UnstructuredGrid" version="2.0" byte_order="LittleEndian">'+'\n')
    vtkfile_beams.write('<UnstructuredGrid>'+'\n')
    vtkfile_beams.write('<Piece NumberOfPoints="'+str(npt_total_vtk)+'"'+' '+'NumberOfCells="'+str(ncell_beams)+'">'+'\n')
    
    # <Points>
    vtkfile_beams.write('<Points>'+'\n')
    vtkfile_beams.write('<DataArray type="Float64" NumberOfComponents="3" format="ascii">'+'\n')
    for i in range(0,npt_total_vtk):
        X,Y,Z = woodVTKvertices[i]
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
        element = woodVTKcell_types_beams[i]
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
    woodVTKcell_types_conns = np.concatenate((3*np.ones(nconnector_t),3*np.ones(nconnector_l),9*np.ones(nconnector_t),9*np.ones(nconnector_l))).astype(int)
    ncell_conns = woodVTKcell_types_conns.shape[0]
    
    Quad_width_vtk = np.tile(Quad_width,(2))

    vtkfile_conns = open (Path('meshes/' + geoName + '/' + geoName + '_conns'+'.vtu'),'w')
    
    vtkfile_conns.write('<VTKFile type="UnstructuredGrid" version="2.0" byte_order="LittleEndian">'+'\n')
    vtkfile_conns.write('<UnstructuredGrid>'+'\n')
    vtkfile_conns.write('<Piece NumberOfPoints="'+str(npt_total_vtk_quad_conn)+'"'+' '+'NumberOfCells="'+str(ncell_conns)+'">'+'\n')
    
    # <Points>
    vtkfile_conns.write('<Points>'+'\n')
    vtkfile_conns.write('<DataArray type="Float64" NumberOfComponents="3" format="ascii">'+'\n')
    for i in range(0,npt_total_vtk_quad_conn):
        X,Y,Z = woodVTKvertices_quad_conn[i]
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
        element = woodVTKcell_types_conns[i]
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
    woodVTKcell_types_conns_vol = np.concatenate((12*np.ones(nconnector_t),12*np.ones(nconnector_l),12*np.ones(nconnector_t),12*np.ones(nconnector_l))).astype(int)
    ncell_conns_vol = woodVTKcell_types_conns_vol.shape[0]
    
    Quad_width_vtk = np.tile(Quad_width,(2))

    vtkfile_conns_vol = open (Path('meshes/' + geoName + '/' + geoName + '_conns_vol'+'.vtu'),'w')
    
    vtkfile_conns_vol.write('<VTKFile type="UnstructuredGrid" version="2.0" byte_order="LittleEndian">'+'\n')
    vtkfile_conns_vol.write('<UnstructuredGrid>'+'\n')
    vtkfile_conns_vol.write('<Piece NumberOfPoints="'+str(npt_total_vtk_hex_conn)+'"'+' '+'NumberOfCells="'+str(ncell_conns_vol)+'">'+'\n')
    
    # <Points>
    vtkfile_conns_vol.write('<Points>'+'\n')
    vtkfile_conns_vol.write('<DataArray type="Float64" NumberOfComponents="3" format="ascii">'+'\n')
    for i in range(0,npt_total_vtk_hex_conn):
        X,Y,Z = woodVTKvertices_hex_conn[i]
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
        element = woodVTKcell_types_conns_vol[i]
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
    
    
def ModelInfo(nbeam_per_grain,fiberlength,NURBS_degree,nctrlpt_per_beam,\
            props_beam,props_connector_t_bot,props_connector_t_reg,\
            props_connector_t_top,props_connector_l,iprops_beam,\
            beam_connectivity,connector_t_bot_connectivity,\
            connector_t_reg_connectivity,connector_t_top_connectivity,\
            connector_l_connectivity,MeshData):
    
    beam_length = fiberlength/nbeam_per_grain
    nelem_beam = beam_connectivity.shape[0]
    nelem_connector_t_bot = connector_t_bot_connectivity.shape[0]
    nelem_connector_t_reg = connector_t_reg_connectivity.shape[0]
    nelem_connector_t_top = connector_t_top_connectivity.shape[0]
    nelem_connector_l = connector_l_connectivity.shape[0]
    
    density_beam = props_beam[0]
    density_connector_t_bot = props_connector_t_bot[0]
    density_connector_t_top = props_connector_t_top[0]
    density_connector_t_reg = props_connector_t_reg[0]
    density_connector_l = props_connector_l[0]
    mass = 0
    
    # mass of beams
    # for i in range(0,nelem_beam):
    mass += nelem_beam*density_beam*props_beam[3]*props_beam[4]*beam_length
    # mass of bot transverse connector
    for i in range(0,nelem_connector_t_bot):
        mass += density_connector_t_bot*props_connector_t_bot[3]*MeshData[i,11]*np.linalg.norm(MeshData[i,5:8] - MeshData[i,8:11])
    # mass of reg transverse connector
    for i in range(nelem_connector_t_bot,nelem_connector_t_bot+nelem_connector_t_reg):
        mass += density_connector_t_reg*props_connector_t_reg[3]*MeshData[i,11]*np.linalg.norm(MeshData[i,5:8] - MeshData[i,8:11])
    # mass of top transverse connector
    for i in range(nelem_connector_t_bot+nelem_connector_t_reg,nelem_connector_t_bot+nelem_connector_t_reg+nelem_connector_t_top):
        mass += density_connector_t_top*props_connector_t_top[3]*MeshData[i,11]*np.linalg.norm(MeshData[i,5:8] - MeshData[i,8:11])
    # mass of longitudinal connector
    for i in range(nelem_connector_t_bot+nelem_connector_t_reg+nelem_connector_t_top,MeshData.shape[0]):
        mass += density_connector_l*props_connector_l[3]*np.linalg.norm(MeshData[i,5:8] - MeshData[i,8:11])    
    return mass


def StlModelFile(geoName):
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


def LogFile(geoName,iter_max,r_min,r_max,nrings,width_heart,width_early,width_late,\
        log_center,box_shape,box_center,box_size,x_min,x_max,y_min,y_max,
        cellsize_early,cellsize_late,cellwallthickness_early,cellwallthickness_late,\
        merge_operation,merge_tol,precrackFlag,precrack_widths,boundaryFlag,\
        fiberlength,theta_min,theta_max,z_min,z_max,long_connector_ratio,\
        NURBS_degree,nctrlpt_per_beam,nconnector_t_precrack,nconnector_l_precrack,\
        nParticles,nbeamElem,\
        STLFlag,\
        startTime,placementTime,voronoiTime,RebuildvorTime,BeamTime,FileTime):

    # get current local time
    current_time = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    
    growth_rule = 'binary'
    
    # Generate log file
    logfile = open(Path('meshes/' + geoName + '/' + geoName + '-report.log'),'w')                                      
    logfile.write('################################################################################\n')
    logfile.write('##             Log Created with RingsPy Mesh Generation Tool V0.1.1           ##\n')
    logfile.write('##             Files generated at local time: ' + current_time + '             ##\n')
    logfile.write('################################################################################\n')
    logfile.write('\n')
    logfile.write('GEOMETRY:\n') 
    logfile.write('geoName:                                 ' + geoName + '\n')  
    logfile.write('\n')
    logfile.write('GENERATION PARAMETERS:\n')
    logfile.write('iter_max:                                ' + str(iter_max) + '\n')
    logfile.write('radial_growth_rule:                      ' + str(growth_rule) + '\n')
    logfile.write('merge_operation:                         ' + merge_operation + '\n')
    logfile.write('boundaryFlag:                            ' + boundaryFlag + '\n')
    logfile.write('\n')
    logfile.write('GEOMETRY PARAMETERS:\n')                   
    logfile.write('r_min:                                   ' + str(r_min) + '\n')
    logfile.write('r_max:                                   ' + str(r_max) + '\n')
    logfile.write('nrings:                                  ' + str(nrings) + '\n')
    logfile.write('width_heart:                             ' + str(width_heart) + '\n')
    logfile.write('width_early:                             ' + str(width_early) + '\n')
    logfile.write('width_late:                              ' + str(width_late) + '\n')
    logfile.write('log_center:                              ' + str(log_center) + '\n')
    logfile.write('cellsize_early:                          ' + str(cellsize_early) + '\n')
    logfile.write('cellsize_late:                           ' + str(cellsize_late) + '\n')
    logfile.write('cellwallthickness_early:                 ' + str(cellwallthickness_early) + '\n')
    logfile.write('cellwallthickness_late:                  ' + str(cellwallthickness_late) + '\n')

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
    logfile.write('fiberlength:                             ' + str(fiberlength) + '\n')
    logfile.write('theta_min:                               ' + str(theta_min) + '\n')
    logfile.write('theta_max:                               ' + str(theta_max) + '\n')
    logfile.write('long_connector_ratio:                    ' + str(long_connector_ratio) + '\n')
    logfile.write('NURBS_degree:                            ' + str(NURBS_degree) + '\n')
    logfile.write('nctrlpt_per_beam:                        ' + str(nctrlpt_per_beam) + '\n')
    logfile.write('\n')
    logfile.write('GENERATION:\n')  
    logfile.write('nParticles:                              ' + str(nParticles) + '\n')
    logfile.write('nbeamElem:                               ' + str(nbeamElem) + '\n')
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
    # logfile.write('Geometry File:                           ./meshes/' + geoName + '/' + geoFile + '\n')
    logfile.write('NURBS beam File:                         ./meshes/' + geoName + '/' + geoName \
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
    if STLFlag in ['on','On','Y','y','Yes','yes']:
        logfile.write('3D Model File:                           ./meshes/' + geoName + '/' + geoName \
            + '.stl\n')   
    logfile.write('Log File:                                ./meshes/' + geoName + '/' + geoName \
        + '-report.log\n')
    logfile.close()