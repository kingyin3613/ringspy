# -*- coding: utf-8 -*-
"""
RingsPy utilities and functions

author: king_yin3613
email: haoyin2022@u.northwestern.edu
"""

import math
import bisect
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from scipy.spatial import cKDTree
from scipy.spatial import ConvexHull
from scipy.spatial import Voronoi
from pathlib import Path
import datetime

from . import __version__

# Accepted spellings of the clipping box shapes
SQUARE_SHAPES = ['square','Square','SQUARE','rectangle','Rectangle','RECTANGLE','cube','Cube','CUBE','s']
TRIANGLE_SHAPES = ['triangle','Triangle','TRIANGLE','triangular','Triangular','TRIANGULAR','tri','Tri','TRI','^']
HEXAGON_SHAPES = ['hexagon','Hexagon','HEXAGON','hexagonal','Hexagonal','HEXAGONAL','hex','Hex','HEX','h']
POLYGON_SHAPES = ['polygon','Polygon','POLYGON','poly']   # any simple polygon, given by its vertices
ON_FLAGS = ['on','On','Y','y','Yes','yes']
ABAQUS_TYPES = ['abaqus','Abaqus','ABQ','abq','ABAQUS','Abq']


def _is_on(flag):
    """Whether an on/off flag is on: a boolean, or one of the ON_FLAGS strings"""
    if isinstance(flag, (bool, np.bool_)):
        return bool(flag)
    return flag in ON_FLAGS


def _flag_text(flag):
    """Flag as written in the report log: 'on'/'off' for booleans, strings unchanged"""
    if isinstance(flag, (bool, np.bool_)):
        return 'on' if flag else 'off'
    return str(flag)


def _rng(rng):
    """Random number source for the cell placement functions.
    None: the global NumPy random state (np.random.seed controls it); an int: a seed for
    np.random.default_rng; otherwise a np.random.Generator (or RandomState) used as is."""
    if rng is None:
        return np.random
    if isinstance(rng, (int, np.integer)) and not isinstance(rng, bool):
        return np.random.default_rng(rng)
    return rng


def cross2D(a, b):
    """
    z-component of the cross product of 2D vectors (or arrays of 2D vectors).
    NumPy 2 no longer accepts 2D vectors in np.cross.
    """
    a = np.asarray(a)
    b = np.asarray(b)
    return a[...,0]*b[...,1] - a[...,1]*b[...,0]


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
    A function to check whether point P(x, y) lies strictly inside the
    triangle formed by the start points of the three boundary lines
    (P must be on the same side of all three edges)
    """
    p = np.asarray(circC[0:2])
    corners = np.array([boundary[1][0] for boundary in boundaries[0:3]])
    edges = np.roll(corners,-1,axis=0) - corners
    sides = cross2D(edges, p - corners)
    within = bool(np.all(sides > 0) or np.all(sides < 0))

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

class _CircleGrid:
    """
    Placed circles binned in a uniform grid, so that an overlap check only
    looks at the circles in the 3x3 neighbouring grid cells. The grid spacing
    is 1.5 times the largest possible sum of two radii: circles further away
    cannot overlap, so check_overlap gives the same answer as with all circles.
    """
    def __init__(self, max_radius):
        self.spacing = 3*max_radius
        if not self.spacing > 0: # degenerate radii: one grid cell, i.e. check all circles
            self.spacing = math.inf
        self.bins = {}

    def _key(self, circC):
        return (math.floor(circC[0]/self.spacing), math.floor(circC[1]/self.spacing))

    def add(self, circC):
        self.bins.setdefault(self._key(circC), []).append(circC)

    def check_overlap(self, circC):
        """ Same as check_overlap(all placed circles, circC) """
        i, j = self._key(circC)
        neighbours = []
        for key in ((i-1,j-1),(i-1,j),(i-1,j+1),(i,j-1),(i,j),(i,j+1),(i+1,j-1),(i+1,j),(i+1,j+1)):
            neighbours.extend(self.bins.get(key, ()))
        return check_overlap(neighbours, circC)

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


def Clipping_Box(box_shape,box_center,box_size,boundaryFlag,polygon=None):
    """
    clipping box for the 2D Voronoi diagram

    box_shape 'polygon' clips with any simple polygon (convex or not) whose vertices
    [[x, y], ...] are given as `polygon`; box_center and box_size are not used for it.
    """
    
    if box_shape in SQUARE_SHAPES: # sqaure box
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

        if _is_on(boundaryFlag):
            boundarylines = patches.Rectangle(boundaries[0][1][0], x_max-x_min, y_max-y_min, linewidth=2, edgecolor='k',facecolor='none')
        else:
            boundarylines = patches.Rectangle(boundaries[0][1][0], x_max-x_min, y_max-y_min, linewidth=2, linestyle='-.', edgecolor='k',facecolor='none')
            
    elif box_shape in TRIANGLE_SHAPES: # Triangle box
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
        
        if _is_on(boundaryFlag):
            boundarylines = patches.RegularPolygon((box_center[0]+0,box_center[1]+np.sqrt(3)/3*box_size),numVertices=3,radius=2*np.sqrt(3)/3*box_size,orientation=0,linewidth=2,edgecolor='k',facecolor='none')
        else:
            boundarylines = patches.RegularPolygon((box_center[0]+0,box_center[1]+np.sqrt(3)/3*box_size),numVertices=3,radius=2*np.sqrt(3)/3*box_size,orientation=0,linewidth=2,linestyle='-.',edgecolor='k',facecolor='none')
        
        x_min = x_0
        x_max = x_2
        y_min = y_0
        y_max = y_1
        
    elif box_shape in HEXAGON_SHAPES: # hexagonal box
        x_0 = box_center[0] - box_size
        x_1 = box_center[0] - box_size/2.0
        x_2 = box_center[0] + box_size/2.0
        x_3 = box_center[0] + box_size
        y_0 = box_center[1] - box_size/2*np.sqrt(3)
        y_1 = box_center[1]
        y_2 = box_center[1] + box_size/2*np.sqrt(3)
        
        boundary_points = np.array([[x_1, y_0], [x_2, y_0], [x_3, y_1], [x_2, y_2], [x_1, y_2], [x_0, y_1]])
        
        l0 = [(x_1, y_0), (x_2, y_0)]
        l1 = [(x_2, y_0), (x_3, y_1)]
        l2 = [(x_3, y_1), (x_2, y_2)]
        l3 = [(x_2, y_2), (x_1, y_2)]
        l4 = [(x_1, y_2), (x_0, y_1)]
        l5 = [(x_0, y_1), (x_1, y_0)]
        
        boundaries = [('bottom', l0), ('bottom-right', l1), ('top-right', l2), ('top', l3), ('top-left', l4), ('bottom-left', l5)]
        
        if _is_on(boundaryFlag):
            boundarylines = patches.RegularPolygon(box_center,numVertices=6,radius=box_size,orientation=np.pi/6,linewidth=2,edgecolor='k',facecolor='none')
        else:
            boundarylines = patches.RegularPolygon(box_center,numVertices=6,radius=box_size,orientation=np.pi/6,linewidth=2,linestyle='-.',edgecolor='k',facecolor='none')
        
        x_min = x_0
        x_max = x_3
        y_min = y_0
        y_max = y_2
        
    elif box_shape in POLYGON_SHAPES: # any simple polygon
        if polygon is None:
            raise ValueError("box_shape 'polygon' needs the vertices of the polygon (box_polygon)")
        boundary_points = polygon_vertices(polygon)
        corners = [(float(x), float(y)) for x, y in boundary_points]
        boundaries = [('edge{:d}'.format(i), [corners[i], corners[(i+1) % len(corners)]]) for i in range(len(corners))]
        x_min, y_min = boundary_points.min(axis=0)
        x_max, y_max = boundary_points.max(axis=0)
        boundarylines = patches.Polygon(boundary_points,closed=True,linewidth=2,linestyle='-' if _is_on(boundaryFlag) else '-.',
                                        edgecolor='k',facecolor='none')

    else:
        raise ValueError('box_shape: {:s} is not supported for current version, please check README for more details.'.format(box_shape))

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
                         iter_max,print_interval,rng=None):
    """
    packing cells by randomly placing new cells in the generation rings

    rng :: None (global np.random state), an int seed, or a np.random.Generator
    """
    rng = _rng(rng)

    # generate radii for rings
    width = np.concatenate(([width_heart],np.tile([width_sparse,width_dense],nrings)))
    noise = rng.normal(1,0.25,len(width))
    width = np.multiply(width,noise)
    radii = np.concatenate(([0],np.cumsum(width)))
    
    # Place circular cells in each ring (can be parallelized in the future)
    # circles: list = list()
    circles = []
    grid = _CircleGrid(max(cellsize_dense,cellsize_sparse)) # placed circles, for fast overlap checks
    for ibin in range(0,nrings*2+1):
        iter = 0 
        while iter < iter_max:
            r = (radii[ibin+1]-radii[ibin])*math.sqrt(rng.random()) + radii[ibin]  # randomly generate a point in the domain
            t = rng.random()*2*math.pi 
            
            x = r*math.cos(t) + generation_center[0]
            y = r*math.sin(t) + generation_center[1]
    
            if (ibin % 2) == 0: # if even, dense cells
                w = cellsize_dense
            else:
                w = cellsize_sparse
    
            circC = [x, y, w]
            
            if grid.check_overlap(circC) and check_isinside_boundcircle(circC,generation_center,radii[ibin],radii[ibin+1]):
                circles.append(circC)
                grid.add(circC)
                
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
                               print_interval,omega=1.0,rng=None):
    """
    packing cells by firstly placing new cells and then performing Lloyd's relaxation in the generation rings

    rng :: None (global np.random state), an int seed, or a np.random.Generator
    """
    rng = _rng(rng)
    # generate radii for rings
    radii = np.concatenate(([width_heart],np.tile([width_sparse,width_dense],nrings)))
    noise = rng.normal(1,0.25,len(radii))
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
    inside_cells = 1e-3*(rng.random((n_nonoverlapped_cells,2)) - [0.5,0.5]) # add cell points in a very small area
    
    sites = np.vstack((PerimeterPointsSites,inside_cells))

    lloyd_iter = 0
    while lloyd_iter < iter_max and n_nonoverlapped_cells > 0: # a heart narrower than one cell has no cells inside to relax
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
        
        circles = [np.zeros((0,2))] # no sub-rings (nsubrings = 0) if the ring is narrower than one cell
        
        for subr, subnpoint in zip(subradii, subnpoints):
            t = np.linspace(0, 2*np.pi, subnpoint, endpoint=False)
            x = subr * np.cos(t)
            y = subr * np.sin(t)
            circles.append(np.c_[x, y])
            
        inside_cells = np.concatenate(circles, axis=0)
        
        noise = (rng.random((len(inside_cells),2)) - [0.5,0.5])*cellsize
        inside_cells = inside_cells + 0.25*noise

        sites = np.vstack((OuterPerimeterPointsSites,inside_cells))
                    
        existing_sites = np.vstack((sites,existing_sites))
        PerimeterPointsSites = np.copy(OuterPerimeterPointsSites)
    # sites are generated around the origin, move them to the generation center
    sites = existing_sites + np.asarray(generation_center)

    return sites, radii


def CellPlacement_Honeycomb(generation_center,r_max,r_min,nrings,box_center,box_size,\
                            width_heart,width_sparse,width_dense,\
                            cellsize_sparse,cellsize_dense,\
                            cellwallthickness_sparse,cellwallthickness_dense,\
                            iter_max,print_interval,rng=None):
    """
    packing cells in hexagonal grids

    rng :: None (global np.random state), an int seed, or a np.random.Generator
    """
    rng = _rng(rng)
    try:
        from hexalattice.hexalattice import create_hex_grid
    except ImportError as e:
        raise ImportError("radial_growth_rule 'regular_hexagonal' needs the optional package hexalattice: "
                          "pip install 'ringspy[hex]'") from e

    # generate radii for rings
    width = np.concatenate(([width_heart],np.tile([width_sparse,width_dense],nrings)))
    noise = rng.normal(1,0.25,len(width))
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


def _inside_clip_box(box_shape,boundaries,x_min,x_max,y_min,y_max,box_center):
    """Return a function telling whether a 2D point lies strictly inside the clipping box"""
    if box_shape in TRIANGLE_SHAPES:
        return lambda p: check_isinside_boundTriangle2D(p,boundaries)
    elif box_shape in SQUARE_SHAPES:
        return lambda p: check_isinside_boundbox2D(np.append(p,0),x_min,x_max,y_min,y_max)
    elif box_shape in HEXAGON_SHAPES:
        L = y_max - y_min
        return lambda p: check_isinside_boundHexagon2D(p,L,box_center)
    elif box_shape in POLYGON_SHAPES:
        return _PolygonInside([boundary[1][0] for boundary in boundaries])
    elif box_shape == 'notched_square': # notched
        x_min = boundaries[0][1][0][0]
        x_max = boundaries[1][1][0][0]
        y_min = boundaries[0][1][0][1]
        y_max = boundaries[2][1][0][1]
        x_indent = boundaries[5][1][0][0]
        y_indent_min = boundaries[6][1][0][1]
        y_indent_max = boundaries[4][1][0][1]
        return lambda p: bool(check_isinside_boundbox2Dindent(np.append(p,0),x_min,x_max,y_min,y_max,x_indent,y_indent_min,y_indent_max))
    else:
        raise ValueError('box_shape: {:s} is not supported for current version, please check README for more details.'.format(box_shape))


def _ray_boundary_intersection(p,direction,boundaries):
    """
    Closest intersection point of the ray p + t*direction (t >= 0) with the boundary line segments
    Ref: https://stackoverflow.com/questions/14307158/how-do-you-check-for-intersection-between-a-line-segment-and-a-line-ray-emanatin
    """
    t_final = np.inf
    for boundary in boundaries: # loop over boundary lines
        q = np.asarray(boundary[1][0],dtype=float)
        s = np.asarray(boundary[1][1],dtype=float) - q
        denom = cross2D(direction,s)
        if denom == 0: # ray parallel to the boundary line
            continue
        t = cross2D((q-p),s)/denom
        u = cross2D((q-p),direction)/denom
        if (-1e-9 <= u <= 1+1e-9) and (0 <= t < t_final):
            t_final = t
    if not np.isfinite(t_final):
        raise RuntimeError('Voronoi ridge starting at {} does not cross the clipping box boundary.'.format(p))
    return p + direction*t_final


def _boundary_ridges_along_perimeter(boundary_points,boundaries):
    """
    Connect the boundary points into a closed line path by sorting them by their
    position along the perimeter of the clipping box (boundary lines are ordered
    around the box). Returns pairs of indices of boundary_points.
    """
    starts = np.array([boundary[1][0] for boundary in boundaries],dtype=float)
    ends = np.array([boundary[1][1] for boundary in boundaries],dtype=float)
    segments = ends - starts
    # closest position of each point on each boundary line
    rel = boundary_points[:,np.newaxis,:] - starts[np.newaxis,:,:]
    u = np.clip(np.einsum('nmk,mk->nm',rel,segments)/np.einsum('mk,mk->m',segments,segments),0,1)
    closest = starts[np.newaxis,:,:] + u[:,:,np.newaxis]*segments[np.newaxis,:,:]
    iboundary = np.linalg.norm(boundary_points[:,np.newaxis,:] - closest,axis=2).argmin(axis=1)
    perimeter_position = iboundary + u[np.arange(len(iboundary)),iboundary]
    order = np.argsort(perimeter_position,kind='stable')
    return np.column_stack((order,np.roll(order,-1)))


def _project_into_convex_box(points,boundaries):
    """
    Move points lying outside a convex clipping box (boundary lines ordered
    counter-clockwise) onto the closest point of its boundary
    """
    starts = np.array([boundary[1][0] for boundary in boundaries],dtype=float)
    ends = np.array([boundary[1][1] for boundary in boundaries],dtype=float)
    segments = ends - starts
    lengths = np.linalg.norm(segments,axis=1)
    # signed distances to the boundary lines (positive inside)
    distance = cross2D(segments[np.newaxis,:,:],points[:,np.newaxis,:] - starts[np.newaxis,:,:])/lengths
    outside = distance.min(axis=1) < -1e-12*lengths.max()
    if not np.any(outside):
        return points
    rel = points[outside][:,np.newaxis,:] - starts[np.newaxis,:,:]
    u = np.clip(np.einsum('nmk,mk->nm',rel,segments)/lengths**2,0,1)
    closest = starts[np.newaxis,:,:] + u[:,:,np.newaxis]*segments[np.newaxis,:,:]
    nearest = np.linalg.norm(rel - u[:,:,np.newaxis]*segments[np.newaxis,:,:],axis=2).argmin(axis=1)
    points = np.copy(points)
    points[outside] = closest[np.arange(len(nearest)),nearest]
    return points


def _segments_intersect(p1,p2,q1,q2):
    """Whether the closed line segments p1-p2 and q1-q2 have a point in common (arrays of segments
    with shapes broadcasting together, the last axis holding x and y)"""
    def orientation(a,b,c):
        return np.sign(cross2D(b - a,c - a))
    def within(a,b,c): # c on the line through a and b: is it between them?
        return (np.minimum(a[...,0],b[...,0]) <= c[...,0]) & (c[...,0] <= np.maximum(a[...,0],b[...,0])) & \
               (np.minimum(a[...,1],b[...,1]) <= c[...,1]) & (c[...,1] <= np.maximum(a[...,1],b[...,1]))
    o1, o2 = orientation(p1,p2,q1), orientation(p1,p2,q2)
    o3, o4 = orientation(q1,q2,p1), orientation(q1,q2,p2)
    return ((o1*o2 < 0) & (o3*o4 < 0)) | \
           ((o1 == 0) & within(p1,p2,q1)) | ((o2 == 0) & within(p1,p2,q2)) | \
           ((o3 == 0) & within(q1,q2,p1)) | ((o4 == 0) & within(q1,q2,p2))


def polygon_vertices(polygon):
    """
    Check the vertices [[x, y], ...] of a clipping polygon and return them as an (n, 2) array in
    counter-clockwise order. The polygon must be simple: at least 3 distinct vertices, a non-zero
    area and no crossing or touching edges. A vertex at the place of the previous one (within 1e-9 of the
    polygon's size) is dropped, also a last vertex at the place of the first one (a closed list, or a
    circle from np.linspace(0, 2*pi, n)).
    """
    p = np.asarray(polygon,dtype=float)
    if p.ndim != 2 or p.shape[1] != 2 or not np.all(np.isfinite(p)):
        raise ValueError('polygon vertices must be a list of [x, y] pairs of numbers, not {!r}'.format(polygon))
    tol = 1e-9*np.ptp(p,axis=0).max() if len(p) else 0.0
    rows = [] # the vertices kept, by their row in the list given (the last of repeated ones: its edge follows)
    for k in range(len(p)):
        if not rows or np.linalg.norm(p[k] - p[rows[-1]]) > tol:
            rows.append(k)
        else:
            rows[-1] = k
    if len(rows) > 1 and np.linalg.norm(p[rows[-1]] - p[rows[0]]) <= tol:
        rows.pop()
    p = p[rows]
    if len(p) < 3:
        raise ValueError('a polygon needs at least 3 distinct vertices')
    n = len(p)
    a, b = p, np.roll(p,-1,axis=0)
    j = np.arange(n)
    block = max(1,2000000//n)
    for start in range(0,n,block): # edges i of the block against all edges j > i+1
        i = np.arange(start,min(start+block,n))[:,np.newaxis]
        pairs = (j > i+1) & ~((i == 0) & (j == n-1)) # adjacent edges share a vertex
        cross = pairs & _segments_intersect(a[i],b[i],a[j],b[j])
        if np.any(cross):
            ii, jj = np.unravel_index(np.argmax(cross),cross.shape)
            # edge k of the list given runs from vertex k to the next one
            raise ValueError('the polygon is not simple: edges {} and {} cross or touch'.format(rows[start+ii],rows[jj]))
    area = 0.5*np.sum(cross2D(p,np.roll(p,-1,axis=0)))
    if abs(area) <= 1e-12*np.ptp(p,axis=0).max()**2:
        raise ValueError('the polygon has no area')
    if area < 0:
        p = p[::-1]
    return p


class _PolygonInside:
    """Inside test for a simple polygon (convex or not); points on the boundary are outside"""

    def __init__(self, polygon):
        self.polygon = np.asarray(polygon,dtype=float)
        self.a = self.polygon
        self.b = np.roll(self.polygon,-1,axis=0)
        self.scale = np.ptp(self.polygon,axis=0).max()
        self.tol = 1e-9*self.scale

    def classify(self, points):
        """(inside by the crossing number, on the boundary within tol) for an array of points"""
        points = np.asarray(points,dtype=float).reshape(-1,2)
        chunk = max(256,2000000//len(self.a)) # points per pass: bounded memory for polygons with many edges
        inside = np.zeros(len(points),dtype=bool)
        on = np.zeros(len(points),dtype=bool)
        xa, ya = self.a[:,0][np.newaxis], self.a[:,1][np.newaxis]
        xb, yb = self.b[:,0][np.newaxis], self.b[:,1][np.newaxis]
        ex, ey = xb - xa, yb - ya
        length2 = ex**2 + ey**2
        for k in range(0,len(points),chunk):
            x, y = points[k:k+chunk,0:1], points[k:k+chunk,1:2]
            straddle = (ya > y) != (yb > y)
            with np.errstate(divide='ignore',invalid='ignore'):
                x_cross = xa + (y - ya)*ex/ey
            inside[k:k+chunk] = np.count_nonzero(straddle & (x < x_cross),axis=1) % 2 == 1
            u = np.clip(((x - xa)*ex + (y - ya)*ey)/length2,0,1)
            on[k:k+chunk] = np.min((x - xa - u*ex)**2 + (y - ya - u*ey)**2,axis=1) <= self.tol**2
        return inside, on

    def many(self, points):
        inside, on = self.classify(points)
        return inside & ~on

    def __call__(self, p):
        return bool(self.many(np.asarray(p,dtype=float)[0:2])[0])

    def distance(self, points):
        """Distance of each point to the polygon's boundary"""
        points = np.asarray(points,dtype=float).reshape(-1,2)
        e = self.b - self.a
        length2 = np.einsum('mk,mk->m',e,e)
        distance = np.empty(len(points))
        chunk = max(256,2000000//len(self.a))
        for k in range(0,len(points),chunk):
            rel = points[k:k+chunk,np.newaxis,:] - self.a[np.newaxis,:,:]
            u = np.clip(np.einsum('nmk,mk->nm',rel,e)/length2,0,1)
            distance[k:k+chunk] = np.linalg.norm(rel - u[:,:,np.newaxis]*e[np.newaxis,:,:],axis=2).min(axis=1)
        return distance

    def segments_within(self, p, q):
        """Whether each segment p[k]-q[k] lies in the polygon, its boundary included (arrays (n, 2)).
        The segments are cut where they cross the polygon's edges; each piece must be inside or on."""
        p, q = np.asarray(p,dtype=float).reshape(-1,2), np.asarray(q,dtype=float).reshape(-1,2)
        d = q - p
        edges = self.b - self.a
        denom = cross2D(d[:,np.newaxis,:],edges[np.newaxis,:,:])
        with np.errstate(divide='ignore',invalid='ignore'):
            t = cross2D(self.a[np.newaxis,:,:] - p[:,np.newaxis,:],edges[np.newaxis,:,:])/denom
            u = cross2D(self.a[np.newaxis,:,:] - p[:,np.newaxis,:],d[:,np.newaxis,:])/denom
        hit = (denom != 0) & (u >= -1e-9) & (u <= 1 + 1e-9) & (t > 1e-9) & (t < 1 - 1e-9)
        k = int(hit.sum(axis=1).max()) if len(p) else 0
        breaks = np.hstack((np.zeros((len(p),1)),np.sort(np.where(hit,t,1.0),axis=1)[:,:k],np.ones((len(p),1))))
        middles = p[:,np.newaxis,:] + 0.5*(breaks[:,:-1] + breaks[:,1:])[:,:,np.newaxis]*d[:,np.newaxis,:]
        inside, on = self.classify(middles)
        return (inside | on).reshape(len(p),k+1).all(axis=1)

    def segments_inside(self, p, q):
        """Whether each segment p[k]-q[k] lies strictly inside the polygon between its ends: it does not touch
        the boundary there (its ends may lie on it)"""
        p, q = np.asarray(p,dtype=float).reshape(-1,2), np.asarray(q,dtype=float).reshape(-1,2)
        d = q - p
        edges = self.b - self.a
        denom = cross2D(d[:,np.newaxis,:],edges[np.newaxis,:,:])
        with np.errstate(divide='ignore',invalid='ignore'):
            t = cross2D(self.a[np.newaxis,:,:] - p[:,np.newaxis,:],edges[np.newaxis,:,:])/denom
            u = cross2D(self.a[np.newaxis,:,:] - p[:,np.newaxis,:],d[:,np.newaxis,:])/denom
        touch = ((denom != 0) & (u >= -1e-9) & (u <= 1 + 1e-9) & (t > 1e-9) & (t < 1 - 1e-9)).any(axis=1)
        inside, on = self.classify(0.5*(p + q))
        return ~touch & inside & ~on


def _polygon_crossings(p,d,tmax,starts,edges):
    """Sorted parameters 0 < t < tmax where p + t*d crosses the edges starts + u*edges (0 <= u <= 1)"""
    denom = cross2D(d,edges)
    with np.errstate(divide='ignore',invalid='ignore'):
        t = cross2D(starts - p,edges)/denom
        u = cross2D(starts - p,d)/denom
    hit = (np.abs(denom) > 1e-14*np.linalg.norm(d)*np.linalg.norm(edges,axis=1)) & \
          (u >= -1e-9) & (u <= 1 + 1e-9) & (t > 0) & (t < tmax)
    return np.sort(t[hit])


def _clip_voronoi_polygon(vor,circles,boundaries,generation_center,inside,boundaryFlag,plot=False):
    """
    _clip_voronoi for box_shape 'polygon' (any simple polygon, also non-convex). Each ridge (and the
    ray of each infinite ridge) is cut into its parts inside the polygon. A ridge can leave and enter
    a non-convex polygon more than once, so besides ridges inside the polygon and cut ridges from an
    inside vertex to the boundary, there can be ridges between two boundary points; they are
    returned with the cut ridges in infinite_ridges_new (only with boundaryFlag on: without the walls
    along the boundary they would be loose pieces). Same return values as _clip_voronoi, except that
    walls ending at the same boundary point share it, and walls ending at a corner end at the corner.
    """
    generation_center = np.asarray(generation_center,dtype=float)
    vertices = vor.vertices
    # a vertex this close to the boundary is a boundary point: a wall from it running along an edge lies
    # within tol of the boundary and is dropped, so the vertex must not stay an inside vertex without it
    vertex_in, vertex_on = inside.classify(vertices)
    vertex_on |= vertex_in & (inside.distance(vertices) <= 4*inside.tol)
    vertex_in &= ~vertex_on
    starts, edges = inside.a, inside.b - inside.a
    eps = 1e-12

    # finite ridges that touch or cross the boundary (vectorized); the others are inside or outside
    ridge_vertices = [np.asarray(simplex) for simplex in vor.ridge_vertices]
    finite = np.array([np.all(simplex >= 0) for simplex in ridge_vertices],dtype=bool)
    pairs = np.array([simplex if f else [0, 0] for simplex, f in zip(ridge_vertices,finite)],dtype=int).reshape(-1,2)
    near = np.zeros(len(pairs),dtype=bool)
    chunk = max(256,2000000//len(edges)) # ridges per pass: bounded memory for polygons with many edges
    for k in range(0,len(pairs),chunk):
        p = vertices[pairs[k:k+chunk,0]][:,np.newaxis,:]
        d = vertices[pairs[k:k+chunk,1]][:,np.newaxis,:] - p
        denom = cross2D(d,edges[np.newaxis,:,:])
        with np.errstate(divide='ignore',invalid='ignore'):
            t = cross2D(starts[np.newaxis,:,:] - p,edges[np.newaxis,:,:])/denom
            u = cross2D(starts[np.newaxis,:,:] - p,d)/denom
        near[k:k+chunk] = np.any((denom != 0) & (u >= -1e-9) & (u <= 1 + 1e-9) & (t >= -1e-9) & (t <= 1 + 1e-9),axis=1)
    near |= vertex_on[pairs[:,0]] | vertex_on[pairs[:,1]] # ends on the boundary within tol

    finite_ridges = []
    cut_vertex = [] # Voronoi vertex inside the polygon of each cut ridge
    cut_first = [] # whether that vertex is the first one of the cut ridge
    cut_boundary = [] # index of the boundary point of each cut ridge
    boundary_boundary = [] # ridges between two boundary points
    boundary_points = []

    def boundary_point(x):
        boundary_points.append(x)
        return len(boundary_points) - 1

    def at_corner(x):
        return np.min(np.sum((inside.polygon - x)**2,axis=1)) <= inside.tol**2

    def inside_parts(p,d,tmax,start,end):
        """Parts of p + t*d (0 <= t <= tmax) inside the polygon, each as a pair of ends; an end is the
        index of the Voronoi vertex `start`/`end` when the part ends there, else a point on the boundary"""
        breaks = np.concatenate(([0.0],_polygon_crossings(p,d,tmax,starts,edges),[tmax]))
        parts = []
        for t0, t1 in zip(breaks[:-1],breaks[1:]):
            if t1 - t0 <= eps*tmax:
                continue
            if inside(p + d*(0.5*(t0+t1))):
                # the same crossing found twice: same part; but with boundaryFlag, a wall touching a corner of
                # the polygon from inside is cut there into two parts ending at the corner (a vertex of the mesh)
                if parts and abs(parts[-1][1] - t0) <= eps*tmax and not (_is_on(boundaryFlag) and at_corner(p + d*t0)):
                    parts[-1][1] = t1
                else:
                    parts.append([t0, t1])
        ends = []
        for t0, t1 in parts:
            a = start if (t0 == 0 and start is not None and vertex_in[start]) else p + d*t0
            b = end if (t1 == tmax and end is not None and vertex_in[end]) else p + d*t1
            ends.append((a, b))
        return ends

    def add(a,b):
        a_vertex, b_vertex = not isinstance(a,np.ndarray), not isinstance(b,np.ndarray)
        if a_vertex and b_vertex:
            finite_ridges.append(np.array([a, b]))
        elif a_vertex or b_vertex:
            cut_vertex.append(a if a_vertex else b)
            cut_first.append(a_vertex)
            cut_boundary.append(boundary_point(b if a_vertex else a))
        elif _is_on(boundaryFlag): # without the walls along the boundary it would be a loose piece
            boundary_boundary.append([boundary_point(a), boundary_point(b)])

    for k, (pointidx, simplex) in enumerate(zip(vor.ridge_points,ridge_vertices)):
        if finite[k]:
            if not near[k]:
                if vertex_in[simplex[0]] and vertex_in[simplex[1]]:
                    finite_ridges.append(simplex)
                continue
            p = vertices[simplex[0]]
            for a, b in inside_parts(p,vertices[simplex[1]] - p,1.0,simplex[0],simplex[1]):
                add(a,b)
        else: # infinite Voronoi ridge: a ray from its finite vertex, facing outwards
            i = simplex[simplex >= 0][0]
            tangent = circles[pointidx[1],0:2] - circles[pointidx[0],0:2]
            tangent = tangent/np.linalg.norm(tangent)
            normal = np.array([-tangent[1], tangent[0]])
            midpoint = circles[pointidx,0:2].mean(axis=0) # midpoint between site points
            if np.dot(midpoint - generation_center, normal) < 0:
                normal = -normal
            p = vertices[i]
            crossings = _polygon_crossings(p,normal,np.inf,starts,edges)
            if len(crossings) == 0:
                continue
            for a, b in inside_parts(p,normal,crossings[-1] + inside.scale,i,None):
                add(a,b) if simplex[0] >= 0 else add(b,a)

    if len(finite_ridges) == 0:
        raise ValueError('No Voronoi ridge has both ends strictly inside the clipping polygon: the polygon holds no '
                         'complete cell wall (it is only a few cells large or one cell wide, or outside the cells); '
                         'please check the generation and clipping box parameters.')
    finite_ridges = np.asarray(finite_ridges,dtype=int).reshape(-1,2)
    nfinite_ridge = finite_ridges.shape[0]

    # unique vertices inside the polygon (vertices with identical coordinates are merged)
    used_vertices = np.concatenate((finite_ridges.ravel(),np.asarray(cut_vertex,dtype=int)))
    voronoi_vertices_in, inverse = np.unique(vertices[used_vertices],axis=0,return_inverse=True)
    inverse = inverse.reshape(-1)
    nvertices_in = voronoi_vertices_in.shape[0]
    finite_ridges_new = inverse[:2*nfinite_ridge].reshape(-1,2)

    # boundary points at the same place (walls meeting on the boundary at a Voronoi vertex, or cut at a
    # corner) become one point; with boundaryFlag on, a point at a corner of the polygon is that corner
    boundary_points = np.asarray(boundary_points,dtype=float).reshape(-1,2)
    corners = inside.polygon if _is_on(boundaryFlag) else np.zeros((0,2))
    candidates = np.vstack((boundary_points,corners))
    close = cKDTree(candidates).query_pairs(r=inside.tol,output_type='ndarray')
    group = np.arange(len(candidates))
    if len(close):
        from scipy.sparse import coo_matrix
        from scipy.sparse.csgraph import connected_components
        graph = coo_matrix((np.ones(len(close)),(close[:,0],close[:,1])),shape=(len(candidates),)*2)
        group = connected_components(graph,directed=False)[1]
    corner_of_group = np.full(len(candidates),-1)
    corner_of_group[group[len(boundary_points):]] = np.arange(len(corners))
    first = np.unique(group[:len(boundary_points)],return_index=True)[1] # first point of each group, in order
    kept = np.sort(first[corner_of_group[group[first]] < 0])
    ncut = len(kept)
    new_index = np.full(len(candidates),-1)
    new_index[group[kept]] = np.arange(ncut)
    corner_groups = corner_of_group >= 0
    new_index[np.nonzero(corner_groups)[0]] = ncut + corner_of_group[corner_groups]
    point_index = new_index[group[:len(boundary_points)]] # index of each boundary point among [kept, corners]
    boundary_points = boundary_points[kept]

    cut_inside = inverse[2*nfinite_ridge:]
    cut_point = nvertices_in + point_index[np.asarray(cut_boundary,dtype=int)]
    cut_ridges = np.where(np.asarray(cut_first,dtype=bool)[:,np.newaxis],
                          np.column_stack((cut_inside,cut_point)),
                          np.column_stack((cut_point,cut_inside))).reshape(-1,2)
    boundary_boundary = nvertices_in + point_index[np.asarray(boundary_boundary,dtype=int).reshape(-1,2)]
    boundary_boundary = boundary_boundary[boundary_boundary[:,0] != boundary_boundary[:,1]] # both ends at one point
    infinite_ridges_new = np.vstack((cut_ridges,boundary_boundary)).astype(int)

    if _is_on(boundaryFlag):
        boundary_points_featured = np.reshape(np.asarray([x[1][0] for x in boundaries],dtype=float),(-1,2))
        boundary_points = np.vstack((boundary_points,boundary_points_featured))
    else:
        boundary_points_featured = np.zeros((0,2))

    if plot:
        _plot_segments(vertices[finite_ridges],'k',linewidth=1)
        if len(infinite_ridges_new) > 0:
            all_points = np.vstack((voronoi_vertices_in,boundary_points))
            _plot_segments(all_points[infinite_ridges_new],'k',linewidth=1)

    return voronoi_vertices_in,boundary_points,boundary_points_featured,ncut,\
        finite_ridges,finite_ridges_new,infinite_ridges_new


def _plot_segments(segments,color,linestyle='-',linewidth=None):
    """
    Draw line segments (array of shape (n,2,2)) in the current axes as one LineCollection,
    styled like plt.plot(segments[:,:,0].T, segments[:,:,1].T, ...), which creates one
    Line2D per segment and is slow to draw for large meshes
    """
    from matplotlib.collections import LineCollection
    if linewidth is None:
        linewidth = plt.rcParams['lines.linewidth']
    solid = linestyle in ['-','solid']
    lines = LineCollection(np.asarray(segments,dtype=float).reshape(-1,2,2),colors=color,linestyles=linestyle,
                           linewidths=linewidth,
                           capstyle=plt.rcParams['lines.solid_capstyle' if solid else 'lines.dash_capstyle'],
                           joinstyle=plt.rcParams['lines.solid_joinstyle' if solid else 'lines.dash_joinstyle'])
    plt.gca().add_collection(lines)
    return lines


def _clip_voronoi(vor,circles,boundaries,generation_center,inside,boundaryFlag,plot=False):
    """
    Clip the Voronoi diagram by the clipping box. Ridges with both vertices inside
    are kept; finite ridges with one vertex outside and infinite ridges are cut at
    the box boundary.

    Returns
    -------
    voronoi_vertices_in :: unique vertices inside the box
    boundary_points :: cut points (one per cut ridge), followed by the box corners if boundaryFlag is on
    boundary_points_featured :: the box corners appended to boundary_points
    ncut :: number of cut points
    finite_ridges :: ridges inside the box (indices of vor.vertices)
    finite_ridges_new :: ridges inside the box (indices of voronoi_vertices_in)
    infinite_ridges_new :: cut ridges (indices of [voronoi_vertices_in, boundary_points])
    """
    if isinstance(inside,_PolygonInside):
        return _clip_voronoi_polygon(vor,circles,boundaries,generation_center,inside,boundaryFlag,plot)
    generation_center = np.asarray(generation_center)
    vertex_in = np.array([inside(vertex) for vertex in vor.vertices],dtype=bool)

    finite_ridges = []
    cut_vertex = [] # index of the inside vertex of each cut ridge
    cut_first = [] # whether the inside vertex is the first vertex of the cut ridge
    cut_direction = [] # direction from the inside vertex towards the box boundary
    for pointidx, simplex in zip(vor.ridge_points, vor.ridge_vertices):
        simplex = np.asarray(simplex)
        if np.all(simplex >= 0):
            in0 = vertex_in[simplex[0]]
            in1 = vertex_in[simplex[1]]
            if in0 and in1:
                finite_ridges.append(simplex)
            elif in0 or in1: # finite Voronoi ridge but one vertex is out of bound
                i, j = (simplex[0], simplex[1]) if in0 else (simplex[1], simplex[0])
                direction = vor.vertices[j] - vor.vertices[i]
                cut_vertex.append(i)
                cut_first.append(in0)
                cut_direction.append(direction/np.linalg.norm(direction))
        else: # infinite Voronoi ridge
            i = simplex[simplex >= 0][0] # index of finite end Voronoi vertex of the infinite ridge
            if vertex_in[i]:
                tangent = circles[pointidx[1],0:2] - circles[pointidx[0],0:2]
                tangent = tangent/np.linalg.norm(tangent)
                normal = np.array([-tangent[1], tangent[0]])
                midpoint = circles[pointidx,0:2].mean(axis=0) # midpoint between site points
                if np.dot(midpoint - generation_center, normal) < 0: # make the ray face outwards
                    normal = -normal
                cut_vertex.append(i)
                cut_first.append(simplex[0] >= 0)
                cut_direction.append(normal)

    if len(finite_ridges) == 0:
        raise ValueError('No Voronoi ridge lies inside the clipping box, please check the generation and clipping box parameters.')
    finite_ridges = np.asarray(finite_ridges)
    nfinite_ridge = finite_ridges.shape[0]
    ncut = len(cut_vertex)

    # unique vertices inside the box (vertices with identical coordinates are merged)
    used_vertices = np.concatenate((finite_ridges.ravel(),np.asarray(cut_vertex,dtype=int)))
    voronoi_vertices_in, inverse = np.unique(vor.vertices[used_vertices],axis=0,return_inverse=True)
    inverse = inverse.reshape(-1)
    nvertices_in = voronoi_vertices_in.shape[0]
    finite_ridges_new = inverse[:2*nfinite_ridge].reshape(-1,2)

    # Find the intersect points of the cut ridges with the boundary lines
    boundary_points = np.array([_ray_boundary_intersection(vor.vertices[i],direction,boundaries) \
                                for i, direction in zip(cut_vertex,cut_direction)]).reshape(-1,2)
    cut_inside = inverse[2*nfinite_ridge:]
    cut_boundary = nvertices_in + np.arange(ncut)
    infinite_ridges_new = np.where(np.asarray(cut_first,dtype=bool)[:,np.newaxis],
                                   np.column_stack((cut_inside,cut_boundary)),
                                   np.column_stack((cut_boundary,cut_inside))).reshape(-1,2).astype(int)

    if plot:
        _plot_segments(vor.vertices[finite_ridges],'k',linewidth=1)
        if ncut > 0:
            _plot_segments(np.stack((vor.vertices[cut_vertex],boundary_points),axis=1),'k',linewidth=1)

    if _is_on(boundaryFlag):
        boundary_points_featured = np.reshape(np.asarray([x[1][0] for x in boundaries],dtype=float),(-1,2))
        boundary_points = np.vstack((boundary_points,boundary_points_featured))
    else:
        boundary_points_featured = np.zeros((0,2))

    return voronoi_vertices_in,boundary_points,boundary_points_featured,ncut,\
        finite_ridges,finite_ridges_new,infinite_ridges_new


def RebuildVoronoi(vor,circles,boundaries,generation_center,x_min,x_max,y_min,y_max,box_center,box_shape,boundaryFlag):
    """Clip Voronoi mesh by the boundaries, rebuild the new Voronoi mesh"""
    inside = _inside_clip_box(box_shape,boundaries,x_min,x_max,y_min,y_max,box_center)
    [voronoi_vertices_in,boundary_points,boundary_points_featured,nboundary_pts_ini,\
     finite_ridges,finite_ridges_new,infinite_ridges_new] = \
        _clip_voronoi(vor,circles,boundaries,generation_center,inside,boundaryFlag,plot=True)

    nvertices_in = voronoi_vertices_in.shape[0]
    nboundary_pts = boundary_points.shape[0]
    nboundary_pts_featured = boundary_points_featured.shape[0]
    nfinite_ridge = finite_ridges.shape[0]

    if _is_on(boundaryFlag):
        # construct the connectivity for a line path consisting of the boundary points
        boundary_ridges_new = _boundary_ridges_along_perimeter(boundary_points,boundaries) + nvertices_in
        boundary_ridges_new = np.vstack((infinite_ridges_new,boundary_ridges_new))
    else:
        boundary_ridges_new = np.copy(infinite_ridges_new)

    voronoi_vertices = np.vstack((voronoi_vertices_in,boundary_points)) # vertical stack "in" vetices and "cross" boundary points
    nvertex = voronoi_vertices.shape[0]
    nboundary_ridge = boundary_ridges_new.shape[0]
    voronoi_ridges = np.vstack((finite_ridges_new,boundary_ridges_new))
    nridge = voronoi_ridges.shape[0]

    return voronoi_vertices,finite_ridges,boundary_points,finite_ridges_new,\
        boundary_ridges_new,nvertex,nvertices_in,nfinite_ridge,nboundary_ridge,\
            nboundary_pts,nboundary_pts_featured,voronoi_ridges,nridge

class _PointIndex:
    """
    Row index of 2D points by their coordinates, for a list of points that only grows.
    lookup(p) gives the first row equal to p, same as np.where(np.all(points == p, axis=1))[0][0]
    """
    def __init__(self, points):
        self.index = {}
        self.npoints = 0
        self.extend(points)

    def extend(self, points):
        for point in map(tuple, np.reshape(points,(-1,2)).tolist()):
            self.index.setdefault(point, self.npoints)
            self.npoints += 1

    def lookup(self, points):
        """Rows of an array of points with shape (..., 2), as an array with shape (...)"""
        points = np.asarray(points)
        rows = [self.index[point] for point in map(tuple, np.reshape(points,(-1,2)).tolist())]
        return np.array(rows,dtype=np.intp).reshape(points.shape[:-1])


def _replace_merged_indices(ridge_lists,pairs,merged_indices,lo,hi):
    """
    Replace the vertex indices of each pair of merged vertices with the index of
    their merged vertex: same result as, for each pair in order,
        ridges[ridges == pair[0]] = merged_index; ridges[ridges == pair[1]] = merged_index
    for all ridge lists. The vertex indices of the pairs are in [lo, hi) and the merged
    indices are not, so a vertex in several pairs gets the merged index of its first pair.
    """
    lut = np.arange(lo,hi)
    replaced = np.zeros(hi-lo,dtype=bool)
    for pair, merged_index in zip(np.reshape(pairs,(-1,2)).tolist(), np.reshape(merged_indices,-1).tolist()):
        for vertex in pair:
            if not replaced[vertex-lo]:
                lut[vertex-lo] = merged_index
                replaced[vertex-lo] = True
    replaced_lists = []
    for ridges in ridge_lists:
        ridges = np.copy(ridges)
        inside = (ridges >= lo) & (ridges < hi)
        ridges[inside] = lut[ridges[inside]-lo]
        replaced_lists.append(ridges)
    return replaced_lists


def _rebuild_voronoi_merge_polygon(vor,circles,boundaries,generation_center,inside,merge_tol,boundaryFlag):
    """
    RebuildVoronoi_merge for box_shape 'polygon' (any simple polygon, also non-convex). As for the other
    shapes, vertices of the clipped mesh closer than merge_tol are fused: two vertices inside the polygon at
    their midpoint (Case 1), two neighbouring boundary points along the perimeter at their midpoint (Case 2)
    or onto a corner, and a vertex inside onto a boundary point or a corner (Case 3); without boundaryFlag
    only Case 1. Here the closest pairs go first, until no pair is left, and the fused vertices are not
    rounded. A fusion is skipped when the fused vertex or one of its walls would leave the polygon (across
    a notch, a neck or a reflex corner), so the walls stay in the polygon and the corners stay in place.
    """
    flag = _is_on(boundaryFlag)
    [voronoi_vertices_in,boundary_points,boundary_points_featured,ncut,\
     finite_ridges,finite_ridges_new,infinite_ridges_new] = \
        _clip_voronoi_polygon(vor,circles,boundaries,generation_center,inside,boundaryFlag)
    nin = voronoi_vertices_in.shape[0]
    position = np.vstack((voronoi_vertices_in,boundary_points))
    kind = np.zeros(len(position),dtype=int) # 0 inside, 1 on an edge, 2 corner (boundaryFlag on)
    kind[nin:] = 1
    kind[nin+ncut:] = 2

    if flag:
        chain = _boundary_ridges_along_perimeter(boundary_points,boundaries) + nin
    else:
        chain = np.zeros((0,2),dtype=int)
    walls = np.vstack((finite_ridges_new,infinite_ridges_new,chain)).astype(int).tolist()
    group = [0]*len(finite_ridges_new) + [1]*len(infinite_ridges_new) + [2]*len(chain) # finite, cut, boundary
    alive_wall = [True]*len(walls)
    incident = [set() for _ in range(len(position))]
    wall_of = {} # (smaller, larger vertex) -> wall
    for w, (a, b) in enumerate(walls):
        key = (min(a,b),max(a,b))
        if key in wall_of: # the same wall twice: keep one
            group[wall_of[key]] = max(group[wall_of[key]],group[w])
            alive_wall[w] = False
            continue
        wall_of[key] = w
        incident[a].add(w)
        incident[b].add(w)
    alive = np.ones(len(position),dtype=bool)
    # a lower bound of each vertex's distance to the boundary: a disk of that radius around it is in the polygon
    clearance = np.where(kind == 0,inside.distance(position),0.0)

    def neighbours(v):
        return {walls[w][0] if walls[w][1] == v else walls[w][1] for w in incident[v]}

    def plan(u,v):
        """(kept vertex, fused vertex, position and clearance of the kept vertex) or None if u and v are not
        fused"""
        if kind[u] < kind[v]:
            u, v = v, u
        if kind[v] == 2 or (not flag and kind[u] > 0):
            return None
        margin = 0.0
        if kind[v] == 1: # two points on the perimeter, neighbours along it (Case 2, or onto a corner)
            w = wall_of.get((min(u,v),max(u,v)))
            if w is None or group[w] != 2:
                return None
            target = position[u] if kind[u] == 2 else 0.5*(position[u] + position[v])
        elif kind[u] > 0: # a vertex inside onto a boundary point or a corner (Case 3)
            target = position[u]
        else:
            target = 0.5*(position[u] + position[v])
            margin = max(clearance[u],clearance[v]) - 0.5*np.linalg.norm(position[u] - position[v])
        others = list((neighbours(u) | neighbours(v)) - {u,v})
        if margin > 0 and (not others or np.linalg.norm(position[others] - target,axis=1).max() < margin):
            return u, v, target, margin # the fused vertex and its walls lie in a disk inside the polygon
        if kind[u] == 0:
            if not inside(target):
                return None
            margin = inside.distance(target)[0]
        if others and not inside.segments_within(np.tile(target,(len(others),1)),position[others]).all():
            return None
        # a moved wall that is not a wall along the boundary must not touch the boundary between its ends (it
        # would lie along an edge over other boundary vertices, or pass through a corner)
        rest = [o for o in others if not any(group[wall_of[(min(x,o),max(x,o))]] == 2
                                             for x in (u,v) if (min(x,o),max(x,o)) in wall_of)]
        if rest and not inside.segments_inside(np.tile(target,(len(rest),1)),position[rest]).all():
            return None
        return u, v, target, margin

    def fuse(u,v,target,margin):
        position[u] = target
        clearance[u] = margin
        alive[v] = False
        for w in list(incident[v]):
            other = walls[w][0] if walls[w][1] == v else walls[w][1]
            incident[v].discard(w)
            del wall_of[(min(v,other),max(v,other))]
            if other == u: # the wall between them
                alive_wall[w] = False
                incident[u].discard(w)
                continue
            existing = wall_of.get((min(u,other),max(u,other)))
            if existing is not None: # u and v were both joined to other: one wall
                group[existing] = max(group[existing],group[w])
                alive_wall[w] = False
                incident[other].discard(w)
            else:
                walls[w] = [u if x == v else x for x in walls[w]]
                wall_of[(min(u,other),max(u,other))] = w
                incident[u].add(w)

    while True:
        rows = np.nonzero(alive)[0]
        pairs = rows[cKDTree(position[rows]).query_pairs(r=merge_tol,output_type='ndarray')].reshape(-1,2)
        order = np.argsort(np.linalg.norm(position[pairs[:,0]] - position[pairs[:,1]],axis=1),kind='stable')
        nfused = 0
        for u, v in pairs[order].tolist():
            if alive[u] and alive[v] and np.linalg.norm(position[u] - position[v]) <= merge_tol:
                fusion = plan(u,v)
                if fusion is not None:
                    fuse(*fusion)
                    nfused += 1
        if nfused == 0:
            break

    live = [w for w in range(len(walls)) if alive_wall[w]]
    walls = np.asarray([walls[w] for w in live],dtype=int).reshape(-1,2)
    group = np.asarray([group[w] for w in live],dtype=int)
    used = np.unique(walls)
    voronoi_vertices = np.unique(position[used],axis=0)
    index_vertices = _PointIndex(voronoi_vertices)
    walls = index_vertices.lookup(position[walls]).reshape(-1,2)
    finite_ridges_new = walls[group == 0]
    if len(finite_ridges_new) == 0:
        raise ValueError('No cell wall with both ends inside the clipping polygon is left after the merge: the '
                         'polygon is too small or too narrow for merge_tol = {}.'.format(merge_tol))
    boundary_ridges_new = np.vstack((walls[group == 1],walls[group == 2]))
    voronoi_ridges = np.vstack((finite_ridges_new,boundary_ridges_new))
    _plot_segments(voronoi_vertices[voronoi_ridges],'g',linestyle='--')

    on_boundary = used[kind[used] > 0]
    boundary_points = np.unique(position[on_boundary],axis=0)
    nvertex = voronoi_vertices.shape[0]
    # vertices on the polygon's perimeter (with boundaryFlag, the ends of the walls along it), the others inside
    nboundary_pts = np.unique(walls[group == 2]).shape[0] if flag else boundary_points.shape[0]
    nvertices_in = nvertex - nboundary_pts
    return voronoi_vertices,finite_ridges,boundary_points,finite_ridges_new,\
           boundary_ridges_new,nvertex,nvertices_in,finite_ridges_new.shape[0],boundary_ridges_new.shape[0],\
           nboundary_pts,boundary_points_featured.shape[0],voronoi_ridges,voronoi_ridges.shape[0]


def RebuildVoronoi_merge(vor,circles,boundaries,generation_center,x_min,x_max,y_min,y_max,box_center,box_shape,merge_tol,boundaryFlag):
    """Clip Voronoi mesh by the boundaries, merge short Voronoi ridges, and rebuild the new Voronoi mesh"""
    inside = _inside_clip_box(box_shape,boundaries,x_min,x_max,y_min,y_max,box_center)
    if isinstance(inside,_PolygonInside):
        return _rebuild_voronoi_merge_polygon(vor,circles,boundaries,generation_center,inside,merge_tol,boundaryFlag)
    [voronoi_vertices_in,boundary_points,boundary_points_featured,nboundary_pts_ini,\
     finite_ridges,finite_ridges_new,infinite_ridges_new] = \
        _clip_voronoi(vor,circles,boundaries,generation_center,inside,boundaryFlag)

    nvertices_in = voronoi_vertices_in.shape[0]
    nboundary_pts = boundary_points.shape[0]
    nboundary_pts_featured = boundary_points_featured.shape[0]
    nfinite_ridge = finite_ridges.shape[0]
    
    # Merge operation
    # Case 1 merge
    voronoi_vertices_in_whole = np.copy(voronoi_vertices_in)
    index_in_whole = _PointIndex(voronoi_vertices_in_whole) # coordinates -> first row in voronoi_vertices_in_whole
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
            fused_indices = nvertices_in + np.reshape(unique_inverse,-1)
            voronoi_vertices_in_fused = unique_pts.astype(float)
        else:
            # Offset the indices of boundary points as we inserted points in vertice_in list
            infinite_ridges_merged[infinite_ridges_merged >= nvertices_in] += len(rows_to_fuse_voronoi_vertices_in_whole)
            # Replace the indices of two close points with the index of new merged point in both finite_ridges_new and infinite_ridges_new lists
            fused_indices = nvertices_in + np.arange(len(rows_to_fuse_voronoi_vertices_in_whole))
        finite_ridges_merged, infinite_ridges_merged = \
            _replace_merged_indices([finite_ridges_merged,infinite_ridges_merged],rows_to_fuse_voronoi_vertices_in_whole,fused_indices,0,nvertices_in)
                
        voronoi_vertices_in_whole = np.vstack((voronoi_vertices_in_whole,voronoi_vertices_in_fused)) # create a whole vertex list
        index_in_whole.extend(voronoi_vertices_in_fused)
        voronoi_vertices_in_new = np.vstack((np.delete(voronoi_vertices_in,rows_to_fuse_voronoi_vertices_in,0),voronoi_vertices_in_fused)) # create a new vertex list
        # update the voronoi_vertices_in
        voronoi_vertices_in = np.copy(voronoi_vertices_in_new)
        nvertices_in = voronoi_vertices_in_whole.shape[0]
        # update the kd-tree
        tree_voronoi_vertices_in = cKDTree(voronoi_vertices_in)
        rows_to_fuse_voronoi_vertices_in = tree_voronoi_vertices_in.query_pairs(r=merge_tol,output_type='ndarray')
        # update the rows_to_fuse_voronoi_vertices_in with the rows-to-fuse in list voronoi_vertices_in_whole
        rows_to_fuse_voronoi_vertices_in_whole = np.copy(rows_to_fuse_voronoi_vertices_in)
        rows_to_fuse_voronoi_vertices_in_whole[:] = index_in_whole.lookup(voronoi_vertices_in_new[rows_to_fuse_voronoi_vertices_in])

    # delete rows with identical indices and reverse indices in finite_ridges_merged
    zero_rows = np.where((finite_ridges_merged[:,0] - finite_ridges_merged[:,1]) == 0)[0]
    finite_ridges_merged = np.delete(finite_ridges_merged,zero_rows,0)
    
    # delete rows with identical indices and reverse indices in infinite_ridges_merged
    zero_rows = np.where((infinite_ridges_merged[:,0] - infinite_ridges_merged[:,1]) == 0)[0]
    infinite_ridges_merged = np.delete(infinite_ridges_merged,zero_rows,0)
    
    if _is_on(boundaryFlag): 
        # construct the connectivity for a line path consisting of the boundary points
        boundary_ridges_new = _boundary_ridges_along_perimeter(boundary_points,boundaries) + nvertices_in # shift the indices with "nvertices_in"
        
        # Case 2 merge
        tree_boundary_points = cKDTree(boundary_points)
        rows_to_fuse_boundary_points = tree_boundary_points.query_pairs(r=merge_tol,output_type='ndarray')
        rows_to_fuse_boundary_points_whole = np.copy(rows_to_fuse_boundary_points)
        boundary_ridges_merged = np.copy(boundary_ridges_new)
        nboundary_pts = boundary_points.shape[0]
        boundary_points_whole = np.copy(boundary_points)
        index_boundary_whole = _PointIndex(boundary_points_whole) # coordinates -> first row in boundary_points_whole
        
        while len(rows_to_fuse_boundary_points) != 0:
            points_new = np.reshape(boundary_points_whole[rows_to_fuse_boundary_points_whole],(-1,2))
            boundary_points_fused = 0.5*(points_new[0::2] + points_new[1::2])
            
            # boundary_points_fused = boundary_points_fused.round(decimals=0)
            boundary_points_fused = np.round(boundary_points_fused/merge_tol).astype(int)*merge_tol # here we have to round to the nearest merge_tol number to avoid infinite loop
            
            # check if there are identical merged points from different pairs of closeby points (usually happens when three points are clustered)
            if len(np.unique(boundary_points_fused,axis=0)) != len(boundary_points_fused):
                [unique_pts,unique_inverse] = np.unique(boundary_points_fused,return_inverse=True,axis=0)
                # Replace the indices of closeby points with the index of new merged point in both finite_ridges_new and infinite_ridges_new lists
                fused_indices = nvertices_in + nboundary_pts + np.reshape(unique_inverse,-1)
                boundary_points_fused = unique_pts.astype(float)
            else:
                # Replace the indices of two close points with the index of new merged point in both boundary_ridges_new and infinite_ridges_new lists
                fused_indices = nvertices_in + nboundary_pts + np.arange(len(rows_to_fuse_boundary_points))
            boundary_ridges_merged, infinite_ridges_merged = \
                _replace_merged_indices([boundary_ridges_merged,infinite_ridges_merged],rows_to_fuse_boundary_points_whole+nvertices_in,\
                                        fused_indices,nvertices_in,nvertices_in+nboundary_pts)
                    
            boundary_points_whole = np.vstack((boundary_points_whole,boundary_points_fused)) # create a whole vertex list        
            index_boundary_whole.extend(boundary_points_fused)
            boundary_points_new = np.vstack((np.delete(boundary_points,rows_to_fuse_boundary_points,0),boundary_points_fused)) # create a new vertex list
            
            
            # update the boundary_points
            boundary_points = np.copy(boundary_points_new)
            nboundary_pts = boundary_points_whole.shape[0]
            # update the kd-tree
            tree_boundary_points = cKDTree(boundary_points)
            rows_to_fuse_boundary_points = tree_boundary_points.query_pairs(r=merge_tol,output_type='ndarray')
            # update the rows_to_fuse in list voronoi_vertices_in with the rows-to-fuse in list voronoi_vertices_in_whole
            rows_to_fuse_boundary_points_whole = np.copy(rows_to_fuse_boundary_points)
            rows_to_fuse_boundary_points_whole[:] = index_boundary_whole.lookup(boundary_points_new[rows_to_fuse_boundary_points])
        
                    
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
                rows_to_replace_voronoi_vertices_in_whole = index_in_whole.lookup(voronoi_vertices_in[rows_to_replace_voronoi_vertices_in])
                # replace indices
                boundary_point_index = index_boundary_whole.lookup(boundary_points[i,:]) + nvertices_in
                for rows in rows_to_replace_voronoi_vertices_in_whole:
                    finite_ridges_merged[finite_ridges_merged == rows] = boundary_point_index
                    infinite_ridges_merged[infinite_ridges_merged == rows] = boundary_point_index

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
    merged_ridges = voronoi_vertices_whole[all_ridges_merged]
    _plot_segments(merged_ridges,'g',linestyle='--')
    
    
    voronoi_vertices = np.unique(np.reshape(voronoi_vertices_whole[all_ridges_merged],(-1,2)),axis=0)
    # update the indices from voronoi_vertices_whole to voronoi_vertices(which is the list with no useless points)
    index_vertices = _PointIndex(voronoi_vertices)
    finite_ridges_merged[:] = index_vertices.lookup(voronoi_vertices_whole[finite_ridges_merged])
    infinite_ridges_merged[:] = index_vertices.lookup(voronoi_vertices_whole[infinite_ridges_merged])
    if _is_on(boundaryFlag):
        boundary_ridges_merged[:] = index_vertices.lookup(voronoi_vertices_whole[boundary_ridges_merged])
    # merged points are rounded to multiples of merge_tol, which can push them slightly
    # across slanted box edges (triangle/hexagon): move them back onto the boundary
    if box_shape in SQUARE_SHAPES + TRIANGLE_SHAPES + HEXAGON_SHAPES:
        voronoi_vertices = _project_into_convex_box(voronoi_vertices,boundaries)
    # delete rows with identical indices and reverse indices in finite_ridges_merged
    zero_rows = np.where((finite_ridges_merged[:,0] - finite_ridges_merged[:,1]) == 0)[0]
    finite_ridges_merged = np.delete(finite_ridges_merged,zero_rows,0)
    # delete rows with identical indices and reverse indices in infinite_ridges_merged
    zero_rows = np.where((infinite_ridges_merged[:,0] - infinite_ridges_merged[:,1]) == 0)[0]
    infinite_ridges_merged = np.delete(infinite_ridges_merged,zero_rows,0)
    
    if _is_on(boundaryFlag): 
        # delete rows with identical indices and reverse indices in boundary_ridges_merged
        zero_rows = np.where((boundary_ridges_merged[:,0] - boundary_ridges_merged[:,1]) == 0)[0]
        boundary_ridges_merged = np.delete(boundary_ridges_merged,zero_rows,0)
    
    # voronoi_vertices = voronoi_vertices_whole
    nvertex = voronoi_vertices.shape[0]
    finite_ridges_new = finite_ridges_merged
    nfinite_ridge = finite_ridges_new.shape[0]
    if _is_on(boundaryFlag): 
        boundary_ridges_new = np.vstack((infinite_ridges_merged,boundary_ridges_merged))
        nboundary_ridge = boundary_ridges_new.shape[0]
        voronoi_ridges = np.vstack((finite_ridges_new,boundary_ridges_new))
        nridge = voronoi_ridges.shape[0]
        
        nboundary_pts = np.unique(boundary_ridges_merged).shape[0] # vertices on the box perimeter
    else:
        boundary_ridges_new = np.copy(infinite_ridges_merged)  
        nboundary_ridge = boundary_ridges_new.shape[0]
        voronoi_ridges = np.vstack((finite_ridges_new,boundary_ridges_new))
        nridge = voronoi_ridges.shape[0]
        # nboundary_pts is the number of cut points
    # all other vertices are inside the box. Not all of them belong to a finite ridge: a vertex
    # near a box corner whose neighbouring vertices are all outside the box has only cut ridges
    nvertices_in = nvertex - nboundary_pts

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
    if npt_per_layer != nvertex:
        raise ValueError('Inconsistent vertex counts: {:d} vertices inside the box + {:d} on its boundary != {:d} vertices.'.format(nvertices_in,nboundary_pts,nvertex))
    npt_per_layer_normal = npt_per_layer - nboundary_pts_featured
    npt_per_layer_vtk = nfinitept_per_layer + nboundary_pt_per_layer
    
    all_midpt_indices = np.vstack((np.reshape(finite_midpt_indices,(2,-1)).T,np.reshape(boundary_midpt_indices,(2,-1)).T))
    all_quarterpt_indices = np.vstack((np.reshape(finite_quarterpt_indices,(2,-1)).T,np.reshape(boundary_quarterpt_indices,(2,-1)).T))
    
    all_pts_2D = np.vstack((voronoi_vertices,finite_ridge_mid,finite_ridge_quarter,boundary_ridge_mid,boundary_ridge_quarter))
    all_ridges = np.hstack((voronoi_ridges,all_midpt_indices,all_quarterpt_indices))
    
    return all_pts_2D,all_ridges,npt_per_layer,npt_per_layer_normal,npt_per_layer_vtk



def VertexandRidgeinfo(all_pts_2D,all_ridges,npt_per_layer,npt_per_layer_normal,\
                       npt_per_layer_vtk,nridge,geoName,radii,generation_center,\
                       cellwallthickness_sparse,cellwallthickness_dense,path='meshes'):
    
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
        if (bisect.bisect(radii,vertex_distance2logcenter[i]) % 2) == 0: # bisect = ring index + 1: even, sparse rings
            vertex_cellwallthickness_2D[i] = cellwallthickness_sparse
        else: # if even, dense cells
            vertex_cellwallthickness_2D[i] = cellwallthickness_dense
#==================================================================================================
     
    # Find the wings (Voronoi ridges) of all vertices: the entries of all_ridges equal to
    # each vertex index, in row-major order (as np.where(all_ridges == i) gives them)
    all_ridges = np.asarray(all_ridges)
    flat_ridges = all_ridges.reshape(-1)
    order = np.argsort(flat_ridges,kind='stable')
    first, last = np.searchsorted(flat_ridges[order],[0,npt_per_layer])
    wing_vertex = flat_ridges[order[first:last]] # vertex of each wing, wings sorted by vertex
    wing_ridge = order[first:last] // all_ridges.shape[1] # ridge index of each wing
    wing_position = order[first:last] % all_ridges.shape[1] # position of the vertex in the ridge, p1 --> 0, p2 --> 1
    if np.any(wing_position > 1) or np.any(all_ridges[wing_ridge,0] == all_ridges[wing_ridge,1]):
        raise ValueError('Inconsistent Voronoi ridges: a vertex index is also used by a ridge mid/quarter point, or a ridge has zero length.')
    nwings = np.bincount(wing_vertex,minlength=npt_per_layer) # total number of wings for each vertex

    # flattened 1D numpy array of ridge indices, indices of the vertices at the other end of each ridge,
    # lengths, widths, angles and directions of wings
    direction = -np.power(-1,wing_position+1) # convert the way to express the direction of each ridge to: 1/-1 --> pointing from/to the current vertex
    ridge_inices = wing_ridge
    another_vertex_inices = all_ridges[wing_ridge,1-wing_position]
    ridge_lengths = all_ridge_lengths[wing_ridge]/2
    ridge_widths = vertex_cellwallthickness_2D[wing_vertex]
    ridge_angles = all_ridge_angles[wing_ridge]+math.pi*(-direction.clip(max=0))
    ridge_directions = direction
    
    flattened_all_vertices_2D = np.vstack((ridge_inices,another_vertex_inices,ridge_lengths,ridge_widths,ridge_angles,ridge_directions))
 
    # Convert ragged info list to a rectangular shape numpy array
    max_wings = int(nwings.max())
    
    all_vertices_info_2D_nparray = np.zeros((npt_per_layer,1+max_wings*5))
    all_vertices_info_2D_nparray[:,0] = nwings
    
    iwing = np.arange(len(wing_vertex)) - np.repeat(np.cumsum(nwings)-nwings,nwings) # position of each wing among the wings of its vertex
    for icol, values in enumerate([ridge_inices,another_vertex_inices,ridge_lengths,ridge_widths,ridge_angles]):
        all_vertices_info_2D_nparray[wing_vertex,icol*max_wings+1+iwing] = values
    
    # Save info to txt files
    all_vertices_2D = np.hstack((all_pts_2D[0:npt_per_layer,:],all_vertices_info_2D_nparray))
    _savetxt(Path(path, geoName, geoName +'-vertex.mesh'), all_vertices_2D, fmt='%.16g', delimiter=' '\
        ,header='Vertex Data Generated with RingsPy Mesh Generation Tool\n\
Number of vertices\n'+ str(npt_per_layer) +
    '\n\
Max number of wings for one vertex\n'+ str(max_wings) +
    '\n\
[xcoord ycoord nwings ridge1 ... ridgen farvertex1 ... farvertexn length1 ... lengthn width1 ... widthn angle1 ... anglen]', comments='')
    
    _savetxt(Path(path, geoName, geoName +'-ridge.mesh'), all_ridges, fmt='%d', delimiter=' '\
        ,header='Ridge Data Generated with RingsPy Mesh Generation Tool\n\
Number of ridges\n'+ str(nridge) +
    '\n\
[vertex1 vertex2 midpt1 midpt2 qrtrpt1 qrtrpt2]', comments='')
    
    return all_vertices_2D, max_wings, flattened_all_vertices_2D, all_ridges


def GenerateBeamElement(NURBS_degree,nsegments,theta_min,theta_max,\
                        z_min,z_max,long_connector_ratio,npt_per_layer,voronoi_vertices,\
                        nvertex,voronoi_ridges,nridge,generation_center,all_vertices_2D,max_wings,\
                        flattened_all_vertices_2D,all_ridges):
    

    nctrlpt_per_elem = NURBS_degree + 1
    nctrlpt_per_beam = 2*NURBS_degree + 1
    
    nlayers = nctrlpt_per_beam*nsegments
    segment_length = (z_max - z_min) / nsegments
    
    nconnector_t_per_beam = int((nctrlpt_per_beam-1)/NURBS_degree+1)
    nconnector_t_per_grain = int(nconnector_t_per_beam*nsegments)
    
    theta = np.linspace(theta_min,theta_max,nlayers-(nsegments-1))
    z_coord = np.linspace(z_min,z_max,nlayers-(nsegments-1))
    connector_l_length = segment_length*long_connector_ratio
    connector_l_angle = (theta_max-theta_min)/nsegments*long_connector_ratio
    
    # Insert repeated layers for the longitudinal connectors
    for i in range(nlayers-(nsegments-1)-(nctrlpt_per_beam-1),1,-(nctrlpt_per_beam-1)): 
        theta = np.insert(theta,i,theta[i-1])
        theta[i:] += connector_l_angle
    for i in range(nlayers-(nsegments-1)-(nctrlpt_per_beam-1),1,-(nctrlpt_per_beam-1)):
        z_coord = np.insert(z_coord,i,z_coord[i-1])
        z_coord[i:] += connector_l_length
    
    # Data for IGA use
    vertices_new = _rotate_layers(voronoi_vertices[np.arange(npt_per_layer)],theta[0:nlayers],z_coord[0:nlayers],generation_center)
    
    # Vertex Data for IGA
    IGAvertices = np.reshape(vertices_new,(-1,3))
    
    # Connectivity for IGA
    npt_total = nlayers*npt_per_layer
    vertex_connectivity = np.linspace(0,npt_total-1,npt_total)
    
    # Beams
    ngrain = nvertex
    nbeam_total = ngrain*nsegments
    i, j = np.meshgrid(np.arange(nsegments),np.arange(ngrain),indexing='ij')
    ivertex = (i*npt_per_layer*nctrlpt_per_beam + j).reshape(-1,1)
    beam_connectivity_original = (ivertex + np.arange(nctrlpt_per_beam)*npt_per_layer).astype(float)
    
    # Rearange beam connectivity such that each row is corresponding to a Bezier beam element
    beam_connectivity = np.copy(beam_connectivity_original)
    for ictrlpt in range(nctrlpt_per_beam-NURBS_degree,1,-NURBS_degree):
        beam_connectivity = np.insert(beam_connectivity,ictrlpt,beam_connectivity[:,ictrlpt-1],axis=1)
    
    beam_connectivity_original = (beam_connectivity_original+1).astype(int) # +1 because of in abaqus index starts from 1
    beam_connectivity = (np.reshape(beam_connectivity,(-1,nctrlpt_per_elem))+1).astype(int) # +1 because of in abaqus index starts from 1
    
    
    nbeamElem = beam_connectivity.shape[0]
    
    
    # Transverse connectors 
    nconnector_t = nridge*nconnector_t_per_grain
    i, j = np.meshgrid(np.arange(nsegments),np.arange(nconnector_t_per_beam),indexing='ij')
    layer_offset = (i*npt_per_layer*nctrlpt_per_beam + j*NURBS_degree*npt_per_layer).reshape(-1,1,1)
    connector_t_connectivity = np.zeros((nconnector_t,2))
    connector_t_connectivity[:] = (np.asarray(voronoi_ridges)[np.newaxis,0:nridge,:] + layer_offset).reshape(-1,2)
                
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
    nconnector_l = nconnector_l_per_layer * (nsegments-1)
    connector_l_connectivity = np.zeros((nconnector_l,2))
    connector_l_vertex_dict = np.zeros(nconnector_l)
    # one connector per wing of each grain, for each layer of longitudinal connectors
    irow = (np.arange(nsegments-1)[:,np.newaxis]*ngrain + np.arange(ngrain)[np.newaxis,:]).reshape(-1)
    irow = np.repeat(irow,np.tile(nwings[0:ngrain],nsegments-1))
    connector_l_connectivity[0:len(irow)] = np.column_stack((beam_connectivity_original[irow,-1],beam_connectivity_original[irow,-1]+npt_per_layer))
    connector_l_vertex_dict[0:len(irow)] = irow
    
    connector_l_vertex_dict = connector_l_vertex_dict.astype(int)
    connector_l_connectivity = connector_l_connectivity.astype(int)
    
    nconnector_total = nconnector_t + nconnector_l
    
    return IGAvertices,vertex_connectivity,beam_connectivity_original,nbeam_total,\
        beam_connectivity,nbeamElem,nctrlpt_per_beam,nlayers,segment_length,connector_t_connectivity,\
        connector_t_bot_connectivity,connector_t_top_connectivity,\
        connector_t_reg_connectivity,connector_l_connectivity,nconnector_t_per_beam,\
        nconnector_t_per_grain,nconnector_t,nconnector_l,nconnector_total,\
        theta,z_coord,connector_l_vertex_dict


def ConnectorMeshFile(geoName,IGAvertices,connector_t_bot_connectivity,\
                      connector_t_reg_connectivity,connector_t_top_connectivity,\
                      connector_l_connectivity,all_vertices_2D,\
                      max_wings,flattened_all_vertices_2D,nsegments,segment_length,\
                      nctrlpt_per_beam,theta,nridge,connector_l_vertex_dict,path='meshes'):
######### txt File stores the connector data for Abaqus analyses ##############
    nelem_connector_t_bot = connector_t_bot_connectivity.shape[0]
    nelem_connector_t_reg = connector_t_reg_connectivity.shape[0]
    nelem_connector_t_top = connector_t_top_connectivity.shape[0]
    nelem_connector_l = connector_l_connectivity.shape[0]
    nelem_total = nelem_connector_t_bot + nelem_connector_t_reg + nelem_connector_t_top + nelem_connector_l
    
    # calculate heights of transverse connectors
    height_connector_t = segment_length/4
    
    # calculate longitudinal connectors
    nlayers_conn_l = nsegments - 1
    conn_l_ridge_index = flattened_all_vertices_2D[0,:].astype(int)
    conn_l_lengths = flattened_all_vertices_2D[2,:]
    # conn_l_widths = flattened_all_vertices_2D[3,:]  
    conn_l_angles = flattened_all_vertices_2D[4,:]
    
    conn_l_ridge_index = np.tile(conn_l_ridge_index,nlayers_conn_l)
    conn_l_lengths = np.tile(conn_l_lengths,nlayers_conn_l)
    conn_l_angles = np.tile(conn_l_angles,nlayers_conn_l)
    
    for i in range(0,nlayers_conn_l):
        conn_l_angles[i*2*nridge:(i+1)*2*nridge] = conn_l_angles[i*2*nridge:(i+1)*2*nridge] - theta[(i+1)*nctrlpt_per_beam]
        
    rot = np.array([[np.cos(conn_l_angles), -np.sin(conn_l_angles)], [np.sin(conn_l_angles), np.cos(conn_l_angles)]]).T
    rot = np.hstack((rot[:,0,:],np.zeros(len(conn_l_angles))[:, np.newaxis])) # convert to 3D rotation matrix, assumes rotation remains still in-plane
    conn_l_tangents = rot
    
# Meshdata = [nodex1 nodey1 nodez1 nodex2 nodey2 nodez2 centerx centery centerz dx1 dy1 dz1 dx2 dy2 dz2 width height]    
    Meshdata = np.zeros((nelem_total,17))
    
    connectivity = np.concatenate((connector_t_bot_connectivity,connector_t_reg_connectivity,connector_t_top_connectivity,connector_l_connectivity))
    Meshdata[:,0:3] = IGAvertices[connectivity[:,0]-1,:]
    Meshdata[:,3:6] = IGAvertices[connectivity[:,1]-1,:]
    ielem_reg = nelem_connector_t_bot
    ielem_top = ielem_reg + nelem_connector_t_reg
    ielem_l = ielem_top + nelem_connector_t_top
    Meshdata[0:ielem_reg,16] = height_connector_t
    Meshdata[ielem_reg:ielem_top,16] = height_connector_t*2
    Meshdata[ielem_top:ielem_l,16] = height_connector_t
    Meshdata[ielem_l:nelem_total,16] = conn_l_lengths[0:nelem_connector_l]
        
    Meshdata[:,6:9] = (Meshdata[:,0:3] + Meshdata[:,3:6])/2
    # Meshdata[0:nelem_connector_t_bot,8] = Meshdata[0:nelem_connector_t_bot,8]
    # Meshdata[nelem_total-nelem_connector_t_top:nelem_total,8] = Meshdata[nelem_total-nelem_connector_t_top:nelem_total,8]
    Meshdata[:,9:12] = Meshdata[:,6:9] - Meshdata[:,0:3]
    Meshdata[:,12:15] = Meshdata[:,6:9] - Meshdata[:,3:6]

    # Add z-variation/random field for cellwall thicknesses/connector widths
    Meshdata[0:nridge,15] = all_vertices_2D[connector_t_bot_connectivity[0:nridge,1]-1,3+max_wings*3] # assign widths to all bot transverse connectors
    Meshdata[:nelem_total-nelem_connector_l,15] = np.tile(Meshdata[0:nridge,15],3*nsegments) # use the same widths for reg and top transverse connectors
    
    conn_l_widths = Meshdata[conn_l_ridge_index,15]
    Meshdata[nelem_total-nelem_connector_l:nelem_total,15] = conn_l_widths

    # add the eccentricity to the centers for longitudinal connectors (new center is correct)
    Meshdata[nelem_total-nelem_connector_l:nelem_total,6:9] += conn_l_tangents*Meshdata[nelem_total-nelem_connector_l:nelem_total,16][:, np.newaxis]/2

    # Replace nodal coordinates with nodal indices
    Meshdata[:,0:2] = connectivity
    Meshdata = np.delete(Meshdata,[2,3,4,5],1)
    
    _savetxt(Path(path, geoName, geoName+'-mesh.txt'), Meshdata, fmt='%.16g', delimiter=' '\
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

    return Meshdata, conn_l_tangents


def insert_precracks(all_pts_2D,all_ridges,nridge,npt_per_layer,\
                     npt_per_layer_normal,npt_per_layer_vtk,\
                     nlayers,precrack_nodes,precrack_widths,\
                     cellsize_sparse,):
    """
    Find the cell walls (Voronoi ridges) crossed by the precrack line segments
    precrack_nodes = [[x1 y1 x2 y2], ...]. Returns the ridge indices; the
    transverse connectors of these ridges are precracked in all segments.
    """
    precrack_midpts = (precrack_nodes[:,0:2]+precrack_nodes[:,2:4])/2.0
    precrack_lengths = np.linalg.norm(precrack_nodes[:,2:4]-precrack_nodes[:,0:2],axis=1)
    ridge_midpts = all_pts_2D[all_ridges[:,2]]
    ridge_midpts_tree = cKDTree(ridge_midpts)
    near_ridges = []
    for midpt, length in zip(precrack_midpts,precrack_lengths):
        near_ridges.append(ridge_midpts_tree.query_ball_point(midpt,length/2+3*cellsize_sparse))

    # Find the ridges (line segments) intersecting the precrack line segments
    precracked_elem = []
    for i in range(0,len(precrack_nodes)):
        p = precrack_nodes[i,0:2]
        normal = (precrack_nodes[i,2:4] - p)/precrack_lengths[i]
        for ridge in near_ridges[i]:
            q = all_pts_2D[all_ridges[ridge,0],:]
            s = all_pts_2D[all_ridges[ridge,1],:] - q
            denom = cross2D(normal,s)
            if denom == 0: # parallel
                continue
            t = cross2D((q-p),s)/denom # position along the precrack
            u = cross2D((q-p),normal)/denom # position along the ridge
            if (u >= 0) and (u <= 1) and (t >= 0) and (t <= precrack_lengths[i]) and (ridge not in precracked_elem):
                precracked_elem.append(ridge)

    nconnector_t_precrack = len(precracked_elem) # number of precracked cell walls (ridges)
    nconnector_l_precrack = 0

    # Visualize the precracks in the preview plot
    plt.plot(precrack_nodes[0,0::2],precrack_nodes[0,1::2],'r-',linewidth=2)
    
    return precracked_elem, nconnector_t_precrack, nconnector_l_precrack
        
        
def _rotate_layers(points,theta,z,origin):
    """
    2D points rotated around origin by each angle in theta, with z coordinates z:
    array of shape (len(theta), len(points), 3). Same values as calling
    rotate_around_point_highperf for each point and angle.
    """
    offset_x, offset_y = origin
    adjusted_x = points[:,0] - offset_x
    adjusted_y = points[:,1] - offset_y
    layers = np.zeros((len(theta),len(points),3))
    for i, radians in enumerate(theta):
        cos_rad = math.cos(radians)
        sin_rad = math.sin(radians)
        layers[i,:,0] = offset_x + cos_rad * adjusted_x + sin_rad * adjusted_y
        layers[i,:,1] = offset_y + -sin_rad * adjusted_x + cos_rad * adjusted_y
        layers[i,:,2] = z[i]
    return layers


def _format_rows(array,row_format,join=True,chunk_rows=20000):
    """
    Text of a 1D or 2D array, one row_format line per row: the same as
    ''.join(row_format % tuple(row) for row in array), but formatted in chunks of
    rows with one % operation each (much faster). join=False returns the chunks.
    """
    array = np.asarray(array)
    ncol = 1 if array.ndim == 1 else array.shape[1]
    values = array.reshape(-1).tolist()
    chunks = []
    for start in range(0,len(values),chunk_rows*ncol):
        part = values[start:start+chunk_rows*ncol]
        chunks.append((row_format*(len(part)//ncol)) % tuple(part))
    return ''.join(chunks) if join else chunks


def _savetxt(fname,X,fmt='%.18e',delimiter=' ',header='',comments='# '):
    """np.savetxt for a 2D array and a single format, writing the same text but faster"""
    X = np.asarray(X)
    with open(fname,'w') as txtfile:
        if len(header) > 0:
            txtfile.write(comments + header.replace('\n','\n' + comments) + '\n')
        txtfile.writelines(_format_rows(X,delimiter.join([fmt]*X.shape[1]) + '\n',join=False))


def _write_vtu_appended(filename,points,cells,types,point_data=None,cell_data=None):
    """
    Write an unstructured grid (.vtu) file with all data arrays stored as raw binary
    appended data. cells is a list of connectivity arrays (one row per cell),
    point_data and cell_data are (name, values) pairs or None.
    """
    connectivity = np.concatenate([np.reshape(c,-1) for c in cells])
    offsets = np.cumsum(np.concatenate([np.full(len(c),np.shape(c)[1]) for c in cells]))
    int_type, int_dtype = ('Int32','<i4') if max(connectivity.max(initial=0),offsets.max(initial=0)) < 2**31 else ('Int64','<i8')
    blocks = []
    def appended(data):
        """Add a data block, return its offset in the appended data"""
        offset = sum(len(block) for block in blocks)
        data = np.ascontiguousarray(data).tobytes()
        blocks.append(np.array(len(data),dtype='<u8').tobytes() + data)
        return offset

    xml = ['<?xml version="1.0"?>\n',
           '<VTKFile type="UnstructuredGrid" version="1.0" byte_order="LittleEndian" header_type="UInt64">\n',
           '<UnstructuredGrid>\n',
           '<Piece NumberOfPoints="{:d}" NumberOfCells="{:d}">\n'.format(len(points),len(types))]
    for tag, data in [('PointData',point_data),('CellData',cell_data)]:
        if data is not None:
            name, values = data
            xml.append('<{} Scalars="{}">\n'.format(tag,name))
            xml.append('<DataArray type="Float32" Name="{}" format="appended" offset="{:d}"/>\n'.format(name,appended(np.asarray(values,dtype='<f4'))))
            xml.append('</{}>\n'.format(tag))
    xml.append('<Points>\n')
    xml.append('<DataArray type="Float64" NumberOfComponents="3" format="appended" offset="{:d}"/>\n'.format(appended(np.asarray(points,dtype='<f8'))))
    xml.append('</Points>\n')
    xml.append('<Cells>\n')
    xml.append('<DataArray type="{}" Name="connectivity" format="appended" offset="{:d}"/>\n'.format(int_type,appended(connectivity.astype(int_dtype))))
    xml.append('<DataArray type="{}" Name="offsets" format="appended" offset="{:d}"/>\n'.format(int_type,appended(offsets.astype(int_dtype))))
    xml.append('<DataArray type="UInt8" Name="types" format="appended" offset="{:d}"/>\n'.format(appended(np.asarray(types,dtype='u1'))))
    xml.append('</Cells>\n')
    xml.append('</Piece>\n')
    xml.append('</UnstructuredGrid>\n')
    xml.append('<AppendedData encoding="raw">\n_')
    with open(filename,'wb') as vtkfile:
        vtkfile.write(''.join(xml).encode('ascii'))
        for block in blocks:
            vtkfile.write(block)
        vtkfile.write(b'\n</AppendedData>\n</VTKFile>\n')


def VisualizationFiles(geoName,NURBS_degree,nlayers,npt_per_layer_vtk,all_pts_2D,\
                       segment_length,theta,z_coord,nsegments,nridge,\
                       voronoi_ridges,generation_center,all_ridges,nvertex,\
                       nconnector_t,nconnector_l,nctrlpt_per_beam,ConnMeshData,\
                       conn_l_tangents,all_vertices_2D,max_wings,\
                       flattened_all_vertices_2D,path='meshes',vtk_format='ascii'):
    """
    Write the Paraview files geoName_vertices.vtu, geoName_beams.vtu, geoName_conns.vtu
    and geoName_conns_vol.vtu. vtk_format = 'ascii' (default) writes text files,
    'binary' writes the data arrays as raw binary appended data (faster to write and
    read, about 3 times smaller).
    """
    if vtk_format not in ['ascii','binary']:
        raise ValueError("vtk_format must be 'ascii' or 'binary', not {!r}".format(vtk_format))
    
    # Calculate model parameters
    ngrain = nvertex
    ninterval_per_beam_vtk = int((nctrlpt_per_beam-1)/2) # 2 layers ----> 1 interval
    nconnector_t_per_beam = int((nctrlpt_per_beam-1)/NURBS_degree+1)
    
    # Vertex Data for VTK use
    vertices_new = _rotate_layers(all_pts_2D[np.arange(npt_per_layer_vtk)],theta[0:nlayers],z_coord[0:nlayers],generation_center)
    VTKvertices = np.reshape(vertices_new,(-1,3))
    
    # Connectivity for VTK use
    # Vertex
    npt_total_vtk = nlayers*npt_per_layer_vtk
    vertex_connectivity_vtk = np.arange(npt_total_vtk)
    
    # Beam (3 points per beam in Paraview)
    nbeam_total_vtk = ngrain*nsegments*ninterval_per_beam_vtk
    i, j, k = np.meshgrid(np.arange(nsegments),np.arange(ninterval_per_beam_vtk),np.arange(ngrain),indexing='ij')
    ivertex = ((i + (i*ninterval_per_beam_vtk+j)*2)*npt_per_layer_vtk + k).reshape(-1)
    beam_connectivity_vtk = np.column_stack((ivertex,ivertex+2*npt_per_layer_vtk,ivertex+npt_per_layer_vtk))
    
    # Transverse connectors 
    i, j = np.meshgrid(np.arange(nsegments),np.arange(nconnector_t_per_beam),indexing='ij')
    layer_offset = (i*npt_per_layer_vtk*nctrlpt_per_beam + j*NURBS_degree*npt_per_layer_vtk).reshape(-1,1,1)
    connector_t_connectivity_vtk = np.zeros((nconnector_t,2),dtype=np.int64)
    connector_t_connectivity_vtk[:] = (np.asarray(voronoi_ridges)[np.newaxis,0:nridge,0:2] + layer_offset).reshape(-1,2)
    
    # Longitudinal connectors: one per wing of each grain, for each layer of longitudinal connectors
    connector_l_connectivity_vtk = np.zeros((nconnector_l,2),dtype=np.int64)
    nwings = all_vertices_2D[:,2].astype(int)
    irow = (np.arange(nsegments-1)[:,np.newaxis]*ngrain + np.arange(ngrain)[np.newaxis,:]).reshape(-1)
    irow = np.repeat(irow,np.tile(nwings[0:ngrain],nsegments-1))
    beam_end = beam_connectivity_vtk[irow+(ninterval_per_beam_vtk-1)*ngrain,1]
    connector_l_connectivity_vtk[0:len(irow)] = np.column_stack((beam_end,beam_end+npt_per_layer_vtk))
    
    # Quad - Beam Wings
    nquad = nridge*2
    nquad_total = nquad*nsegments*ninterval_per_beam_vtk
    i, j = np.meshgrid(np.arange(nsegments),np.arange(ninterval_per_beam_vtk),indexing='ij')
    layer = (i + (i*ninterval_per_beam_vtk+j)*2).reshape(-1,1,1)*npt_per_layer_vtk # [segment and interval, wing, ridge]
    all_ridges = np.asarray(all_ridges)
    vertex = all_ridges[np.newaxis,0:nridge,0:2].transpose(0,2,1)
    midpt = all_ridges[np.newaxis,0:nridge,2:4].transpose(0,2,1)
    qrtrpt = all_ridges[np.newaxis,0:nridge,4:6].transpose(0,2,1)
    n = npt_per_layer_vtk
    quad_connectivity = np.stack((vertex+layer,vertex+layer+2*n,midpt+layer+2*n,midpt+layer,\
                                  vertex+layer+n,qrtrpt+layer+2*n,midpt+layer+n,qrtrpt+layer),axis=-1).reshape(nquad_total,8)

    # Quad - Connector Cross-sections
    nquad_conn = nconnector_t + nconnector_l
    nquad_conn_total = nquad_conn
    
    Quad_center = np.copy(ConnMeshData[:,2:5])
    
    ##### Set normals of longitudinal connectors to be vertical normal = (0,0,1)
    ConnMeshData[nconnector_t:nquad_conn_total,5:7] = 0
    ConnMeshData[nconnector_t:nquad_conn_total,8:10] = 0
    #### (Need to changed for inclined beam axis in the future)
    
    Quad_normal1 = np.copy(ConnMeshData[:,5:8])
    Quad_normal2 = np.copy(ConnMeshData[:,8:11])
    Quad_length1 = np.linalg.norm(Quad_normal1, axis=1)
    Quad_length2 = np.linalg.norm(Quad_normal2, axis=1)
    Quad_width = np.copy(ConnMeshData[:,11])
    Quad_height = np.copy(ConnMeshData[:,12])
    Quad_normal = Quad_normal1/Quad_length1[:, np.newaxis]
    Quad_tangent = np.zeros((nquad_conn_total,3))
    Quad_tangent[0:nconnector_t] = np.tile(np.array([0.0,0.0,1.0]),(nconnector_t,1)) # assume a tangent of (0,0,1) for all conn_t

    Quad_tangent[nconnector_t:nquad_conn] = conn_l_tangents
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
    
    quad_conn_connectivity = np.arange(nquad_conn_total*4)
    quad_conn_connectivity = np.reshape(quad_conn_connectivity, [4,-1]).T
    quad_conn_connectivity = quad_conn_connectivity + npt_total_vtk

    # Hex - Connector Volumes
    # vertices: quad_conn_vertices (the cross-sections), followed by the cross-sections
    # moved to both ends of the connector
    hex_conn_vertices = np.zeros((nquad_conn_total*4*2,3))
    hex_conn_vertices[0:nquad_conn_total,:] = Quad_center - Quad_tangent*Quad_height[:, np.newaxis]/2 - Quad_bitangent*Quad_width[:, np.newaxis]/2 - Quad_normal*Quad_length1[:, np.newaxis]
    hex_conn_vertices[nquad_conn_total:2*nquad_conn_total,:] = Quad_center - Quad_tangent*Quad_height[:, np.newaxis]/2 + Quad_bitangent*Quad_width[:, np.newaxis]/2 - Quad_normal*Quad_length1[:, np.newaxis]
    hex_conn_vertices[2*nquad_conn_total:3*nquad_conn_total,:] = Quad_center + Quad_tangent*Quad_height[:, np.newaxis]/2 + Quad_bitangent*Quad_width[:, np.newaxis]/2 - Quad_normal*Quad_length1[:, np.newaxis]
    hex_conn_vertices[3*nquad_conn_total:4*nquad_conn_total,:] = Quad_center + Quad_tangent*Quad_height[:, np.newaxis]/2 - Quad_bitangent*Quad_width[:, np.newaxis]/2 - Quad_normal*Quad_length1[:, np.newaxis]

    hex_conn_vertices[4*nquad_conn_total:5*nquad_conn_total,:] = Quad_center - Quad_tangent*Quad_height[:, np.newaxis]/2 - Quad_bitangent*Quad_width[:, np.newaxis]/2 + Quad_normal*Quad_length2[:, np.newaxis]
    hex_conn_vertices[5*nquad_conn_total:6*nquad_conn_total,:] = Quad_center - Quad_tangent*Quad_height[:, np.newaxis]/2 + Quad_bitangent*Quad_width[:, np.newaxis]/2 + Quad_normal*Quad_length2[:, np.newaxis]
    hex_conn_vertices[6*nquad_conn_total:7*nquad_conn_total,:] = Quad_center + Quad_tangent*Quad_height[:, np.newaxis]/2 + Quad_bitangent*Quad_width[:, np.newaxis]/2 + Quad_normal*Quad_length2[:, np.newaxis]
    hex_conn_vertices[7*nquad_conn_total:8*nquad_conn_total,:] = Quad_center + Quad_tangent*Quad_height[:, np.newaxis]/2 - Quad_bitangent*Quad_width[:, np.newaxis]/2 + Quad_normal*Quad_length2[:, np.newaxis]

    hex_conn1_connectivity = np.arange(nquad_conn_total*4*2)
    hex_conn1_connectivity = np.reshape(hex_conn1_connectivity, [8,-1]).T
    hex_conn2_connectivity = np.arange(nquad_conn_total*4*2)
    hex_conn2_connectivity = np.reshape(hex_conn2_connectivity, [8,-1]).T
    hex_conn2_connectivity[:,4:8] = hex_conn2_connectivity[:,4:8] + nquad_conn_total*4
    hex_conn_connectivity = np.vstack((hex_conn1_connectivity,hex_conn2_connectivity))
    hex_conn_connectivity = hex_conn_connectivity + npt_total_vtk

# =============================================================================
    # Paraview Visualization File
    collocation_flag_vtk = np.concatenate((np.ones(ngrain),np.zeros(npt_per_layer_vtk*NURBS_degree-ngrain)))
    collocation_flag_vtk = np.concatenate((np.tile(collocation_flag_vtk, nconnector_t_per_beam-1),np.concatenate((np.ones(ngrain),np.zeros(npt_per_layer_vtk-ngrain)))))
    collocation_flag_vtk = np.tile(collocation_flag_vtk, nsegments)

    VTKcell_types_vertices = (np.ones(npt_total_vtk)).astype(int)
    VTKcell_types_beams = np.concatenate((21*np.ones(nbeam_total_vtk),23*np.ones(nquad_total))).astype(int)
    VTKcell_types_conns = np.concatenate((3*np.ones(nconnector_t),3*np.ones(nconnector_l),9*np.ones(nconnector_t),9*np.ones(nconnector_l))).astype(int)
    VTKcell_types_conns_vol = np.concatenate((12*np.ones(nconnector_t),12*np.ones(nconnector_l),12*np.ones(nconnector_t),12*np.ones(nconnector_l))).astype(int)
    Quad_width_vtk = np.tile(Quad_width,(2))

    files = [ # file name, point arrays, connectivity arrays (one row per cell), cell types, point data, cell data
        ('_vertices',[VTKvertices],[vertex_connectivity_vtk.reshape(-1,1)],VTKcell_types_vertices,
         ('IGAcollocation_flag',collocation_flag_vtk),None),
        ('_beams',[VTKvertices],[beam_connectivity_vtk,quad_connectivity],VTKcell_types_beams,None,None),
        ('_conns',[VTKvertices,quad_conn_vertices],[connector_t_connectivity_vtk,connector_l_connectivity_vtk,quad_conn_connectivity],
         VTKcell_types_conns,None,('Connector_width',Quad_width_vtk)),
        ('_conns_vol',[VTKvertices,quad_conn_vertices,hex_conn_vertices],[hex_conn_connectivity],
         VTKcell_types_conns_vol,None,('Connector_width',Quad_width_vtk)),
        ]

    if vtk_format == 'binary':
        for suffix, points, cells, types, point_data, cell_data in files:
            _write_vtu_appended(Path(path, geoName, geoName + suffix + '.vtu'),np.vstack(points),cells,types,point_data,cell_data)
        return

    # ASCII files (the points shared by several files are formatted once)
    points_format = ' %11.8e  %11.8e  %11.8e\n'
    shared_points_text = {id(VTKvertices): _format_rows(VTKvertices,points_format),
                          id(quad_conn_vertices): _format_rows(quad_conn_vertices,points_format)}
    for suffix, points, cells, types, point_data, cell_data in files:
        with open(Path(path, geoName, geoName + suffix + '.vtu'),'w') as vtkfile:
            write = vtkfile.write
            write('<VTKFile type="UnstructuredGrid" version="2.0" byte_order="LittleEndian">\n')
            write('<UnstructuredGrid>\n')
            write('<Piece NumberOfPoints="'+str(sum(p.shape[0] for p in points))+'" NumberOfCells="'+str(len(types))+'">\n')

            # <Points>
            write('<Points>\n')
            write('<DataArray type="Float64" NumberOfComponents="3" format="ascii">\n')
            for p in points:
                if id(p) in shared_points_text:
                    write(shared_points_text[id(p)])
                else:
                    vtkfile.writelines(_format_rows(p,points_format,join=False))
            write('</DataArray>\n')
            write('</Points>\n')
            # </Points>

            # <PointData>
            if point_data is not None:
                name, values = point_data
                write('<PointData Tensors="" Vevtors="" Scalars="'+name+'">\n')
                write('<DataArray type="Float32" Name="'+name+'" format="ascii">\n')
                write(_format_rows(values,'%11.8e\n'))
                write('</DataArray>\n')
                write('</PointData>\n')
            # </PointData>

            # <Cells>
            write('<Cells>\n')

            # Cell connectivity
            write('<DataArray type="Int32" Name="connectivity" format="ascii">\n')
            if suffix == '_vertices':
                write(_format_rows(cells[0],'%d\n'))
                write('\n')
            else:
                for connectivity in cells: # rows separated by newlines, no newline after the last row
                    write(_format_rows(connectivity,' '.join(['%d']*connectivity.shape[1])+'\n')[:-1])
                    write('\n')
            write('</DataArray>\n')

            # Cell offsets
            write('<DataArray type="Int32" Name="offsets" format="ascii">\n')
            offsets = np.cumsum(np.concatenate([np.full(len(connectivity),connectivity.shape[1]) for connectivity in cells]))
            write(_format_rows(offsets,'%d\n'))
            write('</DataArray>\n')

            # Cell type
            write('<DataArray type="UInt8" Name="types" format="ascii">\n')
            write(_format_rows(types,'%d\n'))
            write('</DataArray>\n')

            write('</Cells>\n')
            # </Cells>

            # <CellData>
            if cell_data is not None:
                name, values = cell_data
                write('<CellData Tensors="" '+('Vevtors' if suffix == '_conns' else 'Vectors')+'="" Scalars="'+name+'">\n')
                write('<DataArray type="Float32" Name="'+name+'" format="ascii">\n')
                write(_format_rows(values,'%11.8e\n'))
                write('</DataArray>\n')
                write('</CellData>\n')
            # </CellData>

            write('</Piece>\n')
            write('</UnstructuredGrid>\n')
            write('</VTKFile>\n')
    
    
def BezierExtraction(NURBS_degree,nbeam_total):
    nctrlpt_per_beam = 2*NURBS_degree + 1
    knotVec_per_beam = np.concatenate((np.zeros(NURBS_degree),(np.linspace(0,1,int((nctrlpt_per_beam-1)/NURBS_degree+1))),np.ones(NURBS_degree)))
    knotVec = np.tile(knotVec_per_beam, (nbeam_total,1))
    return knotVec
        

def BezierBeamFile(geoName,NURBS_degree,nctrlpt_per_beam,\
                   nconnector_t_per_beam,npatch,knotVec,path='meshes'):
############################# Abaqus txt File #################################
    txtfile = open (Path(path, geoName, geoName + 'IGA.txt'),'w')
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
    
            
def AbaqusFile(geoName,NURBS_degree,npatch,IGAvertices,beam_connectivity,\
               connector_t_bot_connectivity,connector_t_reg_connectivity,\
               connector_t_top_connectivity,connector_l_connectivity,nbeamElem,\
               nconnector_t,nconnector_l,nconnector_t_precrack,nconnector_l_precrack,\
               segment_length,long_connector_ratio,\
               x_max,x_min,y_max,y_min,z_coord,box_shape,box_size,\
               cellwallthickness_sparse,cellwallthickness_dense,\
               merge_operation,merge_tol,\
               precrackFlag='off',precrack_elem=None,path='meshes'):
    
    """Write the Abaqus input file (geoName.inp) with beam (B32) and connector (T3D2) elements"""
    if precrack_elem is None:
        precrack_elem = []
    height_connector_t = segment_length/4 # height of bot/top transverse connectors [mm]
    area_connector_l = cellwallthickness_sparse*cellwallthickness_dense # sectional area of longitudinal connectors [mm^2]
    
    z_min = np.min(z_coord)
    z_max = np.max(z_coord)
    # boundary_conditions = ['Hydrostatic']
    boundary_conditions = ['Bottom','Top','Left','Right','Front','Back']

    if NURBS_degree not in [2,3]:
        raise ValueError('NURBS degree {} is not supported. Check variable NURBS_degree.'.format(NURBS_degree))
    nelnode_connector = 2
    
    numnode = IGAvertices.shape[0]
    nelem = beam_connectivity.shape[0]
    nnode = beam_connectivity.shape[1]
    nelem_connector_t_bot = connector_t_bot_connectivity.shape[0]
    nelem_connector_t_reg = connector_t_reg_connectivity.shape[0]
    nelem_connector_t_top = connector_t_top_connectivity.shape[0]
    nelem_connector_l = connector_l_connectivity.shape[0]
    # precracked connectors are taken out of the bot/reg/top transverse connectors, they are not extra elements
    nelem_total = nelem + nelem_connector_t_bot + nelem_connector_t_reg + nelem_connector_t_top + nelem_connector_l
    
    # transverse connectors of precracked ridges (cell walls) are cracked in every segment
    nsegments = len(z_coord)//(2*NURBS_degree+1)
    nridge = nelem_connector_t_bot//nsegments
    precrack_ridges = [int(ridge) for ridge in precrack_elem] if _is_on(precrackFlag) else []
    def precracked_rows(connectivity): # rows of the bot/reg/top transverse connector arrays, which repeat every nridge rows
        return np.isin(np.arange(connectivity.shape[0]) % max(nridge,1),precrack_ridges)
    
    xaxis = 1
    yaxis = 2
    zaxis = 3

    # Generate a .inp file which can be directly imported and played in Abaqus
    meshinpfile = open(Path(path, geoName, geoName + '.inp'),'w')
    meshinpfile.write('*HEADING'+'\n')
    meshinpfile.write('** Job name: {0} Model name: Model-{1}\n'.format(geoName,geoName))
    meshinpfile.write('** Generated with RingsPy Mesh Generation Tool V{}\n'.format(__version__))
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
    meshinpfile.write(_format_rows(np.column_stack((np.arange(1,numnode+1),IGAvertices[:,0:3])),'%d, %#.9e, %#.9e,  %#.9e\n'))
    def write_elements(ids,connectivity):
        """Write one line per element: id, node1, node2, ..."""
        connectivity = np.reshape(connectivity,(len(ids),-1))
        meshinpfile.write(_format_rows(np.column_stack((ids,connectivity)).astype(int),'%d'+', %d'*connectivity.shape[1]+'\n'))
    # beam element connectivity 
    count = 1
    meshinpfile.write('*Element, type=B32,elset=AllBeams\n') 
    write_elements(np.arange(count,count+nelem),beam_connectivity[:,0:nnode])
    count += nelem

    # bottom/regular/top transverse connector element connectivity, precracked connectors are collected
    connector_t_precracked_index = []
    connector_t_precracked_connectivity = []
    for elset, connectivity in [('AllBotConns',connector_t_bot_connectivity),\
                                ('AllRegConns',connector_t_reg_connectivity),\
                                ('AllTopConns',connector_t_top_connectivity)]:
        meshinpfile.write('*Element, type=T3D2,elset={:s}\n'.format(elset))
        ids = np.arange(count,count+connectivity.shape[0])
        cracked = precracked_rows(connectivity)
        write_elements(ids[~cracked],connectivity[~cracked,0:nelnode_connector])
        connector_t_precracked_index.extend(ids[cracked].tolist())
        connector_t_precracked_connectivity.extend(connectivity[cracked,0:nelnode_connector])
        count += connectivity.shape[0]
    # precracked transverse connector element connectivity
    if connector_t_precracked_index:
        meshinpfile.write('*Element, type=T3D2,elset=AllPrecrackTConns\n')
        write_elements(connector_t_precracked_index,connector_t_precracked_connectivity)
    # longitudinal connector element connectivity (longitudinal connectors are not precracked)
    meshinpfile.write('*Element, type=T3D2,elset=AllLongConns\n')
    write_elements(np.arange(count,count+nelem_connector_l),connector_l_connectivity[:,0:nelnode_connector])
    count += nelem_connector_l

    # Ghost mesh for easy visualization
    count_offset = 10**(int(math.log10(count))+1) # set an offset with the order of magnitude of the max number + 1
    count = count_offset + 1
    
    if NURBS_degree == 2:
        meshinpfile.write('*ELEMENT, TYPE=B32, ELSET=VisualBeams\n')
        write_elements(np.arange(count,count+nelem),beam_connectivity[:,0:nnode])
        count += nelem

    visual_connector_precracked_index = []
    for elset, connectivity in [('VisualBotConns',connector_t_bot_connectivity),\
                                ('VisualRegConns',connector_t_reg_connectivity),\
                                ('VisualTopConns',connector_t_top_connectivity)]:
        meshinpfile.write('*ELEMENT, TYPE=T3D2, ELSET={:s}\n'.format(elset))
        ids = np.arange(count,count+connectivity.shape[0])
        cracked = precracked_rows(connectivity)
        write_elements(ids[~cracked],connectivity[~cracked,0:nelnode_connector])
        visual_connector_precracked_index.extend(ids[cracked].tolist())
        count += connectivity.shape[0]
    # precracked transverse connector element connectivity
    if visual_connector_precracked_index:
        meshinpfile.write('*ELEMENT, TYPE=T3D2, ELSET=VisualPrecrackTConns\n')
        write_elements(visual_connector_precracked_index,connector_t_precracked_connectivity)
    # longitudinal connector element connectivity
    meshinpfile.write('*ELEMENT, TYPE=T3D2, ELSET=VisualLongConns\n')
    write_elements(np.arange(count,count+nelem_connector_l),connector_l_connectivity[:,0:nelnode_connector])
    count += nelem_connector_l
            
    area_connector_t = height_connector_t*cellwallthickness_dense # sectional area of transverse connectors [mm^2]
    meshinpfile.write('** Section: Section-1\n')
    meshinpfile.write('*Solid Section, elset=VisualBotConns, material=VisualTConns\n')
    meshinpfile.write('{:8.4E}\n'.format(area_connector_t))
    meshinpfile.write('** Section: Section-2\n')
    meshinpfile.write('*Solid Section, elset=VisualRegConns, material=VisualTConns\n')
    meshinpfile.write('{:8.4E}\n'.format(area_connector_t))
    meshinpfile.write('** Section: Section-3\n')
    meshinpfile.write('*Solid Section, elset=VisualTopConns, material=VisualTConns\n')
    meshinpfile.write('{:8.4E}\n'.format(area_connector_t))
    meshinpfile.write('** Section: Section-11\n')
    meshinpfile.write('*Solid Section, elset=AllBotConns, material=VisualTConns\n')
    meshinpfile.write('{:8.4E}\n'.format(area_connector_t))
    meshinpfile.write('** Section: Section-12\n')
    meshinpfile.write('*Solid Section, elset=AllRegConns, material=VisualTConns\n')
    meshinpfile.write('{:8.4E}\n'.format(area_connector_t))
    meshinpfile.write('** Section: Section-13\n')
    meshinpfile.write('*Solid Section, elset=AllTopConns, material=VisualTConns\n')
    meshinpfile.write('{:8.4E}\n'.format(area_connector_t))
    if _is_on(precrackFlag) and connector_t_precracked_index:
        meshinpfile.write('** Section: Section-6\n')
        meshinpfile.write('*Solid Section, elset=VisualPrecrackTConns, material=VisualTConns\n')
        meshinpfile.write('{:8.4E}\n'.format(area_connector_t))
        meshinpfile.write('** Section: Section-16\n')
        meshinpfile.write('*Solid Section, elset=AllPrecrackTConns, material=VisualTConns\n')
        meshinpfile.write('{:8.4E}\n'.format(area_connector_t))
    meshinpfile.write('** Section: Section-4\n')
    meshinpfile.write('*Solid Section, elset=VisualLongConns, material=VisualLConns\n')
    meshinpfile.write('{:8.4E}\n'.format(area_connector_l))
    meshinpfile.write('** Section: Section-14\n')
    meshinpfile.write('*Solid Section, elset=AllLongConns, material=VisualLConns\n')
    meshinpfile.write('{:8.4E}\n'.format(area_connector_l))
    if NURBS_degree == 2:
        meshinpfile.write('** Section: Section-5  Profile: Profile-1\n')
        meshinpfile.write('*Beam Section, elset=VisualBeams, material=VisualBeams, temperature=GRADIENTS, section=RECT\n')
        meshinpfile.write('{:8.4E}, {:8.4E}\n'.format(cellwallthickness_dense,cellwallthickness_dense))
        meshinpfile.write('1.,0.,0.\n')
        meshinpfile.write('** Section: Section-15  Profile: Profile-1\n')
        meshinpfile.write('*Beam Section, elset=AllBeams, material=VisualBeams, temperature=GRADIENTS, section=RECT\n')
        meshinpfile.write('{:8.4E}, {:8.4E}\n'.format(cellwallthickness_dense,cellwallthickness_dense))
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
        if _is_on(merge_operation):
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


def ReadSavedSites(sites_path, radial_growth_rule):
    """Search in the current directory and all directories above it 
    for a file of a particular name.

    Arguments:
    ---------
    sites_path :: pathlib obj, the path of saved file.
    radial_growth_rule :: string, the filename to look for.

    Returns
    -------
    sites and radii arrays loaded from the first directory containing both
    the sites file and the matching radii file (name with '_sites' -> '_radii')

    Raises
    ------
    FileNotFoundError, if no directory up to the root contains the files
    """

    sitesfile = radial_growth_rule
    radiifile = radial_growth_rule.strip().replace("_sites", '_radii')
    sites_path = Path(sites_path).resolve()

    while True:
        sites_attempt = sites_path / sitesfile
        radii_attempt = sites_path / radiifile
        if sites_attempt.exists() and radii_attempt.exists():
            print('Sites info from file: {:s} has been loaded.'.format(str(sites_attempt)))
            print('Radii info from file: {:s} has been loaded.'.format(str(radii_attempt)))
            return np.load(sites_attempt), np.load(radii_attempt)
        if sites_path == sites_path.parent: # reached the root (e.g. / or C:\)
            break
        sites_path = sites_path.parent

    raise FileNotFoundError('Could not find files: {:s} and {:s}, please check if the existing site files are under the same directory with the input script.'.format(sitesfile,radiifile))


def StlModelFile(geoName,path='meshes'):

    """Generate the 3D model file (geoName.stl file, binary) for the generated geometry
    from the connector volume file (geoName_conns_vol.vtu).

    Arguments:
    ---------
    geoName :: string, geometry name.
    path :: folder containing the geoName folder.
    """

    try:
        # only the VTK modules used here ('import vtk' loads all of VTK, which takes seconds)
        from vtkmodules.vtkIOXML import vtkXMLUnstructuredGridReader
        from vtkmodules.vtkFiltersGeometry import vtkDataSetSurfaceFilter
        from vtkmodules.vtkFiltersCore import vtkTriangleFilter
        from vtkmodules.vtkIOGeometry import vtkSTLWriter
    except ImportError:
        try:
            from vtk import vtkXMLUnstructuredGridReader, vtkDataSetSurfaceFilter, vtkTriangleFilter, vtkSTLWriter
        except ImportError as e:
            raise ImportError("STL export (stlFlag = 'on') needs the optional package vtk: pip install 'ringspy[vtk]'") from e

    file_name = str(Path(path, geoName, geoName + '_conns_vol'+'.vtu'))
    stlfilename = str(Path(path, geoName, geoName + '.stl'))

    # create a new 'XML Unstructured Grid Reader'
    reader = vtkXMLUnstructuredGridReader()
    reader.SetFileName(file_name)
    reader.Update()
    output = reader.GetOutputPort()

    # Extract the outer (polygonal) surface.
    surface = vtkDataSetSurfaceFilter()
    surface.SetInputConnection(0, output)
    surface.Update()

    # Split the quads of the surface into triangles. vtkSTLWriter would do the same
    # split itself, but about 100 times slower.
    triangles = vtkTriangleFilter()
    triangles.SetInputConnection(surface.GetOutputPort())
    triangles.Update()

    # Write the stl file to disk
    stlWriter = vtkSTLWriter()
    stlWriter.SetFileTypeToBinary() # ASCII STL files are ~3x larger
    stlWriter.SetFileName(stlfilename)
    stlWriter.SetInputConnection(triangles.GetOutputPort())
    stlWriter.Write()


def ModelInfo(box_shape,boundary_points,z_min,z_max,skeleton_density,MeshData,polygon=None):

    """Calculate the material properties of generated geometry

    Arguments:
    ---------
    polygon :: vertices of the clipping polygon (box_shape 'polygon'): the bulk volume is its area times
        (z_max - z_min), with or without boundaryFlag; for the other shapes it is the convex hull of the
        boundary points (without the box corners when boundaryFlag is off) times (z_max - z_min)

    Returns
    -------
    total mass, bulk volume, bulk density, and porosity of the model
    """

    mass = 0
    
    for i in range(0,MeshData.shape[0]):
        mass += skeleton_density*MeshData[i,11]*MeshData[i,12]*np.linalg.norm(MeshData[i,5:8] - MeshData[i,8:11])

    
    if polygon is not None: # area of the clipping polygon (the convex hull is too big if it is not convex)
        polygon = np.asarray(polygon,dtype=float)
        area = 0.5*abs(np.sum(cross2D(polygon,np.roll(polygon,-1,axis=0))))
    else:
        area = ConvexHull(boundary_points).volume
    
    volume = (z_max-z_min)*area

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
        startTime,placementTime,voronoiTime,RebuildvorTime,BeamTime,FileTime,path='meshes',seed=None,box_polygon=None):

    """Generate the log file (geoName.log file) for the generation procedure.

    Arguments:
    ---------
    geoName :: string, geometry name.
    -
    seed :: int seed of the cell placement, or None for the global NumPy random state.

    Returns
    -------
    None or error message, if generation failed (TBD)
    """
    
    # get current local time
    current_time = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    
    # Generate log file
    prefix = str(Path(path, geoName, geoName))
    logfile = open(Path(path, geoName, geoName + '-report.log'),'w')
    logfile.write('################################################################################\n')
    logfile.write('##             Log Created with RingsPy Mesh Generation Tool V' + __version__ + '           ##\n')
    logfile.write('##             Files generated at local time: ' + current_time + '             ##\n')
    logfile.write('################################################################################\n')
    logfile.write('\n')
    logfile.write('GEOMETRY:\n') 
    logfile.write('geoName:                                 ' + geoName + '\n')  
    logfile.write('\n')
    logfile.write('GENERATION PARAMETERS:\n')
    logfile.write('iter_max:                                ' + str(iter_max) + '\n')
    logfile.write('radial_growth_rule:                      ' + str(radial_growth_rule) + '\n')
    logfile.write('seed:                                    ' + str(seed) + '\n')
    logfile.write('merge_operation:                         ' + _flag_text(merge_operation) + '\n')
    logfile.write('boundaryFlag:                            ' + _flag_text(boundaryFlag) + '\n')
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

    if _is_on(merge_operation): 
        logfile.write('merge_tol:                               ' + str(merge_tol) + '\n')

    logfile.write('\n')
    logfile.write('CLIPPING PARAMETERS:\n')
    logfile.write('box_shape:                               ' + str(box_shape) + '\n')
    logfile.write('box_center:                              ' + str(box_center) + '\n')
    logfile.write('box_size:                                ' + str(box_size) + '\n')
    if box_polygon is not None:
        logfile.write('box_polygon:                             ' + str([[float(x), float(y)] for x, y in box_polygon]) + '\n')
    logfile.write('precrackFlag:                            ' + _flag_text(precrackFlag) + '\n')
    if _is_on(precrackFlag):
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
    logfile.write('Cross-sectional Image File:              ' + prefix \
        + '.png\n')
    logfile.write('NURBS Beam File:                         ' + prefix \
        + 'IGA.txt\n')
    logfile.write('Connector Data File:                     ' + prefix \
        + '-mesh.txt\n')
    logfile.write('Grain-ridge Data File:                   ' + prefix \
        + '-vertex.mesh\n')
    logfile.write('Ridge Data File:                         ' + prefix \
        + '-ridge.mesh\n')
    logfile.write('Paraview Vertices File:                  ' + prefix
                  
        + '_vertices'+'.vtu\n')
    logfile.write('Paraview Beams File:                     ' + prefix \
        + '_beams'+'.vtu\n')
    logfile.write('Paraview Connectors File:                ' + prefix \
        + '_conns'+'.vtu\n')
    logfile.write('Paraview Connectors (Volume) File:       ' + prefix \
        + '_conns_vol'+'.vtu\n')
    if _is_on(stlFlag):
        logfile.write('3D Model File:                           ' + prefix \
            + '.stl\n')
    if _is_on(inpFlag):
        if inpType in ABAQUS_TYPES:
            logfile.write('Abaqus Input File:                       ' + prefix \
            + '.inp\n') 
    logfile.write('Log File:                                ' + prefix \
        + '-report.log\n')
    logfile.close()