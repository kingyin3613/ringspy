# -*- coding: utf-8 -*-
"""
Synthetic 3D images (voxel volumes) of RingsPy cellular geometries

The cells of a 2D RingsPy tessellation are extruded along a fiber direction field (a global
tilt of the fibers, a twist around the generation center and a smooth random wobble, so that
the fibers grow roughly in the same longitudinal direction but are not aligned neatly) and
rasterised slice by slice into a volume of voxels, as a CT scan of the specimen would give.
The volume is meant to test and demonstrate the reconstruction of a lattice from a 3D image
(ringspy.reconstruct) against a known geometry.

author: king_yin3613
email: haoyin2022@u.northwestern.edu
"""

import math
from dataclasses import dataclass, field
from typing import Callable, Optional, Tuple

import numpy as np
from scipy.spatial import Voronoi, cKDTree

from . import MeshGenTools as rpgen


@dataclass
class SyntheticVolume:
    """A synthetic voxel image and the geometry it was made from"""
    volume: np.ndarray                 # (nx, ny, nz) array; bool (True = cell wall) or float32 grayscale
    voxel_size: Tuple[float, float, float]
    origin: Tuple[float, float, float]  # (x, y, z) of the corner of the first voxel
    sites: np.ndarray                  # cell sites at z_min, one row (x, y, cell size) per cell
    wall_thickness: np.ndarray         # wall thickness of each cell (mm)
    sites_at: Callable[[float], np.ndarray] = field(repr=False)  # sites_at(z): the (x, y) sites at height z
    tilt: Tuple[float, float] = (0.0, 0.0)
    twist: float = 0.0
    wobble_amplitude: float = 0.0
    cell_ends: Optional[list] = None   # heights (mm) at which each cell ends and the next one begins (fiber_length)

    @property
    def shape(self):
        return self.volume.shape

    @property
    def box(self):
        """(x_min, x_max, y_min, y_max, z_min, z_max) of the volume"""
        nx, ny, nz = self.volume.shape
        x0, y0, z0 = self.origin
        vx, vy, vz = self.voxel_size
        return (x0, x0 + nx*vx, y0, y0 + ny*vy, z0, z0 + nz*vz)

    def vertices_at(self, z, margin=0.0):
        """The Voronoi vertices (cell corners) of the cells at height z that lie inside the
        volume, at least `margin` away from its sides (the tapering of the cells near their ends is
        not accounted for)"""
        x_min, x_max, y_min, y_max, _, _ = self.box
        sites = self.sites_at(z)
        vertices = Voronoi(sites).vertices
        inside = (vertices[:, 0] > x_min + margin) & (vertices[:, 0] < x_max - margin) & \
                 (vertices[:, 1] > y_min + margin) & (vertices[:, 1] < y_max - margin)
        return vertices[inside]

    @property
    def porosity(self):
        """Void fraction of the volume (of the thresholded image for a grayscale volume)"""
        solid = self.volume if self.volume.dtype == bool else self.volume > 0.5
        return 1.0 - solid.mean()


def _wobble_field(amplitude, wavelength, rng, nwaves=6):
    """A smooth random in-plane displacement field w(x, y, z) -> (dx, dy): a sum of random plane waves
    with wave numbers around 2*pi/wavelength and a standard deviation of about `amplitude`"""
    if amplitude <= 0:
        return lambda points, z: np.zeros((len(points), 2))
    directions = rng.normal(size=(2, nwaves, 3))
    directions /= np.linalg.norm(directions, axis=2, keepdims=True)
    wavenumbers = 2*math.pi/wavelength*rng.uniform(0.7, 1.3, size=(2, nwaves, 1))
    k = directions*wavenumbers                                     # (2, nwaves, 3)
    phases = rng.uniform(0, 2*math.pi, size=(2, nwaves))
    scale = amplitude*math.sqrt(2.0/nwaves)                        # std of the sum ~ amplitude

    def field_at(points, z):
        xyz = np.column_stack((points[:, 0], points[:, 1], np.full(len(points), z)))
        w = np.empty((len(points), 2))
        for i in range(2):
            w[:, i] = scale*np.sin(xyz @ k[i].T + phases[i]).sum(axis=1)
        return w
    return field_at


def _resolve_geometry(geometry, rng):
    """Sites, wall thicknesses, box and generation center from a Config (or TOML file/dict) or
    from a dict of arrays {'sites', 'wall_thickness', 'box', 'generation_center'}"""
    from .config import Config
    from .pipeline import place_cells
    if isinstance(geometry, dict) and 'sites' in geometry:
        sites = np.asarray(geometry['sites'], dtype=float)
        if sites.ndim != 2 or sites.shape[1] < 2:
            raise ValueError("sites must be an array of (x, y[, size]) rows")
        if sites.shape[1] == 2:
            sites = np.hstack((sites, np.zeros((len(sites), 1))))
        thickness = np.broadcast_to(np.asarray(geometry['wall_thickness'], dtype=float), (len(sites),)).copy()
        box = tuple(float(v) for v in geometry['box'])
        if len(box) != 6:
            raise ValueError("box must be (x_min, x_max, y_min, y_max, z_min, z_max)")
        center = tuple(geometry.get('generation_center', (0.0, 0.0)))
        return sites, thickness, box, center
    if isinstance(geometry, dict):
        geometry = Config.from_dict(geometry)
    elif not isinstance(geometry, Config):
        geometry = Config.from_toml(geometry)
    c = geometry.resolved()
    sites, radii = place_cells(c, rng)
    x_min, x_max, y_min, y_max, *_ = rpgen.Clipping_Box(c.box_shape, c.box_center, c.box_size, c.boundaryFlag,
                                                        polygon=c.box_polygon)
    # the wall thickness rule of VertexandRidgeinfo, applied to the cells: the thickness of the growth rings
    distance = np.linalg.norm(sites[:, 0:2] - np.asarray(c.generation_center, dtype=float), axis=1)
    profile = rpgen.RingProfile(radii, c.cellwallthickness_sparse, c.cellwallthickness_dense, c.transition,
                                c.transition_width)
    thickness = np.array(profile(distance), dtype=float).reshape(-1)
    return sites, thickness, (x_min, x_max, y_min, y_max, c.z_min, c.z_max), tuple(c.generation_center)


def voxelize(geometry, voxel_size, tilt=(0.0, 0.0), twist=0.0, wobble_amplitude=0.0, wobble_wavelength=None,
             seed=None, noise=0.0, blur=0.0, margin=0.0, progress=False, fiber_length=None, fiber_length_cov=0.3,
             fiber_end_offset=None, taper_length=None):
    """
    Rasterise a RingsPy cell geometry into a synthetic 3D image, slice by slice.

    geometry :: a ringspy.Config (or TOML file name or dict of parameters): its cells are placed with its
        growth rule and `seed`, the volume covers the bounding box of its clipping box, from z_min to
        z_max, and the cell walls are as thick as cellwallthickness_sparse/dense. Or a dict with the
        arrays 'sites' (rows (x, y[, size])), 'wall_thickness' (per cell, or one number), the 'box'
        (x_min, x_max, y_min, y_max, z_min, z_max) and optionally 'generation_center'.
    voxel_size :: edge length of the (cubic) voxels in mm, or (vx, vy, vz).
    tilt :: (tx, ty), the slopes dx/dz and dy/dz of the fibers (the tangent of the grain angles): all
        fibers lean this way.
    twist :: twist angle (radian) of the cross-section around the generation center from z_min to z_max,
        as theta_max - theta_min of the mesh generation.
    wobble_amplitude, wobble_wavelength :: a smooth random in-plane displacement field (mm) of the given
        standard deviation and wavelength (default: a third of the box width) bends the fibers, so they
        are not straight and not parallel. The fiber slopes it adds are about 2*pi*amplitude/wavelength.
    seed :: seed of the cell placement and of the wobble.
    noise, blur :: when either is > 0 the volume is a float32 grayscale image, 1 in the walls and 0 in the
        lumens, blurred with a Gaussian of `blur` voxels and with Gaussian noise of standard deviation
        `noise` added, like a scan. Otherwise it is boolean (True = cell wall).
    margin :: extra width (mm) of cells kept around the box when placing sites (the Voronoi cells of the
        sites near the box sides are correct only with their neighbours present; the growth rules place
        cells well beyond the box, so 0 is fine for them).
    fiber_length, fiber_length_cov :: the cells end and new ones begin at random heights: cell lengths from a
        lognormal distribution of this mean (mm) and coefficient of variation, the first cell of each site
        crossing the bottom face at a random point of its length. A cell ends with a transverse wall as thick
        as its walls, and the next cell is centered up to `fiber_end_offset` (mm, default a quarter of the
        cell size) away. Default: no cell ends inside the volume.
    taper_length :: the cells taper to their ends (and widen from their beginnings) over this length (mm,
        default 0.15*fiber_length; 0: blunt ends): the walls of the cell move in linearly, so its lumen
        closes like a cone, and the neighbouring cells fill the space.

    Returns a SyntheticVolume: the (nx, ny, nz) array with the longitudinal (fiber) direction along the
    last axis, voxel (i, j, k) covering x in [x_min + i*vx, x_min + (i+1)*vx] etc., and the geometry.
    """
    rng = np.random.default_rng(seed)
    sites, thickness, box, center = _resolve_geometry(geometry, rng)
    x_min, x_max, y_min, y_max, z_min, z_max = box
    voxel = np.broadcast_to(np.asarray(voxel_size, dtype=float), (3,)).copy()
    if np.any(voxel <= 0):
        raise ValueError("voxel_size must be positive")
    nx, ny, nz = [max(int(round((hi - lo)/v)), 1) for lo, hi, v in
                  [(x_min, x_max, voxel[0]), (y_min, y_max, voxel[1]), (z_min, z_max, voxel[2])]]
    height = max(z_max - z_min, voxel[2])
    tilt = np.asarray(tilt, dtype=float).reshape(2)
    center = np.asarray(center, dtype=float)
    if wobble_wavelength is None:
        wobble_wavelength = max(x_max - x_min, y_max - y_min)/3.0
    wobble = _wobble_field(wobble_amplitude, wobble_wavelength, rng)

    sites_xy = sites[:, 0:2]
    wobble_ref = wobble(sites_xy, z_min)

    # cells ending inside the volume: end heights, and the offset of the next cell
    cell_ends = None
    end_site = np.zeros(0, dtype=int)
    end_z = np.zeros(0)
    end_jump = np.zeros((0, 2))
    if fiber_length is not None:
        if fiber_length <= 0 or fiber_length_cov < 0:
            raise ValueError("fiber_length must be positive and fiber_length_cov not negative")
        if fiber_end_offset is None:
            fiber_end_offset = 0.5*float(np.median(sites[:, 2])) if np.any(sites[:, 2] > 0) else 0.0
        sigma = math.sqrt(math.log(1 + fiber_length_cov**2))
        mu = math.log(fiber_length) - sigma**2/2
        cell_ends = []
        for i in range(len(sites)):
            z = z_min + rng.lognormal(mu, sigma)*(1 - rng.random())
            ends = []
            while z < z_max:
                ends.append(z)
                z += rng.lognormal(mu, sigma)
            cell_ends.append(np.array(ends))
        end_site = np.concatenate([np.full(len(e), i) for i, e in enumerate(cell_ends)]).astype(int)
        end_z = np.concatenate(cell_ends)
        radius = fiber_end_offset*np.sqrt(rng.random(len(end_z)))
        angle = rng.uniform(0, 2*math.pi, len(end_z))
        end_jump = np.column_stack((radius*np.cos(angle), radius*np.sin(angle)))
        if taper_length is None:
            taper_length = 0.15*fiber_length
    # the cell closes when the walls have moved in by the half distance to its neighbours
    half_spacing = cKDTree(sites_xy).query(sites_xy, k=min(2, len(sites)))[0][:, -1]/2

    def sites_at(z):
        """(x, y) of the cell sites at height z: twisted, tilted and wobbled, and moved at every cell end below z"""
        angle = twist*(z - z_min)/height
        rotated = center + (sites_xy - center) @ np.array([[math.cos(angle), -math.sin(angle)],
                                                            [math.sin(angle), math.cos(angle)]]).T
        moved = rotated + (z - z_min)*tilt + wobble(sites_xy, z) - wobble_ref
        below = end_z < z
        if below.any():
            np.add.at(moved, end_site[below], end_jump[below])
        return moved

    # pixel centers of a slice
    xs = x_min + (np.arange(nx) + 0.5)*voxel[0]
    ys = y_min + (np.arange(ny) + 0.5)*voxel[1]
    pixels = np.column_stack([a.reshape(-1) for a in np.meshgrid(xs, ys, indexing='ij')])
    k_near = min(4, len(sites))
    if k_near < 2:
        raise ValueError("at least 2 cells are needed")

    volume = np.zeros((nx, ny, nz), dtype=bool)
    for k in range(nz):
        z = z_min + (k + 0.5)*voxel[2]
        moved = sites_at(z)
        distance, index = cKDTree(moved).query(pixels, k=k_near)
        power = distance**2
        if len(end_z) and taper_length > 0:
            # tapering cells: a negative power weight moves their walls in, up to the site at the very end
            closing = np.zeros(len(sites))
            nearness = np.maximum(0.0, 1 - np.abs(end_z - z)/taper_length)
            np.maximum.at(closing, end_site, nearness)
            power = power + (2*half_spacing[index]*closing[index])**2
            order = np.argsort(power, axis=1, kind='stable')
            power = np.take_along_axis(power, order, axis=1)
            index = np.take_along_axis(index, order, axis=1)
        solid = np.zeros(len(pixels), dtype=bool)
        for j in range(1, k_near):
            # distance from the pixel to the bisector (cell wall) between its nearest cell and the j-th nearest
            separation = np.linalg.norm(moved[index[:, j]] - moved[index[:, 0]], axis=1)
            to_wall = (power[:, j] - power[:, 0])/(2*np.maximum(separation, 1e-300))
            half_wall = 0.5*(thickness[index[:, 0]] + thickness[index[:, j]])/2
            solid |= to_wall < half_wall
        if len(end_z):  # the transverse walls closing the cells that end at this height
            closing = np.zeros(len(sites), dtype=bool)
            closing[end_site[np.abs(end_z - z) <= thickness[end_site]/2]] = True
            solid |= closing[index[:, 0]]
        volume[:, :, k] = solid.reshape(nx, ny)
        if progress and (k % 50 == 0 or k == nz - 1):
            print('slice {:d}/{:d} rasterised'.format(k + 1, nz))

    if noise > 0 or blur > 0:
        from scipy.ndimage import gaussian_filter
        image = volume.astype(np.float32)
        if blur > 0:
            image = gaussian_filter(image, blur)
        if noise > 0:
            image += rng.normal(0, noise, size=image.shape).astype(np.float32)
        volume = image

    return SyntheticVolume(volume=volume, voxel_size=tuple(voxel), origin=(x_min, y_min, z_min), sites=sites,
                           wall_thickness=thickness, sites_at=sites_at, tilt=tuple(tilt), twist=twist,
                           wobble_amplitude=wobble_amplitude, cell_ends=cell_ends)
