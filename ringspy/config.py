# -*- coding: utf-8 -*-
"""
Input parameters of a RingsPy mesh generation, loadable from a TOML file

author: king_yin3613
email: haoyin2022@u.northwestern.edu
"""

import dataclasses
import difflib
import os
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Optional, Tuple, Union, get_args, get_origin

from .MeshGenTools import (ABAQUS_TYPES, HEXAGON_SHAPES, POLYGON_SHAPES, SQUARE_SHAPES, TRIANGLE_SHAPES,
                           polygon_vertices)

if sys.version_info >= (3, 11):
    import tomllib
else:
    import tomli as tomllib

GROWTH_RULES = ['binary', 'binary_lloyd', 'regular_hexagonal']
BOX_SHAPES = SQUARE_SHAPES + TRIANGLE_SHAPES + HEXAGON_SHAPES + POLYGON_SHAPES
TRUE_STRINGS = ['on', 'y', 'yes', 'true']
FALSE_STRINGS = ['off', 'n', 'no', 'false']


@dataclass
class Config:
    """Input parameters of a mesh generation.

    The names are those of the template scripts and of the report log; see the README for
    their meaning. The defaults generate the ``wood_cube`` example. Lengths are in mm.
    Flags take booleans; the strings 'on'/'off' (and yes/no, y/n, true/false) also work.
    """

    geoName: str = 'wood_cube'
    path: str = 'meshes'                     # output folder, the files go to path/geoName/

    # cell placement
    radial_growth_rule: str = 'binary'       # binary, binary_lloyd, regular_hexagonal, or a *_sites.npy file
    iter_max: int = 500
    print_interval: int = 500
    seed: Optional[int] = None               # None: global NumPy random state (np.random.seed), as before v1.2
    sites_dir: Optional[str] = None          # where to look for a *_sites.npy file (and its parents);
                                             # default: the folder of the TOML file, else the current folder

    # radial cell growth
    r_min: float = 0
    r_max: float = 4
    nrings: int = 4
    width_heart: Optional[float] = None      # default 0.3*(r_max-r_min)/nrings
    width_sparse: Optional[float] = None     # default 0.7*(r_max-r_min)/nrings
    width_dense: Optional[float] = None      # default 0.3*(r_max-r_min)/nrings
    generation_center: Tuple[float, float] = (0, 0)
    cellsize_sparse: float = 0.02
    cellsize_dense: float = 0.01
    cellwallthickness_sparse: float = 0.006
    cellwallthickness_dense: float = 0.010

    # clipping box
    boundaryFlag: bool = True
    box_shape: str = 'square'                # square, triangle, hexagon or polygon
    box_center: Tuple[float, float] = (0, 0)
    box_size: float = 4.0
    box_polygon: Optional[List[Tuple[float, float]]] = None  # vertices [[x, y], ...] of the polygon for
                                             # box_shape 'polygon' (any simple polygon, also non-convex)

    # longitudinal direction
    nsegments: int = 2
    theta_min: float = 0                     # radian
    theta_max: float = 0.05                  # radian
    z_min: float = 0
    z_max: Optional[float] = None            # default box_size (polygon: z_min + its largest width)
    long_connector_ratio: float = 0.02

    # material
    skeleton_density: float = 1.5e-9         # tonne/mm3

    # mesh operations
    merge_operation: bool = True
    merge_tol: float = 0.01
    precrackFlag: bool = False
    precrack_widths: float = 0.1
    precrack_nodes: Optional[List[Tuple[float, float, float, float]]] = None  # [[x1, y1, x2, y2], ...];
                                             # default: one horizontal crack from the left edge of the box

    # output files
    stlFlag: bool = True
    inpFlag: bool = True
    inpType: str = 'Abaqus'
    save_sites: bool = False                 # also save geoName_sites.npy and geoName_radii.npy

    def __post_init__(self):
        for f in dataclasses.fields(self):
            setattr(self, f.name, _check(f.name, getattr(self, f.name), f.type))
        if not (self.radial_growth_rule in GROWTH_RULES or self.radial_growth_rule.endswith('.npy')):
            raise ValueError("radial_growth_rule must be one of {} or a .npy file name, not {!r}"
                             .format(', '.join(GROWTH_RULES), self.radial_growth_rule))
        if self.box_shape not in BOX_SHAPES:
            raise ValueError("box_shape must be square, triangle, hexagon or polygon, not {!r}".format(self.box_shape))
        if self.box_shape in POLYGON_SHAPES:
            if self.box_polygon is None:
                raise ValueError("box_shape 'polygon' needs box_polygon, the vertices [[x, y], ...] of the polygon")
            try:
                polygon_vertices(self.box_polygon)
            except ValueError as e:
                raise ValueError("box_polygon: {}".format(e)) from None
            if self.precrackFlag and self.precrack_nodes is None:
                raise ValueError("with box_shape 'polygon', give the precracks as precrack_nodes")
        elif self.box_polygon is not None:
            raise ValueError("box_polygon is only used with box_shape = 'polygon', not {!r}".format(self.box_shape))
        if self.inpFlag and self.inpType not in ABAQUS_TYPES:
            raise ValueError("inpType {!r} is not supported, use 'Abaqus' (or set inpFlag = false)".format(self.inpType))
        if self.seed is not None and self.seed < 0:
            raise ValueError("seed must be a non-negative integer or None")
        for name in ['iter_max', 'print_interval', 'nrings', 'nsegments']:
            if getattr(self, name) < 1:
                raise ValueError("{} must be at least 1".format(name))
        for name in ['r_max', 'box_size', 'cellsize_sparse', 'cellsize_dense']:
            if getattr(self, name) <= 0:
                raise ValueError("{} must be positive".format(name))

    def resolved(self):
        """Copy with the parameters left at None (ring widths, z_max) replaced by their default values"""
        fractions = {'width_heart': 0.3, 'width_sparse': 0.7, 'width_dense': 0.3}
        values = {name: frac*(self.r_max - self.r_min)/self.nrings  # same expression as the template scripts
                  for name, frac in fractions.items() if getattr(self, name) is None}
        if self.z_max is None:
            if self.box_shape in POLYGON_SHAPES:
                extent = [max(v[i] for v in self.box_polygon) - min(v[i] for v in self.box_polygon) for i in (0, 1)]
                values['z_max'] = self.z_min + max(extent)
            else:
                values['z_max'] = self.box_size
        return dataclasses.replace(self, **values)

    def replace(self, **changes):
        """Copy with some parameters changed, e.g. ``config.replace(seed=1, stlFlag=False)``"""
        _check_names(changes)
        return dataclasses.replace(self, **changes)

    def to_dict(self):
        """Parameters as a dict (tuples become lists, as in TOML)"""
        return {k: _plain(v) for k, v in dataclasses.asdict(self).items()}

    @classmethod
    def from_dict(cls, params):
        """Config from a dict, with a helpful error for unknown names"""
        _check_names(params)
        return cls(**params)

    @classmethod
    def from_toml(cls, file):
        """Config from a TOML file; a relative sites_dir is relative to the file's folder"""
        file = Path(file)
        with open(file, 'rb') as f:
            params = tomllib.load(f)
        for key, value in params.items():
            if isinstance(value, dict):
                raise ValueError("{}: [{}] tables are not supported, write the parameters at the top level"
                                 .format(file, key))
        sites_dir = params.get('sites_dir', '.')
        if isinstance(sites_dir, str):
            params['sites_dir'] = str(file.parent / sites_dir)  # unchanged if sites_dir is absolute
        try:
            return cls.from_dict(params)
        except (TypeError, ValueError) as e:
            raise ValueError("{}: {}".format(file, e)) from None


def load_config(file):
    """Load a Config from a TOML file"""
    return Config.from_toml(file)


def _check_names(params, cls=None):
    """Raise a ValueError for a parameter name that Config (or cls) does not have, suggesting the closest name"""
    names = [f.name for f in dataclasses.fields(Config if cls is None else cls)]
    for key in params:
        if key not in names:
            close = difflib.get_close_matches(key, names, n=1)
            hint = " (did you mean {!r}?)".format(close[0]) if close else ""
            raise ValueError("unknown parameter {!r}{}".format(key, hint))


def _plain(value):
    if isinstance(value, (tuple, list)):
        return [_plain(v) for v in value]
    return value


def _is_number(value):
    return isinstance(value, (int, float)) and not isinstance(value, bool)


def _check(name, value, tp):
    """Check value against the type annotation tp, and convert it (flags to bool, sequences to tuples)"""
    origin, args = get_origin(tp), get_args(tp)
    if origin is Union:  # Optional[X]
        if value is None:
            return None
        return _check(name, value, args[0])
    if tp is bool:
        return to_bool(name, value)
    if tp is int:
        if isinstance(value, int) and not isinstance(value, bool):
            return value
        if hasattr(value, '__index__') and not isinstance(value, bool):  # numpy integers
            return int(value)
        raise TypeError("{} must be an integer, not {!r}".format(name, value))
    if tp is float:
        if _is_number(value):
            return value
        if hasattr(value, 'dtype') and getattr(value, 'shape', None) == ():  # numpy scalars
            return value.item()
        raise TypeError("{} must be a number, not {!r}".format(name, value))
    if tp is str:
        if isinstance(value, str):
            return value
        if isinstance(value, os.PathLike):
            return os.fspath(value)
        raise TypeError("{} must be a string, not {!r}".format(name, value))
    if origin is tuple:
        if isinstance(value, (str, bytes)) or not hasattr(value, '__len__') or len(value) != len(args):
            raise TypeError("{} must be a list of {} numbers, not {!r}".format(name, len(args), value))
        return tuple(_check(name, v, a) for v, a in zip(value, args))
    if origin is list:
        if isinstance(value, (str, bytes)) or not hasattr(value, '__len__'):
            raise TypeError("{} must be a list, not {!r}".format(name, value))
        return [_check(name, v, args[0]) for v in value]
    raise TypeError("unsupported type {} for {}".format(tp, name))


def to_bool(name, value):
    """On/off flag as a bool: accepts booleans and the strings on/off, yes/no, y/n, true/false"""
    if isinstance(value, bool):
        return value
    if hasattr(value, 'dtype') and value.dtype == bool:  # numpy bool
        return bool(value)
    if isinstance(value, str):
        if value.strip().lower() in TRUE_STRINGS:
            return True
        if value.strip().lower() in FALSE_STRINGS:
            return False
    raise ValueError("{} must be true/false (or 'on'/'off'), not {!r}".format(name, value))


# =============================================================================
# Reconstruction from a 3D image
# =============================================================================

TRACKING_MODES = ['labels', 'direction', 'straight']
SOLID_PHASES = ['bright', 'dark']
VTK_FORMATS = ['ascii', 'binary']
REFERENCE_SLICES = ['middle', 'bottom', 'top']
LONGITUDINAL_AXES = ['auto', '0', '1', '2', 'x', 'y', 'z']


@dataclass
class ReconstructionConfig:
    """Input parameters of a lattice reconstruction from a 3D image (ringspy.reconstruct).

    The image is a 3D array of voxels (a .npy/.npz file, or a .tif stack), either boolean (True in the
    cell walls) or grayscale (thresholded into walls and lumens). Lengths are in mm. Parameters left at
    None are derived from the image (their values are written to the report log).
    """

    geoName: str = 'wood_ct'
    path: str = 'meshes'                     # output folder, the files go to path/geoName/

    # the 3D image
    volume: Optional[str] = None             # .npy/.npz/.tif file of the 3D image (relative to the TOML file's
                                             # folder); or pass the array to reconstruct()
    voxel_size: List[float] = field(default_factory=lambda: [0.001, 0.001, 0.001])  # voxel edges (mm) along the
                                             # array axes; one number for cubic voxels
    longitudinal_axis: str = 'auto'          # array axis along the fibers: auto, 0, 1, 2 (or x, y, z)
    solid: str = 'bright'                    # phase of the cell walls: bright (as in X-ray CT) or dark
    threshold: Optional[float] = None        # gray value separating walls and lumens (None: Otsu's threshold;
                                             # not used for a boolean image)
    smooth_sigma: float = 0.0                # Gaussian smoothing (voxels) before thresholding
    crop: Optional[List[int]] = None         # [x0, x1, y0, y1, z0, z1]: voxel index ranges kept, x and y across
                                             # the fibers, z along them (after longitudinal_axis is applied)
    origin: Tuple[float, float, float] = (0.0, 0.0, 0.0)  # coordinates of the corner of the first kept voxel
    min_cell_area: Optional[float] = None    # lumens smaller than this (mm2) count as wall (None: 9 voxels)
    min_wall_area: Optional[float] = None    # wall specks smaller than this (mm2) count as lumen (None: 9 voxels)

    # the network of cell corners and walls, and how the corners are followed along the fibers
    reference_slice: str = 'middle'          # cross-section giving the 2D mesh: middle, bottom, top, a slice
                                             # index, or a height in mm (with a decimal point)
    tracking: str = 'labels'                 # labels (follow the cells), direction (fiber direction field
                                             # only) or straight (plain extrusion of the reference slice)
    link_stride: int = 1                     # label every n-th cross-section when following the cells
    direction_scale: Optional[float] = None  # averaging scale (mm) of the fiber direction field (None: 2 cell diameters)
    max_jump: Optional[float] = None         # largest move (mm) of a corner between two layers (None: 1.5 cell diameters)
    merge_operation: bool = True             # merge cell walls shorter than merge_tol
    merge_tol: Optional[float] = None        # None: 1.5 voxels
    boundaryFlag: bool = True                # walls along the sides of the cross-section
    wall_thickness: Optional[float] = None   # None: measured in the image, wall by wall and layer by layer

    # longitudinal direction
    nsegments: int = 2
    long_connector_ratio: float = 0.02

    # material
    skeleton_density: float = 1.5e-9         # tonne/mm3

    # output files
    stlFlag: bool = True
    inpFlag: bool = True
    inpType: str = 'Abaqus'
    vtk_format: str = 'ascii'                # ascii or binary .vtu files
    save_tracks: bool = True                 # also save geoName-grains.npz: the followed cell corners of every layer

    def __post_init__(self):
        if _is_number(self.voxel_size):
            self.voxel_size = [self.voxel_size]*3
        if isinstance(self.longitudinal_axis, int) and not isinstance(self.longitudinal_axis, bool):
            self.longitudinal_axis = str(self.longitudinal_axis)
        if _is_number(self.reference_slice):
            self.reference_slice = repr(self.reference_slice)
        for f in dataclasses.fields(self):
            setattr(self, f.name, _check(f.name, getattr(self, f.name), f.type))
        if len(self.voxel_size) != 3 or min(self.voxel_size) <= 0:
            raise ValueError("voxel_size must be a positive number or a list of 3 positive numbers, not {!r}"
                             .format(self.voxel_size))
        if self.longitudinal_axis.strip().lower() not in LONGITUDINAL_AXES:
            raise ValueError("longitudinal_axis must be 'auto', 0, 1, 2 or x, y, z, not {!r}".format(self.longitudinal_axis))
        if self.solid not in SOLID_PHASES:
            raise ValueError("solid must be 'bright' or 'dark', not {!r}".format(self.solid))
        if self.tracking not in TRACKING_MODES:
            raise ValueError("tracking must be one of {}, not {!r}".format(', '.join(TRACKING_MODES), self.tracking))
        if self.vtk_format not in VTK_FORMATS:
            raise ValueError("vtk_format must be 'ascii' or 'binary', not {!r}".format(self.vtk_format))
        text = self.reference_slice.strip().lower()
        if text not in REFERENCE_SLICES:
            try:
                float(text)
            except ValueError:
                raise ValueError("reference_slice must be middle, bottom, top, a slice index or a height (mm), not {!r}"
                                 .format(self.reference_slice)) from None
        if self.crop is not None and len(self.crop) != 6:
            raise ValueError("crop must be [x0, x1, y0, y1, z0, z1] (6 voxel indices), not {!r}".format(self.crop))
        if self.inpFlag and self.inpType not in ABAQUS_TYPES:
            raise ValueError("inpType {!r} is not supported, use 'Abaqus' (or set inpFlag = false)".format(self.inpType))
        for name in ['nsegments', 'link_stride']:
            if getattr(self, name) < 1:
                raise ValueError("{} must be at least 1".format(name))
        if self.smooth_sigma < 0:
            raise ValueError("smooth_sigma must not be negative")
        for name in ['min_cell_area', 'min_wall_area', 'direction_scale', 'max_jump', 'merge_tol', 'wall_thickness']:
            value = getattr(self, name)
            if value is not None and value <= 0:
                raise ValueError("{} must be positive (or left out)".format(name))

    def replace(self, **changes):
        """Copy with some parameters changed, e.g. ``config.replace(nsegments=4, stlFlag=False)``"""
        _check_names(changes, type(self))
        return dataclasses.replace(self, **changes)

    def to_dict(self):
        """Parameters as a dict (tuples become lists, as in TOML)"""
        return {k: _plain(v) for k, v in dataclasses.asdict(self).items()}

    @classmethod
    def from_dict(cls, params):
        """ReconstructionConfig from a dict, with a helpful error for unknown names"""
        _check_names(params, cls)
        return cls(**params)

    @classmethod
    def from_toml(cls, file):
        """ReconstructionConfig from a TOML file; a relative `volume` file is relative to the file's folder"""
        file = Path(file)
        with open(file, 'rb') as f:
            params = tomllib.load(f)
        for key, value in params.items():
            if isinstance(value, dict):
                raise ValueError("{}: [{}] tables are not supported, write the parameters at the top level"
                                 .format(file, key))
        volume = params.get('volume')
        if isinstance(volume, str) and not Path(volume).is_absolute():
            params['volume'] = str(file.parent / volume)
        try:
            return cls.from_dict(params)
        except (TypeError, ValueError) as e:
            raise ValueError("{}: {}".format(file, e)) from None
