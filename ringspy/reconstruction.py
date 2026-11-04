# -*- coding: utf-8 -*-
"""
Reconstruction of a RingsPy lattice from a 3D image of a fibrous cellular material

Given a 3D reconstruction of a specimen at the scale of its cells (e.g. a CT scan of a wood cube,
a stack of cross-sectional images), this module finds the cell walls and cell corners in the
cross-sections, follows the corners along the fibers (which grow roughly in the same longitudinal
direction but are not aligned neatly: they lean, wobble and twist) and builds the same lattice as
the mesh generation: grains (beams) along the cell corners, transverse connectors across the cell
walls and longitudinal connectors between the segments, with the measured wall thicknesses.

The steps of ringspy.reconstruct():
  1. load the volume, threshold it (Otsu) into cell walls and lumens, find the longitudinal axis
  2. label the lumens of every cross-section, keeping the label of each cell from slice to slice
  3. in a reference cross-section, build the network of cell corners (junctions of 3 or more
     cells) and cell walls (walls between two cells), as the Voronoi tessellation of the mesh
     generation; each corner is identified by the cells that meet there
  4. the fiber direction field of the volume (structure tensor) and the cell labels follow every
     corner to the heights of all layers of control points of the beams
  5. the extrusion and output functions of MeshGenTools write the usual files

author: king_yin3613
email: haoyin2022@u.northwestern.edu
"""

import datetime
import math
import os
import shutil
import time
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
from scipy import ndimage
from scipy.spatial import cKDTree

from . import MeshGenTools as rpgen
from . import __version__
from . import lattice as rplat
from .synthetic import SyntheticVolume

NURBS_DEGREE = 2

# pseudo labels of the four sides of the (rectangular) cross-section
SIDE_X_MIN, SIDE_X_MAX, SIDE_Y_MIN, SIDE_Y_MAX = -1, -2, -3, -4
SIDE_NAMES = {SIDE_X_MIN: 'x_min', SIDE_X_MAX: 'x_max', SIDE_Y_MIN: 'y_min', SIDE_Y_MAX: 'y_max'}
TRACKING_MODES = ['labels', 'direction', 'straight']
STATUS_NAMES = ['matched by its cells', 'matched by 3 of its cells', 'moved with its neighbours',
                'moved along the fiber direction']


# =============================================================================
# Loading, thresholding and orienting the volume
# =============================================================================

def load_volume(file):
    """Load a 3D image: a .npy or .npz file (the first array), or a .tif/.tiff stack (needs tifffile)"""
    file = Path(file)
    suffix = file.suffix.lower()
    if suffix == '.npy':
        volume = np.load(file)
    elif suffix == '.npz':
        with np.load(file) as data:
            volume = data[data.files[0]]
    elif suffix in ('.tif', '.tiff'):
        try:
            import tifffile
        except ImportError as e:
            raise ImportError("reading TIFF stacks needs the optional package tifffile: pip install 'ringspy[image]'") from e
        volume = tifffile.imread(file)
    else:
        raise ValueError("unsupported volume file {}: use a .npy, .npz or .tif/.tiff file".format(file))
    volume = np.asarray(volume)
    if volume.ndim == 4 and volume.shape[-1] in (3, 4):  # RGB(A) stack: use the first channel
        volume = volume[..., 0]
    if volume.ndim != 3:
        raise ValueError("the volume must be a 3D array, {} has the shape {}".format(file, volume.shape))
    return volume


def otsu_threshold(values, nbins=256):
    """Otsu's threshold of an image: the gray value separating the two phases with the smallest
    intra-class variance (image histogram of nbins bins)"""
    values = np.asarray(values).reshape(-1)
    values = values[np.isfinite(values)]
    if values.size == 0:
        raise ValueError("no finite gray values to threshold")
    hist, edges = np.histogram(values, bins=nbins)
    hist = hist.astype(float)
    centers = 0.5*(edges[:-1] + edges[1:])
    w0 = np.cumsum(hist)
    w1 = w0[-1] - w0
    s0 = np.cumsum(hist*centers)
    s1 = s0[-1] - s0
    with np.errstate(divide='ignore', invalid='ignore'):
        between = w0*w1*(s0/w0 - s1/w1)**2
    between[~np.isfinite(between)] = -1.0
    return float(edges[int(np.argmax(between)) + 1])


def binarize(volume, solid='bright', threshold=None, smooth_sigma=0.0):
    """
    The cell walls of a 3D image as a boolean volume (True = wall).

    solid :: 'bright' if the walls are the bright phase (as in X-ray CT: the walls attenuate more
        than the empty lumens), 'dark' if they are the dark one.
    threshold :: gray value separating the phases; None: Otsu's threshold. Not used for a boolean image.
    smooth_sigma :: Gaussian smoothing (voxels) before thresholding, to suppress noise.

    Returns the boolean volume and the threshold used (None for a boolean image).
    """
    if solid not in ('bright', 'dark'):
        raise ValueError("solid must be 'bright' or 'dark', not {!r}".format(solid))
    if volume.dtype == bool:
        return (volume if solid == 'bright' else ~volume), None
    image = volume
    if smooth_sigma > 0:
        image = ndimage.gaussian_filter(np.asarray(volume, dtype=np.float32), smooth_sigma)
    if threshold is None:
        sample = image if image.size <= 20_000_000 else image[::2, ::2, ::2]
        threshold = otsu_threshold(sample)
    walls = image > threshold if solid == 'bright' else image < threshold
    return walls, float(threshold)


def detect_longitudinal_axis(walls, max_voxels=4_000_000):
    """The axis along which the image changes least, i.e. the fiber direction of the cells
    (the walls of the tube-like cells run along it)"""
    f = max(1, int(math.ceil((walls.size/max_voxels)**(1.0/3))))
    v = walls[::f, ::f, ::f]
    changes = [np.mean(np.diff(v, axis=axis) != 0) if v.shape[axis] > 1 else np.inf for axis in range(3)]
    return int(np.argmin(changes))


def _parse_axis(longitudinal_axis):
    """Axis index or None (auto) from 'auto', 'x'/'y'/'z', '0'/'1'/'2' or an int"""
    if longitudinal_axis is None:
        return None
    text = str(longitudinal_axis).strip().lower()
    if text == 'auto':
        return None
    names = {'x': 0, 'y': 1, 'z': 2, '0': 0, '1': 1, '2': 2}
    if text not in names:
        raise ValueError("longitudinal_axis must be 'auto', 0, 1, 2 or x, y, z, not {!r}".format(longitudinal_axis))
    return names[text]


# =============================================================================
# Cell labels, consistent from slice to slice
# =============================================================================

def _remove_small_components(mask, min_pixels):
    """mask without its connected components (4-connectivity) of fewer than min_pixels pixels"""
    if min_pixels <= 1 or not mask.any():
        return mask
    labels, n = ndimage.label(mask)
    if n == 0:
        return mask
    sizes = np.bincount(labels.reshape(-1))
    small = sizes < min_pixels
    small[0] = False
    return mask & ~small[labels]


def clean_slice(walls2d, min_lumen_pixels=0, min_wall_pixels=0):
    """A cross-section of the walls without lumens (void components) smaller than min_lumen_pixels
    (they become wall) and without wall specks smaller than min_wall_pixels (they become void)"""
    void = _remove_small_components(~walls2d, min_lumen_pixels)
    return _remove_small_components(~void, min_wall_pixels)


def _link_regions(current, ncurrent, history, next_label, min_overlap, dead=()):
    """
    Give the regions of the current cross-section (labels 1..ncurrent) the labels of the regions
    they overlap most in the previous cross-sections (history, most recent last; each previous
    label goes to one current region only, labels in `dead` to none); regions overlapping none get
    new labels.
    """
    assigned = np.zeros(ncurrent + 1, dtype=np.int64)
    taken = set(dead)
    area = np.bincount(current.reshape(-1), minlength=ncurrent + 1)
    for previous in reversed(history):
        pending = assigned == 0
        pending[0] = False
        if not pending.any():
            break
        mask = pending[current] & (previous > 0)
        if not mask.any():
            continue
        c = current[mask].astype(np.int64)
        p = previous[mask].astype(np.int64)
        pmax = int(p.max()) + 1
        keys, counts = np.unique(c*pmax + p, return_counts=True)
        cs, ps = keys // pmax, keys % pmax
        for idx in np.argsort(-counts, kind='stable'):
            ci, pi = int(cs[idx]), int(ps[idx])
            if assigned[ci] or pi in taken:
                continue
            if counts[idx] < min_overlap*area[ci]:
                continue
            assigned[ci] = pi
            taken.add(pi)
    new = np.nonzero(assigned[1:] == 0)[0] + 1
    assigned[new] = np.arange(next_label, next_label + len(new))
    return assigned[current], next_label + len(new)


def link_cell_labels(walls, k_ref, needed, min_lumen_pixels=0, min_wall_pixels=0, stride=1, memory=8,
                     min_overlap=0.25, closed_fraction=0.5):
    """
    Label the lumens (cells) of the cross-sections, keeping the label of a cell from slice to slice.

    The lumens of the reference slice k_ref are labelled first (4-connected void regions); going
    up and down from it, every slice's regions take the labels of the regions they overlap most
    in the previous `memory` labelled slices (so a label survives a slice in which two cells touch
    through a hole in their wall, or a cell is split by a speck), and new cells get new labels.
    A cell whose lumen closes (most of its pixels have become wall in the next slice: the cell
    ends) keeps its label for no later region, so the cell beginning above it is a new cell.
    Only every `stride`-th slice is labelled (the cells move little from slice to slice), always
    including the slices in `needed`.

    Returns {k: (labels, walls2d)} for the slices in `needed` (labels: 0 = wall, cells from 1), the
    number of labels, and {label: [first, last]}, the first and last labelled slice holding each cell
    (where a cell begins and ends inside the volume).
    """
    nz = walls.shape[2]
    needed = set(int(k) for k in needed) | {int(k_ref)}
    if min(needed) < 0 or max(needed) >= nz:
        raise ValueError("slice indices must be within 0 and {:d}".format(nz - 1))
    stride = max(int(stride), 1)
    walls_ref = clean_slice(walls[:, :, k_ref], min_lumen_pixels, min_wall_pixels)
    labels_ref, n = ndimage.label(~walls_ref)
    labels_ref = labels_ref.astype(np.int64)
    slices = {k_ref: (labels_ref, walls_ref)}
    next_label = n + 1
    ranges = {}

    def note(labels, k):
        for label in np.unique(labels).tolist():
            if label > 0:
                span = ranges.setdefault(label, [k, k])
                span[0] = min(span[0], k)
                span[1] = max(span[1], k)
    note(labels_ref, k_ref)
    for step in (1, -1):
        history = [labels_ref]
        dead = set()
        last = max(needed) if step > 0 else min(needed)
        ks = list(range(k_ref + step*stride, last, step*stride))
        if last != k_ref:
            ks.append(last)
        for k in sorted(needed):
            if (step > 0 and k_ref < k < last and k not in ks) or (step < 0 and last < k < k_ref and k not in ks):
                ks.append(k)
        ks = sorted(set(ks), reverse=(step < 0))
        for k in ks:
            walls2d = clean_slice(walls[:, :, k], min_lumen_pixels, min_wall_pixels)
            current, ncurrent = ndimage.label(~walls2d)
            labels, next_label = _link_regions(current, ncurrent, history, next_label, min_overlap, dead)
            # cells whose lumen has closed end here: their labels are not given to the cells beginning above
            previous = history[-1]
            for label in np.setdiff1d(np.unique(previous), np.unique(labels)).tolist():
                if label > 0 and label not in dead and walls2d[previous == label].mean() >= closed_fraction:
                    dead.add(label)
            if k in needed:
                slices[k] = (labels, walls2d)
            note(labels, k)
            history.append(labels)
            del history[:-memory]
    return slices, next_label - 1, ranges


# =============================================================================
# The network of cell corners and cell walls of a cross-section
# =============================================================================

@dataclass
class SliceNetwork:
    """
    Cell corners (points) and cell walls (ridges) of a cross-section, found as the junctions and
    the boundaries of the regions of the image closest to each lumen (the Voronoi-like partition
    of the cross-section into cells). The four sides of the cross-section count as cells with the
    pseudo labels SIDE_X_MIN .. SIDE_Y_MAX, so a wall ending on a side ends at a boundary point,
    and the corners of the cross-section are its corners.
    """
    points: np.ndarray                 # (n, 2) coordinates (mm)
    label_sets: List[frozenset]        # the cells (labels) meeting at each point
    kind: np.ndarray                   # 0: inside the cross-section, 1: on a side, 2: corner of the cross-section
    ridges: np.ndarray                 # (m, 2) point indices
    pairs: List[Tuple[int, int]]       # the two cells (a < b) each ridge separates
    thickness: np.ndarray              # (m,) wall thickness (mm)
    cell_area: Dict[int, float]        # area of each cell (mm2)
    box: Tuple[float, float, float, float]  # x_min, x_max, y_min, y_max
    index: Dict[frozenset, int] = field(init=False, repr=False)
    by_label: Dict[int, List[int]] = field(init=False, repr=False)
    pair_index: Dict[Tuple[int, int], List[int]] = field(init=False, repr=False)

    def __post_init__(self):
        self.index = {}
        self.by_label = {}
        for i, labels in enumerate(self.label_sets):
            self.index.setdefault(labels, i)
            for label in labels:
                self.by_label.setdefault(label, []).append(i)
        self.pair_index = {}
        for r, pair in enumerate(self.pairs):
            self.pair_index.setdefault(pair, []).append(r)

    @property
    def sides(self):
        """The sides (pseudo labels) each point lies on"""
        return [frozenset(label for label in labels if label < 0) for labels in self.label_sets]

    @property
    def median_cell_diameter(self):
        """Median equivalent diameter of the cells (mm)"""
        areas = np.array([a for label, a in self.cell_area.items() if label > 0])
        return float(2*math.sqrt(np.median(areas)/math.pi)) if len(areas) else 0.0


def slice_network(labels2d, walls2d, voxel, origin):
    """
    Build the SliceNetwork of a cross-section from its cell labels (0 = wall) and its walls.

    Every pixel is given the label of the nearest lumen; the corners of the pixel grid where 3 or
    more labels meet are junctions (adjacent junction pixels form one cell corner), and the
    boundary between two labels is a cell wall, as thick as twice the distance of its midline
    pixels to the lumens. voxel = (vx, vy) and origin = (x, y) of the corner of the first pixel.
    """
    nx, ny = labels2d.shape
    vx, vy = float(voxel[0]), float(voxel[1])
    x0, y0 = float(origin[0]), float(origin[1])
    x1, y1 = x0 + nx*vx, y0 + ny*vy
    void = labels2d > 0
    if not void.any():
        raise ValueError("the cross-section has no lumen (void) pixels: check `solid` and the threshold")
    if void.all():
        raise ValueError("the cross-section has no cell wall (solid) pixels: check `solid` and the threshold")

    # nearest lumen of every pixel, and half the wall thickness at every wall pixel
    half_thickness, (ii, jj) = ndimage.distance_transform_edt(~void, sampling=(vx, vy), return_indices=True)
    nearest = labels2d[ii, jj]
    padded = np.pad(nearest, 1).astype(np.int64)
    padded[0, :] = SIDE_X_MIN
    padded[-1, :] = SIDE_X_MAX
    padded[:, 0] = SIDE_Y_MIN
    padded[:, -1] = SIDE_Y_MAX

    # the 2x2 blocks of pixels around each corner (i, j) of the pixel grid, at x0 + i*vx, y0 + j*vy
    blocks = np.stack((padded[:-1, :-1], padded[1:, :-1], padded[:-1, 1:], padded[1:, 1:]), axis=-1)
    blocks.sort(axis=-1)
    distinct = 1 + np.count_nonzero(np.diff(blocks, axis=-1), axis=-1)
    junction = distinct >= 3
    clusters, ncluster = ndimage.label(junction, structure=np.ones((3, 3), dtype=int))
    if ncluster == 0:
        raise ValueError("no cell corner (junction of 3 cells) found in the cross-section")
    ids = np.arange(1, ncluster + 1)
    I, J = np.nonzero(junction)
    cluster_of = clusters[I, J]
    points = np.column_stack((x0 + ndimage.mean(I.astype(float), labels=cluster_of, index=ids)*vx,
                              y0 + ndimage.mean(J.astype(float), labels=cluster_of, index=ids)*vy))
    members = np.unique(np.column_stack((np.repeat(cluster_of, 4), blocks[I, J].reshape(-1))), axis=0)
    first = np.searchsorted(members[:, 0], ids, side='left')
    last = np.searchsorted(members[:, 0], ids, side='right')
    label_sets = [frozenset(members[a:b, 1].tolist()) for a, b in zip(first, last)]
    kind = np.array([min(sum(1 for label in labels if label < 0), 2) for labels in label_sets])
    for i, labels in enumerate(label_sets):  # boundary points lie exactly on their sides
        if SIDE_X_MIN in labels:
            points[i, 0] = x0
        if SIDE_X_MAX in labels:
            points[i, 0] = x1
        if SIDE_Y_MIN in labels:
            points[i, 1] = y0
        if SIDE_Y_MAX in labels:
            points[i, 1] = y1

    # wall pixels: grid corners between exactly two labels, grouped by the pair of labels
    I2, J2 = np.nonzero(distinct == 2)
    a2, b2 = blocks[I2, J2, 0], blocks[I2, J2, 3]
    padded_half = np.pad(half_thickness, 1)
    t2 = 2*(padded_half[I2, J2] + padded_half[I2 + 1, J2] + padded_half[I2, J2 + 1] + padded_half[I2 + 1, J2 + 1])/4
    order = np.lexsort((b2, a2))
    a2, b2, I2, J2, t2 = a2[order], b2[order], I2[order], J2[order], t2[order]
    change = np.ones(len(a2), dtype=bool)
    change[1:] = (a2[1:] != a2[:-1]) | (b2[1:] != b2[:-1])
    starts = np.nonzero(change)[0]
    stops = np.append(starts[1:], len(a2))

    pair_rows = {(int(a2[s0]), int(b2[s0])): (s0, s1) for s0, s1 in zip(starts, stops)}
    counts = np.bincount(labels2d.reshape(-1))
    cell_area = {int(label): float(counts[label]*vx*vy) for label in range(1, len(counts)) if counts[label] > 0}
    median_diameter = 2*math.sqrt(np.median(list(cell_area.values()))/math.pi) if cell_area else 0.0
    points = _refine_corners(points, label_sets, kind, pair_rows, x0 + I2*vx, y0 + J2*vy, (x0, x1, y0, y1),
                             radius=float(np.clip(median_diameter/3, 3*max(vx, vy), 12*max(vx, vy))),
                             max_shift=2*max(vx, vy))

    # the corners each pair of cells meets at
    pair_clusters = {}
    for i, labels in enumerate(label_sets):
        sorted_labels = sorted(labels)
        for u in range(len(sorted_labels)):
            for v in range(u + 1, len(sorted_labels)):
                pair_clusters.setdefault((sorted_labels[u], sorted_labels[v]), []).append(i)

    support_radius = 2.5*max(vx, vy)
    found = {}
    for s0, s1 in zip(starts, stops):
        pair = (int(a2[s0]), int(b2[s0]))
        corners = pair_clusters.get(pair)
        if corners is None or len(corners) < 2:
            continue
        thickness = _wall_thickness(t2[s0:s1], x0 + I2[s0:s1]*vx, y0 + J2[s0:s1]*vy, points[corners]) \
            if pair[0] > 0 else np.nan  # NaN: along a side
        if len(corners) == 2:
            segments = [(corners[0], corners[1])]
        else:  # the wall passes several corners: connect them in order along the wall, where there are wall pixels
            pts = points[corners]
            center = pts.mean(axis=0)
            axis = np.linalg.svd(pts - center, full_matrices=False)[2][0]
            ordered = [corners[k] for k in np.argsort((pts - center) @ axis, kind='stable')]
            segments = list(zip(ordered[:-1], ordered[1:]))
            wall_pixels = np.column_stack((x0 + I2[s0:s1]*vx, y0 + J2[s0:s1]*vy))
            mids = np.array([(points[a] + points[b])/2 for a, b in segments])
            distance = cKDTree(wall_pixels).query(mids)[0]
            segments = [seg for seg, d in zip(segments, distance) if d <= support_radius]
        for a, b in segments:
            if a == b or np.linalg.norm(points[a] - points[b]) <= 1e-12:
                continue
            key = (min(a, b), max(a, b))
            if key not in found:
                found[key] = (pair, thickness)
    if not found:
        raise ValueError("no cell wall between two cell corners found in the cross-section")
    ridges = np.array(list(found.keys()), dtype=int).reshape(-1, 2)
    pairs = [value[0] for value in found.values()]
    thickness = np.array([value[1] for value in found.values()], dtype=float)
    if np.all(np.isnan(thickness)):
        thickness[:] = max(vx, vy)
    else:
        typical = float(np.nanmedian(thickness))
        lengths = np.linalg.norm(points[ridges[:, 0]] - points[ridges[:, 1]], axis=1)
        thickness[np.isnan(thickness) | (lengths < thickness)] = typical  # walls along the sides, and walls too short to measure

    return SliceNetwork(points=points, label_sets=label_sets, kind=kind, ridges=ridges, pairs=pairs,
                        thickness=thickness, cell_area=cell_area, box=(x0, x1, y0, y1))


def _wall_thickness(estimates, wall_x, wall_y, corners):
    """
    Thickness of a wall from the estimates (twice the distance to the lumens) at its midline pixels
    (wall_x, wall_y): the median over the pixels away from the corners the wall ends at, where the
    midline runs through the fillet of the junction and reads thick; for a wall too short for that,
    the 30th percentile of all its pixels.
    """
    rough = float(np.percentile(estimates, 30))
    if len(corners) and len(estimates) >= 6:
        pixels = np.column_stack((wall_x, wall_y))
        distance = np.min(np.linalg.norm(pixels[:, np.newaxis, :] - np.asarray(corners)[np.newaxis, :, :], axis=2), axis=1)
        away = distance > 0.75*rough
        if np.count_nonzero(away) >= 3:
            return float(np.median(estimates[away]))
    return rough


def _refine_corners(points, label_sets, kind, pair_rows, wall_x, wall_y, box, radius, max_shift):
    """
    Sub-pixel positions of the cell corners: a corner is where the midlines of its walls meet. A line
    is fitted to the wall pixels (wall_x, wall_y; pair_rows: the rows of each pair of cells) within
    `radius` of the corner, for each of its walls, and the corner is moved to the point closest to
    these lines (and to the side of the box, for a boundary point), if that is within max_shift of
    the pixel position. The corners of the box stay where they are.
    """
    x0, x1, y0, y1 = box
    side_lines = {SIDE_X_MIN: ((x0, 0.0), (0.0, 1.0)), SIDE_X_MAX: ((x1, 0.0), (0.0, 1.0)),
                  SIDE_Y_MIN: ((0.0, y0), (1.0, 0.0)), SIDE_Y_MAX: ((0.0, y1), (1.0, 0.0))}
    refined = points.copy()
    identity = np.eye(2)
    for i, labels in enumerate(label_sets):
        if kind[i] == 2:
            continue
        center = points[i]
        A = np.zeros((2, 2))
        b = np.zeros(2)
        nlines = 0
        sorted_labels = sorted(labels)
        for u in range(len(sorted_labels)):
            for v in range(u + 1, len(sorted_labels)):
                if sorted_labels[u] < 0:
                    continue
                rows = pair_rows.get((sorted_labels[u], sorted_labels[v]))
                if rows is None:
                    continue
                s0, s1 = rows
                dx = wall_x[s0:s1] - center[0]
                dy = wall_y[s0:s1] - center[1]
                near = dx*dx + dy*dy <= radius*radius
                if np.count_nonzero(near) < 3:
                    continue
                pts = np.column_stack((wall_x[s0:s1][near], wall_y[s0:s1][near]))
                mean = pts.mean(axis=0)
                cov = (pts - mean).T @ (pts - mean)
                w, vec = np.linalg.eigh(cov)
                if w[1] <= 1e-12 or w[0] > 0.25*w[1]:  # no clear line (a wall pixel blob)
                    continue
                direction = vec[:, 1]
                projector = identity - np.outer(direction, direction)
                A += projector
                b += projector @ mean
                nlines += 1
        for label in sorted_labels:
            if label in side_lines:
                point, direction = side_lines[label]
                projector = identity - np.outer(direction, direction)
                A += projector
                b += projector @ np.asarray(point)
                nlines += 1
        if nlines < 2 or abs(np.linalg.det(A)) < 1e-9:
            continue
        x = np.linalg.solve(A, b)
        if np.linalg.norm(x - center) <= max_shift:
            refined[i] = x
    for i, labels in enumerate(label_sets):  # boundary points exactly on their sides
        if SIDE_X_MIN in labels:
            refined[i, 0] = x0
        if SIDE_X_MAX in labels:
            refined[i, 0] = x1
        if SIDE_Y_MIN in labels:
            refined[i, 1] = y0
        if SIDE_Y_MAX in labels:
            refined[i, 1] = y1
    return refined


def merge_short_ridges(net, tol):
    """
    The network with its ridges shorter than tol merged into one point (the merge_operation of
    the mesh generation): the merged point lies at the mean of the inside points, stays on the side
    for a boundary point and at the corner for a corner; points on different sides, and two
    corners, are not merged. The merged point belongs to all the cells of the points merged.
    """
    n = len(net.points)
    parent = np.arange(n)

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    points = net.points.copy()
    kind = net.kind.copy()
    sets = list(net.label_sets)
    sides = net.sides
    count = np.ones(n)
    lengths = np.linalg.norm(points[net.ridges[:, 0]] - points[net.ridges[:, 1]], axis=1)
    for r in np.argsort(lengths, kind='stable'):
        if lengths[r] >= tol:
            break
        u, v = find(net.ridges[r, 0]), find(net.ridges[r, 1])
        if u == v or np.linalg.norm(points[u] - points[v]) >= tol:
            continue
        if kind[u] < kind[v]:
            u, v = v, u
        if kind[u] == 2 and kind[v] == 2:
            continue
        if kind[u] == 2 and kind[v] == 1 and not sides[v] <= sides[u]:
            continue
        if kind[u] == 1 and kind[v] == 1 and sides[u] != sides[v]:
            continue
        if kind[u] == 0 or (kind[u] == 1 and kind[v] == 1):
            points[u] = (points[u]*count[u] + points[v]*count[v])/(count[u] + count[v])
        count[u] += count[v]
        sets[u] = sets[u] | sets[v]
        sides[u] = sides[u] | sides[v]
        parent[v] = u

    roots = np.array([find(i) for i in range(n)])
    kept = np.unique(roots)
    new_index = np.full(n, -1)
    new_index[kept] = np.arange(len(kept))
    ridge_of = new_index[roots[net.ridges]]
    found = {}
    for r, (a, b) in enumerate(ridge_of):
        if a == b:
            continue
        key = (min(a, b), max(a, b))
        if key not in found:
            found[key] = r
    rows = list(found.values())
    return SliceNetwork(points=points[kept], label_sets=[sets[i] for i in kept], kind=kind[kept],
                        ridges=np.array(list(found.keys()), dtype=int).reshape(-1, 2),
                        pairs=[net.pairs[r] for r in rows], thickness=net.thickness[rows],
                        cell_area=net.cell_area, box=net.box)


@dataclass
class OrderedNetwork:
    """The reference network in the order of the mesh generation: inside vertices, then boundary points,
    then the corners of the box; finite ridges (both ends inside), then the cut ridges (one end on the
    boundary) and the ridges between boundary points (along the sides when the boundary walls are on)"""
    vertices: np.ndarray
    label_sets: List[frozenset]
    kind: np.ndarray
    sides: List[frozenset]
    nvertices_in: int
    nboundary_pts: int
    nboundary_pts_featured: int
    finite_ridges: np.ndarray
    boundary_ridges: np.ndarray
    ridges: np.ndarray
    pairs: List[Tuple[int, int]]
    thickness: np.ndarray
    box: Tuple[float, float, float, float]

    @property
    def nvertex(self):
        return len(self.vertices)

    @property
    def nridge(self):
        return len(self.ridges)


def order_network(net, boundary_on=True):
    """
    Order a SliceNetwork as the mesh generation orders a clipped Voronoi diagram (see OrderedNetwork).
    Without the boundary walls (boundary_on False), the ridges along the sides, the ridges between two
    boundary points, the corners of the box and the boundary points left without a ridge are dropped.
    """
    kind = net.kind
    ridges = net.ridges
    ends_kind = kind[ridges]
    keep_ridge = np.ones(len(ridges), dtype=bool)
    if not boundary_on:
        keep_ridge = ~((ends_kind[:, 0] >= 1) & (ends_kind[:, 1] >= 1))
    ridges = ridges[keep_ridge]
    used = np.zeros(len(kind), dtype=bool)
    used[ridges.reshape(-1)] = True
    keep_vertex = (kind == 0) | used
    if not boundary_on:
        keep_vertex &= kind < 2
    perm = np.concatenate([np.nonzero(keep_vertex & (kind == k))[0] for k in (0, 1, 2)])
    inverse = np.full(len(kind), -1)
    inverse[perm] = np.arange(len(perm))
    nvertices_in = int(np.count_nonzero(keep_vertex & (kind == 0)))
    nboundary = int(np.count_nonzero(keep_vertex & (kind == 1)))
    nfeatured = int(np.count_nonzero(keep_vertex & (kind == 2)))
    if boundary_on and nfeatured != 4:
        raise ValueError("expected the 4 corners of the cross-section among the cell corners, found {:d}".format(nfeatured))
    new_ridges = inverse[ridges]
    inside = new_ridges < nvertices_in
    finite = inside[:, 0] & inside[:, 1]
    cut = inside[:, 0] ^ inside[:, 1]
    other = ~(finite | cut)
    if not finite.any():
        raise ValueError("no cell wall with both ends inside the cross-section: the cross-section holds no "
                         "complete cell (too small a crop, or the cells were not separated by the threshold)")
    order = np.concatenate([np.nonzero(finite)[0], np.nonzero(cut)[0], np.nonzero(other)[0]])
    rows = np.nonzero(keep_ridge)[0][order]
    sides = net.sides
    return OrderedNetwork(vertices=net.points[perm].copy(), label_sets=[net.label_sets[i] for i in perm],
                          kind=kind[perm].copy(), sides=[sides[i] for i in perm], nvertices_in=nvertices_in,
                          nboundary_pts=nboundary + nfeatured, nboundary_pts_featured=nfeatured,
                          finite_ridges=new_ridges[finite], boundary_ridges=np.vstack((new_ridges[cut], new_ridges[other])),
                          ridges=new_ridges[order], pairs=[net.pairs[r] for r in rows],
                          thickness=net.thickness[rows].copy(), box=net.box)


# =============================================================================
# The fiber direction field
# =============================================================================

class DirectionField:
    """
    The local fiber direction of the volume from its structure tensor: at every point, the unit
    vector (pointing to +z) along which the image changes least, i.e. the direction of the
    tube-like cells. Computed on a coarse grid (the volume block-averaged to at most max_voxels
    voxels), with gradients at the scale of one coarse voxel and the tensor averaged over `scale`
    (mm, about two cell diameters so that the walls of a few cells are averaged).

    mean_direction, grain_angle_deg (its angle to the z axis), azimuth_deg (the direction it leans
    to, from the x axis) and dispersion_deg (the mean angle of the local directions to the mean
    direction) describe how the fibers are aligned.
    """

    def __init__(self, walls, voxel, origin, scale, max_voxels=3_000_000, sigma=1.0):
        shape = np.array(walls.shape)
        factor = max(1, int(math.ceil((walls.size/max_voxels)**(1.0/3)))) if walls.size > max_voxels else 1
        while factor > 1 and np.any(shape//factor < 3):
            factor -= 1
        cropped = (shape//factor)*factor
        v = walls[:cropped[0], :cropped[1], :cropped[2]].astype(np.float32)
        if factor > 1:
            v = v.reshape(cropped[0]//factor, factor, cropped[1]//factor, factor, cropped[2]//factor, factor).mean(axis=(1, 3, 5))
        if min(v.shape) < 3:
            raise ValueError("the volume is too small (fewer than 3 voxels along an axis) for a fiber direction field")
        self.factor = factor
        self.spacing = np.asarray(voxel, dtype=float)*factor
        self.origin = np.asarray(origin, dtype=float)
        self.scale = float(scale)

        gradients = [ndimage.gaussian_filter(v, sigma, order=tuple(int(a == axis) for a in range(3)))/self.spacing[axis]
                     for axis in range(3)]
        rho = np.maximum(self.scale/self.spacing, 1.0)
        tensor = np.empty(v.shape + (3, 3), dtype=np.float32)
        for a in range(3):
            for b in range(a, 3):
                tensor[..., a, b] = ndimage.gaussian_filter(gradients[a]*gradients[b], rho)
                tensor[..., b, a] = tensor[..., a, b]
        del gradients
        flat = tensor.reshape(-1, 3, 3)
        direction = np.empty((flat.shape[0], 3), dtype=np.float32)
        anisotropy = np.empty(flat.shape[0], dtype=np.float32)
        chunk = 400_000
        for s in range(0, flat.shape[0], chunk):
            w, vec = np.linalg.eigh(flat[s:s + chunk].astype(np.float64))
            direction[s:s + chunk] = vec[:, :, 0]  # eigenvector of the smallest eigenvalue: least change along it
            anisotropy[s:s + chunk] = (w[:, 1] - w[:, 0])/(w[:, 1] + w[:, 0] + 1e-30)
        del tensor, flat
        direction *= np.where(direction[:, 2] < 0, -1.0, 1.0)[:, np.newaxis]
        weights = np.maximum(anisotropy, 1e-6).astype(np.float64)
        mean = (direction.astype(np.float64)*weights[:, np.newaxis]).sum(axis=0)
        self.mean_direction = mean/np.linalg.norm(mean)
        cosine = np.clip(direction.astype(np.float64) @ self.mean_direction, -1.0, 1.0)
        self.dispersion_deg = float(math.degrees(np.average(np.arccos(cosine), weights=weights)))
        self.grain_angle_deg = float(math.degrees(math.acos(min(1.0, self.mean_direction[2]))))
        self.azimuth_deg = float(math.degrees(math.atan2(self.mean_direction[1], self.mean_direction[0])))
        self.mean_anisotropy = float(anisotropy.mean())
        self.direction = direction.reshape(v.shape + (3,))
        self.anisotropy = anisotropy.reshape(v.shape)

    def at(self, xyz):
        """Unit fiber directions at the points xyz (n, 3), interpolated on the coarse grid"""
        xyz = np.asarray(xyz, dtype=float).reshape(-1, 3)
        coords = ((xyz - self.origin)/self.spacing - 0.5).T
        t = np.stack([ndimage.map_coordinates(self.direction[..., c], coords, order=1, mode='nearest')
                      for c in range(3)], axis=1).astype(float)
        norm = np.linalg.norm(t, axis=1, keepdims=True)
        norm[norm == 0] = 1.0
        return t/norm

    def advect(self, xy, z_from, z_to, step=None):
        """Move the points xy (n, 2) from height z_from to z_to along the fiber directions"""
        xy = np.array(xy, dtype=float, copy=True).reshape(-1, 2)
        if z_to == z_from or len(xy) == 0:
            return xy
        if step is None:
            step = 2*self.spacing[2]
        n = max(1, int(math.ceil(abs(z_to - z_from)/step)))
        dz = (z_to - z_from)/n
        z = z_from
        for _ in range(n):
            t = self.at(np.column_stack((xy, np.full(len(xy), z + dz/2))))
            tz = np.maximum(t[:, 2], 0.2)  # fibers are never flatter than ~78 degrees from the plane here
            xy += t[:, 0:2]/tz[:, np.newaxis]*dz
            z += dz
        return xy


# =============================================================================
# Following the cell corners through the layers
# =============================================================================

def _project_to_sides(points, ref, rows=None):
    """Put the boundary points of `ref` (rows: which ones, default all) back on their sides, and the
    corners at the corners"""
    x0, x1, y0, y1 = ref.box
    rows = range(len(points)) if rows is None else rows
    for i in rows:
        sides = ref.sides[i]
        if not sides:
            continue
        if SIDE_X_MIN in sides:
            points[i, 0] = x0
        if SIDE_X_MAX in sides:
            points[i, 0] = x1
        if SIDE_Y_MIN in sides:
            points[i, 1] = y0
        if SIDE_Y_MAX in sides:
            points[i, 1] = y1
    return points


def _match_corners(ref, net, predicted, max_jump):
    """
    Position of every reference vertex in the slice network `net`: the junction of the same cells
    (exact match: status 0), or the junction sharing at least 3 of its cells closest to the predicted
    position (status 1), within max_jump of the prediction. A junction takes one vertex only (the
    one matching it best, then the closest one): two cell corners that have come together at this
    height stay two nodes. Unmatched vertices keep the predicted position with status 3; the corners
    of the box keep their place (status 0).
    """
    n = len(predicted)
    new = predicted.copy()
    status = np.full(n, 3, dtype=int)
    proposals = []  # (-shared cells, distance, vertex, junction)
    for i in range(n):
        if ref.kind[i] == 2:
            new[i] = ref.vertices[i]
            status[i] = 0
            continue
        labels = ref.label_sets[i]
        candidates = Counter()
        for label in labels:
            for j in net.by_label.get(label, ()):
                candidates[j] += 1
        if not candidates:
            continue
        best = max(candidates.values())
        if best < 3:
            continue
        js = [j for j, c in candidates.items() if c == best]
        distance = np.linalg.norm(net.points[js] - predicted[i], axis=1)
        for j, d in sorted(zip(js, distance), key=lambda jd: jd[1]):
            if d <= max_jump:
                proposals.append((-best, d, i, j))
    taken = set()
    for shared, d, i, j in sorted(proposals):
        if status[i] <= 1 or j in taken:
            continue
        new[i] = net.points[j]
        status[i] = 0 if -shared == len(ref.label_sets[i]) else 1
        taken.add(j)
    return new, status


def track_corners(ref, networks, slice_z, k_ref, field=None, mode='labels', max_jump=np.inf):
    """
    Follow the vertices of the reference network through the cross-sections in `networks`
    ({k: SliceNetwork}), from the reference slice k_ref up and down.

    From slice to slice, each vertex is first moved along the fiber direction field (if given) and
    then, in mode 'labels', put at the junction of the same cells in the new slice (see
    _match_corners); a vertex whose junction is not found (a cell ended, a wall was not resolved)
    moves with its neighbours in the network instead. Mode 'direction' only moves the vertices
    along the fiber directions; mode 'straight' keeps the reference positions.

    Returns {k: positions (n, 2)} and {k: status (n,)}, the status of each vertex in each slice
    being the index into STATUS_NAMES.
    """
    if mode not in TRACKING_MODES:
        raise ValueError("tracking must be one of {}, not {!r}".format(', '.join(TRACKING_MODES), mode))
    n = ref.nvertex
    positions = {k_ref: ref.vertices.copy()}
    status = {k_ref: np.zeros(n, dtype=int)}
    neighbours = [[] for _ in range(n)]
    for a, b in ref.ridges:
        neighbours[a].append(b)
        neighbours[b].append(a)
    ks = sorted(networks)
    chains = [[k for k in ks if k > k_ref], [k for k in reversed(ks) if k < k_ref]]
    for chain in chains:
        k_prev = k_ref
        for k in chain:
            previous = positions[k_prev]
            if mode == 'straight':
                new, st = ref.vertices.copy(), np.full(n, 3, dtype=int)
            else:
                predicted = field.advect(previous, slice_z[k_prev], slice_z[k]) if field is not None else previous.copy()
                predicted = _project_to_sides(predicted, ref)
                if mode == 'direction':
                    new, st = predicted, np.full(n, 3, dtype=int)
                    new[ref.kind == 2] = ref.vertices[ref.kind == 2]
                    st[ref.kind == 2] = 0
                else:
                    new, st = _match_corners(ref, networks[k], predicted, max_jump)
                    # vertices without a junction move with their matched neighbours
                    matched = st <= 1
                    displacement = new - previous
                    for _ in range(2):
                        for i in np.nonzero(~matched)[0]:
                            rows = [j for j in neighbours[i] if matched[j] or st[j] == 2]
                            if len(rows) >= 2:
                                new[i] = previous[i] + displacement[rows].mean(axis=0)
                                st[i] = 2
                        displacement = new - previous
                    new = _project_to_sides(new, ref, np.nonzero(st >= 2)[0])
            positions[k] = new
            status[k] = st
            k_prev = k
    return positions, status


def fit_midpoints(layer_vertices, nsegments, degree=NURBS_DEGREE):
    """
    Move the middle control point of every quadratic Bezier element of the beams so that the beam
    passes through the position followed in the image at the middle layer (a quadratic Bezier
    curve passes through its end control points but not through the middle one: P(1/2) =
    (P0 + 2*P1 + P2)/4). For other degrees the control points are left at the followed positions.
    """
    layer_vertices = np.array(layer_vertices, dtype=float, copy=True)
    if degree != 2:
        return layer_vertices
    nctrlpt_per_beam = 2*degree + 1
    for seg in range(nsegments):
        base = seg*nctrlpt_per_beam
        for mid in (1, 3):
            layer_vertices[base + mid] = 2*layer_vertices[base + mid] - (layer_vertices[base + mid - 1] +
                                                                           layer_vertices[base + mid + 1])/2
    return layer_vertices


def fit_midpoints_grid(layer_vertices):
    """fit_midpoints for positions on the half-layers of a grid of quadratic beam elements (odd half-layers
    hold the middle control points of the elements)"""
    layer_vertices = np.array(layer_vertices, dtype=float, copy=True)
    for m in range(1, len(layer_vertices) - 1, 2):
        layer_vertices[m] = 2*layer_vertices[m] - (layer_vertices[m - 1] + layer_vertices[m + 1])/2
    return layer_vertices


def fiber_ends_from_image(ref, label_ranges, nz, z_min, vz, h, nelem):
    """
    The grid lines (1 .. nelem-1, spacing h from z_min) where the fibers of each vertex of the reference
    network end: wherever one of the cells meeting at the vertex begins or ends inside the volume
    (label_ranges: the first and last of the nz slices holding each cell, from link_cell_labels).
    """
    ends = []
    for labels in ref.label_sets:
        ks = set()
        for label in labels:
            if label < 0 or label not in label_ranges:
                continue
            first, last = label_ranges[label]
            heights = []
            if first > 0:
                heights.append(z_min + first*vz)
            if last < nz - 1:
                heights.append(z_min + (last + 1)*vz)
            for z in heights:
                k = int(round((z - z_min)/h))
                if 0 < k < nelem:
                    ks.add(k)
        ends.append(np.array(sorted(ks), dtype=int))
    return ends


def connector_tangents(IGAvertices, connectivity, npt_per_layer, nctrlpt_per_beam):
    """Unit tangents (along the grains) of the transverse connectors with the given 1-based node
    pairs: the mean of the directions of the two grains at their nodes (between the neighbouring
    control points of the same segment)"""
    nodes = np.asarray(connectivity, dtype=int) - 1
    layer = nodes // npt_per_layer
    grain = nodes % npt_per_layer
    first = (layer // nctrlpt_per_beam)*nctrlpt_per_beam
    lo = np.maximum(layer - 1, first)
    hi = np.minimum(layer + 1, first + nctrlpt_per_beam - 1)
    t = (IGAvertices[hi*npt_per_layer + grain] - IGAvertices[lo*npt_per_layer + grain]).sum(axis=1)
    norm = np.linalg.norm(t, axis=1, keepdims=True)
    norm[norm == 0] = 1.0
    return t/norm


# =============================================================================
# The reconstruction
# =============================================================================

@dataclass
class ReconstructionResult:
    """What reconstruct() made: the output folder and files, and the properties of the lattice and of the image"""
    config: object                   # the ReconstructionConfig, with the values derived from the image filled in
    output_dir: Path
    files: List[Path]
    shape: Tuple[int, int, int]      # voxels of the (oriented, cropped) volume
    voxel_size: Tuple[float, float, float]
    origin: Tuple[float, float, float]
    threshold: Optional[float]
    longitudinal_axis: int           # axis of the input array taken as the fiber direction
    nParticles: int                  # cells in the reference cross-section
    nbeamElem: int
    mass: float
    bulk_volume: float
    bulk_density: float
    porosity: float                  # of the lattice (from the connector volumes)
    image_porosity: float            # void fraction of the thresholded image
    grain_angle: float               # mean angle of the fibers to the longitudinal axis (degrees)
    grain_azimuth: float             # direction the fibers lean to, from the x axis (degrees)
    grain_dispersion: float          # mean angle of the local fiber directions to the mean direction (degrees)
    grain_inclination_mean: float    # mean inclination of the grains followed through the image (degrees)
    grain_inclination_max: float
    tracking_counts: Dict[str, int]  # how the vertices were placed in the layers (STATUS_NAMES)
    median_cell_diameter: float
    median_wall_thickness: float
    layer_vertices: np.ndarray       # (nlayers, nvertex, 2) in-plane positions of the vertices in every layer
    z_coord: np.ndarray              # heights of the layers
    vertices: np.ndarray             # (nvertex, 2) of the reference cross-section
    ridges: np.ndarray               # (nridge, 2)
    fiber_ends: Optional[list] = None  # grid lines where the fibers of each grain end (fiber_ends != 'none')
    nfibers: int = 0


def reconstruct(config, **changes):
    """
    Reconstruct the lattice of a cellular specimen from a 3D image and write the mesh files.

    config :: a ReconstructionConfig, a TOML file name or a dict of parameters, or the volume itself
        (a 3D numpy array, or a ringspy.synthetic.SyntheticVolume) with the parameters as keyword
        arguments: ``reconstruct(volume, voxel_size=0.002, geoName='scan')``.
    changes :: parameters to override.

    The files go to config.path/config.geoName/, which is deleted first if it exists.
    Returns a ReconstructionResult.
    """
    from .config import ReconstructionConfig
    volume = None
    if isinstance(config, np.ndarray):
        volume = config
        config = ReconstructionConfig(**changes)
        changes = {}
    elif isinstance(config, SyntheticVolume):
        volume = config.volume
        params = {'voxel_size': list(config.voxel_size), 'origin': list(config.origin), 'longitudinal_axis': '2'}
        params.update(changes)
        config = ReconstructionConfig(**params)
        changes = {}
    elif isinstance(config, (str, os.PathLike)) and Path(config).suffix.lower() in ('.npy', '.npz', '.tif', '.tiff'):
        config = ReconstructionConfig(volume=os.fspath(config), **changes)
        changes = {}
    elif isinstance(config, (str, os.PathLike)):
        config = ReconstructionConfig.from_toml(config)
    elif isinstance(config, dict):
        config = ReconstructionConfig.from_dict(config)
    elif not isinstance(config, ReconstructionConfig):
        raise TypeError("config must be a ReconstructionConfig, a TOML file name, a dict or a volume, not {}"
                        .format(type(config).__name__))
    if changes:
        config = config.replace(**changes)
    if volume is None and config.volume is None:
        raise ValueError("no volume to reconstruct: set `volume` to a .npy/.npz/.tif file, or pass the array")

    import matplotlib.pyplot as plt
    open_figures = set(plt.get_fignums())
    try:
        return _reconstruct(config, volume)
    finally:
        for num in set(plt.get_fignums()) - open_figures:
            plt.close(num)


def _reconstruct(c, volume):
    startTime = time.time()
    geoName, path = c.geoName, c.path
    output_dir = Path(path, geoName)
    shutil.rmtree(output_dir, ignore_errors=True)
    os.makedirs(output_dir, exist_ok=True)

    # ==================================================================
    # Load and threshold the volume, put the fiber direction along the last axis, crop
    if volume is None:
        volume = load_volume(c.volume)
        print('Volume {} loaded: shape {}, dtype {}'.format(c.volume, volume.shape, volume.dtype))
    volume = np.asarray(volume)
    if volume.ndim != 3:
        raise ValueError("the volume must be a 3D array, not of shape {}".format(volume.shape))
    walls, threshold = binarize(volume, c.solid, c.threshold, c.smooth_sigma)
    del volume
    axis = _parse_axis(c.longitudinal_axis)
    if axis is None:
        axis = detect_longitudinal_axis(walls)
        print('Longitudinal (fiber) axis of the volume: {:d}'.format(axis))
    walls = np.moveaxis(walls, axis, 2)
    voxel = np.asarray(c.voxel_size, dtype=float).reshape(3)
    voxel = np.array([voxel[a] for a in [i for i in range(3) if i != axis] + [axis]])
    if c.crop is not None:
        x0, x1, y0, y1, z0, z1 = c.crop
        walls = walls[x0:x1, y0:y1, z0:z1]
    walls = np.ascontiguousarray(walls)
    nx, ny, nz = walls.shape
    if min(nx, ny) < 4 or nz < 1:
        raise ValueError("the (cropped) volume is too small: shape {}".format(walls.shape))
    origin = np.asarray(c.origin, dtype=float).reshape(3)
    x_min, y_min, z_min = origin
    x_max, y_max, z_max = origin + voxel*np.array([nx, ny, nz])
    image_porosity = 1.0 - float(walls.mean())
    loadTime = time.time()
    print('Volume thresholded: {:d} x {:d} x {:d} voxels, porosity {:.3f}, in {:.3f} seconds'
          .format(nx, ny, nz, image_porosity, loadTime - startTime))

    # ==================================================================
    # Layers of control points and the cross-sections they need
    nctrlpt_per_beam = 2*NURBS_DEGREE + 1
    fibers_on = c.fiber_ends != 'none'
    if fibers_on:  # a grid of 2*nsegments beam elements, nodes on its half-layers; fibers end on its lines
        nelem = 2*c.nsegments
        nlayers = 2*nelem + 1
        h = (z_max - z_min)/nelem
        z_coord = z_min + np.arange(nlayers)*h/2
        segment_length = 2*h
        theta = np.zeros(nlayers)
    else:          # the regular segments, with the layers at the segment interfaces repeated
        theta, z_coord, nlayers, segment_length = rpgen.LayerHeights(NURBS_DEGREE, c.nsegments, 0.0, 0.0, z_min, z_max,
                                                                     c.long_connector_ratio)
        nelem, h = 2*c.nsegments, segment_length/2
    layer_slices = np.clip(np.floor((z_coord - z_min)/voxel[2]).astype(int), 0, nz - 1)
    k_ref = _reference_slice(c.reference_slice, nz, z_min, voxel[2])
    needed = set(layer_slices.tolist()) | {k_ref}
    slice_z = {k: z_min + (k + 0.5)*voxel[2] for k in needed}

    # ==================================================================
    # Cell labels of the cross-sections, consistent from slice to slice
    pixel_area = voxel[0]*voxel[1]
    min_lumen_px = int(round(c.min_cell_area/pixel_area)) if c.min_cell_area is not None else 9
    min_wall_px = int(round(c.min_wall_area/pixel_area)) if c.min_wall_area is not None else 9
    slices, nlabels, label_ranges = link_cell_labels(walls, k_ref, needed, min_lumen_px, min_wall_px, stride=c.link_stride)
    labelTime = time.time()
    print('{:d} cells labelled in {:d} cross-sections in {:.3f} seconds'.format(nlabels, len(slices), labelTime - loadTime))

    # ==================================================================
    # The network of cell corners and walls of every needed cross-section; the reference one is the 2D mesh
    networks = {k: slice_network(labels2d, walls2d, voxel[0:2], origin[0:2]) for k, (labels2d, walls2d) in slices.items()}
    ref_net = networks[k_ref]
    median_cell_diameter = ref_net.median_cell_diameter
    nParticles = sum(1 for label in ref_net.cell_area if label > 0)
    merge_tol = c.merge_tol if c.merge_tol is not None else 1.5*max(voxel[0], voxel[1])
    if c.merge_operation:
        ref_net = merge_short_ridges(ref_net, merge_tol)
    ref = order_network(ref_net, c.boundaryFlag)
    if c.wall_thickness is not None:
        ref.thickness[:] = c.wall_thickness
    median_wall_thickness = float(np.median(ref.thickness))
    networkTime = time.time()
    print('Reference cross-section (slice {:d}): {:d} cells, {:d} cell corners, {:d} cell walls, in {:.3f} seconds'
          .format(k_ref, nParticles, ref.nvertex, ref.nridge, networkTime - labelTime))

    # ==================================================================
    # Fiber direction field and the positions of the cell corners in every layer
    direction_scale = c.direction_scale if c.direction_scale is not None else max(2*median_cell_diameter, 3*max(voxel))
    field = None
    if c.tracking != 'straight':
        field = DirectionField(walls, voxel, origin, direction_scale)
        print('Fiber direction field: mean grain angle {:.2f} deg (azimuth {:.1f} deg), dispersion {:.2f} deg'
              .format(field.grain_angle_deg, field.azimuth_deg, field.dispersion_deg))
    max_jump = c.max_jump if c.max_jump is not None else 1.5*max(median_cell_diameter, 2*max(voxel))
    positions, status = track_corners(ref, networks, slice_z, k_ref, field, c.tracking, max_jump)
    layer_vertices = np.array([positions[k] for k in layer_slices])
    layer_status = np.array([status[k] for k in layer_slices])
    tracked_vertices = layer_vertices.copy()
    layer_vertices = fit_midpoints_grid(layer_vertices) if fibers_on else fit_midpoints(layer_vertices, c.nsegments, NURBS_DEGREE)
    followed = layer_status[layer_slices != k_ref][:, ref.kind < 2]  # the layers of the reference slice are trivial
    tracking_counts = {name: int(np.count_nonzero(followed == s)) for s, name in enumerate(STATUS_NAMES)}
    chord = tracked_vertices[-1] - tracked_vertices[0]
    height = max(z_coord[-1] - z_coord[0], 1e-300)
    inclination = np.degrees(np.arctan(np.linalg.norm(chord, axis=1)/height))[ref.kind < 2]
    trackTime = time.time()
    print('Cell corners followed through {:d} layers in {:.3f} seconds'.format(nlayers, trackTime - networkTime))

    # ==================================================================
    # Wall thickness of every transverse connector, measured in its layer
    nconnector_t_per_beam = int((nctrlpt_per_beam - 1)/NURBS_DEGREE + 1)

    def widths_in_layer(layer):
        if c.wall_thickness is not None:
            return np.full(ref.nridge, float(c.wall_thickness))
        return _widths_in_slice(ref, networks[layer_slices[layer]])

    # ==================================================================
    # The lattice: the same steps and files as the mesh generation
    boundary_points = np.array([[x_min, y_min], [x_max, y_min], [x_max, y_max], [x_min, y_max]])
    [all_pts_2D, all_ridges, npt_per_layer, npt_per_layer_normal, npt_per_layer_vtk] = \
        rpgen.RidgeMidQuarterPts(ref.vertices, ref.nvertex, ref.nvertices_in, ref.ridges, ref.finite_ridges,
                                 ref.boundary_ridges, len(ref.finite_ridges), len(ref.boundary_ridges),
                                 ref.nboundary_pts, ref.nboundary_pts_featured)
    [all_vertices_2D, max_wings, flattened_all_vertices_2D, all_ridges] = \
        rpgen.VertexandRidgeinfo(all_pts_2D, all_ridges, npt_per_layer, npt_per_layer_normal, npt_per_layer_vtk,
                                 ref.nridge, geoName, np.array([0.0]), (0.0, 0.0), median_wall_thickness,
                                 median_wall_thickness, path=path, ridge_widths=ref.thickness)
    layer_points = rpgen.LayerPoints(layer_vertices, all_ridges, npt_per_layer_vtk)
    fiber_ends = None
    fibers = None
    if fibers_on:
        if c.fiber_ends == 'image':
            fiber_ends = fiber_ends_from_image(ref, label_ranges, nz, z_min, voxel[2], h, nelem)
        else:
            fiber_ends = rplat.fiber_end_layers(ref.nvertex, nelem, h, c.fiber_length, c.fiber_length_cov, rng=c.seed)
        widths = np.array([widths_in_layer(2*k) for k in range(nelem + 1)])
        lat = rplat.build_lattice(all_ridges, ref.nvertex, ref.nridge, all_vertices_2D, max_wings,
                                  flattened_all_vertices_2D, layer_points, z_min, z_max, nelem, c.long_connector_ratio,
                                  fiber_ends, ridge_widths=widths)
        nbeamElem, npatch, nconnector_total = lat.nbeamElem, lat.npatch, lat.nconnector_t + lat.nconnector_l
        rplat.write_connector_file(geoName, lat, path=path)
        [mass, bulk_volume, bulk_density, porosity] = \
            rpgen.ModelInfo('square', boundary_points, z_min, z_max, c.skeleton_density, lat.Meshdata, polygon=boundary_points)
        rplat.write_iga_file(geoName, lat, path=path)
        rplat.write_vtu_files(geoName, lat, layer_points, all_ridges, ref.nvertex, path=path, vtk_format=c.vtk_format)
        fibers = {'fiber_ends': c.fiber_ends, 'nelem': nelem}
        if c.fiber_ends == 'random':
            fibers.update({'fiber_length': c.fiber_length, 'fiber_length_cov': c.fiber_length_cov, 'seed': c.seed})
        fibers.update(rplat.fiber_statistics(fiber_ends, nelem, h))
    else:
        widths = np.empty((c.nsegments, nconnector_t_per_beam, ref.nridge))
        for seg in range(c.nsegments):
            for j in range(nconnector_t_per_beam):
                widths[seg, j] = widths_in_layer(seg*nctrlpt_per_beam + j*NURBS_DEGREE)
        connector_t_widths = np.concatenate([widths[:, 0].reshape(-1), widths[:, 1:-1].reshape(-1), widths[:, -1].reshape(-1)])
        [IGAvertices, vertex_connectivity, beam_connectivity_original, nbeam_total,
         beam_connectivity, nbeamElem, nctrlpt_per_beam, nlayers, segment_length, connector_t_connectivity,
         connector_t_bot_connectivity, connector_t_top_connectivity,
         connector_t_reg_connectivity, connector_l_connectivity, nconnector_t_per_beam,
         nconnector_t_per_grain, nconnector_t, nconnector_l, nconnector_total,
         theta, z_coord, connector_l_vertex_dict] = \
            rpgen.GenerateBeamElement(NURBS_DEGREE, c.nsegments, 0.0, 0.0, z_min, z_max, c.long_connector_ratio,
                                      npt_per_layer, ref.vertices, ref.nvertex, ref.ridges, ref.nridge, (0.0, 0.0),
                                      all_vertices_2D, max_wings, flattened_all_vertices_2D, all_ridges,
                                      layer_vertices=layer_vertices)
        [ConnMeshData, conn_l_tangents] = \
            rpgen.ConnectorMeshFile(geoName, IGAvertices, connector_t_bot_connectivity, connector_t_reg_connectivity,
                                    connector_t_top_connectivity, connector_l_connectivity, all_vertices_2D, max_wings,
                                    flattened_all_vertices_2D, c.nsegments, segment_length, nctrlpt_per_beam, theta,
                                    ref.nridge, connector_l_vertex_dict, path=path, layer_vertices=layer_vertices,
                                    connector_t_widths=connector_t_widths)
        [mass, bulk_volume, bulk_density, porosity] = \
            rpgen.ModelInfo('square', boundary_points, z_min, z_max, c.skeleton_density, ConnMeshData, polygon=boundary_points)
        knotVec = rpgen.BezierExtraction(NURBS_DEGREE, nbeam_total)
        npatch = beam_connectivity_original.shape[0]
        rpgen.BezierBeamFile(geoName, NURBS_DEGREE, nctrlpt_per_beam, nconnector_t_per_beam, npatch, knotVec, path=path)
        tangents = connector_tangents(IGAvertices, np.vstack((connector_t_bot_connectivity, connector_t_reg_connectivity,
                                                               connector_t_top_connectivity)), npt_per_layer, nctrlpt_per_beam)
        rpgen.VisualizationFiles(geoName, NURBS_DEGREE, nlayers, npt_per_layer_vtk, all_pts_2D, segment_length, theta,
                                 z_coord, c.nsegments, ref.nridge, ref.ridges, (0.0, 0.0), all_ridges, ref.nvertex,
                                 nconnector_t, nconnector_l, nctrlpt_per_beam, ConnMeshData, conn_l_tangents,
                                 all_vertices_2D, max_wings, flattened_all_vertices_2D, path=path, vtk_format=c.vtk_format,
                                 layer_points=layer_points, connector_t_tangents=tangents)
    meshTime = time.time()
    print('{:d} beam elements of {:d} fibers and {:d} connectors generated in {:.3f} seconds'
          .format(nbeamElem, npatch, nconnector_total, meshTime - trackTime))

    _preview_figure(Path(output_dir, geoName + '.png'), slices[k_ref][1], ref, tracked_vertices, z_coord, inclination,
                    voxel, origin, k_ref, field)
    if c.save_tracks:
        end_pairs = np.array([(g, int(k)) for g, ks in enumerate(fiber_ends) for k in ks], dtype=int).reshape(-1, 2) \
            if fiber_ends is not None else np.zeros((0, 2), dtype=int)
        np.savez(Path(output_dir, geoName + '-grains.npz'), layer_vertices=tracked_vertices,
                 control_points=layer_vertices, z_coord=z_coord, status=layer_status, vertices=ref.vertices,
                 ridges=ref.ridges, wall_thickness=ref.thickness, kind=ref.kind, fiber_ends=end_pairs,
                 z_grid=z_min + np.arange(nelem + 1)*h)
    if c.stlFlag:
        rpgen.StlModelFile(geoName, path=path)
    if c.inpFlag:
        if fibers_on:
            rpgen.AbaqusFile(geoName, NURBS_DEGREE, lat.npatch, lat.nodes, lat.beam_connectivity, lat.connector_t_bot,
                             lat.connector_t_reg, lat.connector_t_top, lat.connector_l, nbeamElem, lat.nconnector_t,
                             lat.nconnector_l, 0, 0, segment_length, c.long_connector_ratio, x_max, x_min, y_max, y_min,
                             lat.nodes[:, 2], 'square', max(x_max - x_min, y_max - y_min), median_wall_thickness,
                             median_wall_thickness, c.merge_operation, merge_tol, 'off', [], path=path,
                             connector_t_ridge=lat.connector_t_ridge)
        else:
            rpgen.AbaqusFile(geoName, NURBS_DEGREE, npatch, IGAvertices, beam_connectivity,
                             connector_t_bot_connectivity, connector_t_reg_connectivity,
                             connector_t_top_connectivity, connector_l_connectivity, nbeamElem,
                             nconnector_t, nconnector_l, 0, 0, segment_length, c.long_connector_ratio,
                             x_max, x_min, y_max, y_min, z_coord, 'square', max(x_max - x_min, y_max - y_min),
                             median_wall_thickness, median_wall_thickness, c.merge_operation, merge_tol,
                             'off', [], path=path)
    fileTime = time.time()
    print('Files generated in {:.3f} seconds'.format(fileTime - meshTime))

    resolved = c.replace(threshold=threshold if threshold is not None else c.threshold, longitudinal_axis=str(axis),
                         merge_tol=merge_tol, direction_scale=direction_scale, max_jump=max_jump,
                         min_cell_area=min_lumen_px*pixel_area, min_wall_area=min_wall_px*pixel_area)
    result = ReconstructionResult(
        config=resolved, output_dir=output_dir, files=[], shape=(nx, ny, nz), voxel_size=tuple(voxel.tolist()),
        origin=tuple(origin.tolist()), threshold=threshold, longitudinal_axis=axis, nParticles=nParticles,
        nbeamElem=nbeamElem, mass=mass, bulk_volume=bulk_volume, bulk_density=bulk_density, porosity=porosity,
        image_porosity=image_porosity,
        grain_angle=field.grain_angle_deg if field is not None else 0.0,
        grain_azimuth=field.azimuth_deg if field is not None else 0.0,
        grain_dispersion=field.dispersion_deg if field is not None else 0.0,
        grain_inclination_mean=float(inclination.mean()) if len(inclination) else 0.0,
        grain_inclination_max=float(inclination.max()) if len(inclination) else 0.0,
        tracking_counts=tracking_counts, median_cell_diameter=median_cell_diameter,
        median_wall_thickness=median_wall_thickness, layer_vertices=tracked_vertices, z_coord=z_coord,
        vertices=ref.vertices, ridges=ref.ridges, fiber_ends=fiber_ends, nfibers=npatch)
    _write_log(result, k_ref, segment_length, nctrlpt_per_beam, nconnector_total, nlabels,
               [startTime, loadTime, labelTime, networkTime, trackTime, meshTime, fileTime], fibers)
    result.files = sorted(output_dir.iterdir())
    return result


def _reference_slice(reference_slice, nz, z_min, vz):
    """Slice index of the reference cross-section: 'middle', 'bottom', 'top', an index, or a height z (mm, with a decimal point)"""
    text = str(reference_slice).strip().lower()
    if text == 'middle':
        return nz//2
    if text == 'bottom':
        return 0
    if text == 'top':
        return nz - 1
    try:
        value = float(text)
    except ValueError:
        raise ValueError("reference_slice must be 'middle', 'bottom', 'top', a slice index or a height, not {!r}"
                         .format(reference_slice)) from None
    if text.lstrip('-').isdigit():
        k = int(value) if value >= 0 else nz + int(value)
    else:
        k = int(math.floor((value - z_min)/vz))
    if not 0 <= k < nz:
        raise ValueError("reference_slice {!r} is outside the volume ({:d} slices)".format(reference_slice, nz))
    return k


def _widths_in_slice(ref, net):
    """Wall thickness of every ridge of the reference network measured in the cross-section `net`
    (the wall between the same two cells); the reference thickness where it is not found"""
    widths = ref.thickness.copy()
    for r, pair in enumerate(ref.pairs):
        if pair[0] < 0:
            continue
        rows = net.pair_index.get(pair)
        if rows:
            widths[r] = float(np.mean(net.thickness[rows]))
    return widths


def _preview_figure(file, walls2d, ref, layer_vertices, z_coord, inclination, voxel, origin, k_ref, field):
    """geoName.png: the reference cross-section with its cell corners and walls, and the top view of the grains"""
    import matplotlib.pyplot as plt
    from matplotlib.collections import LineCollection
    x0, x1, y0, y1 = ref.box
    fig, axes = plt.subplots(1, 2, figsize=(13, 6))
    ax = axes[0]
    ax.imshow(walls2d.T, origin='lower', cmap='Greys', vmin=0, vmax=1.6, extent=[x0, x1, y0, y1], interpolation='nearest')
    segments = ref.vertices[ref.ridges]
    lines = LineCollection(segments, array=ref.thickness, cmap='viridis', linewidths=1.2)
    ax.add_collection(lines)
    fig.colorbar(lines, ax=ax, fraction=0.046, pad=0.04, label='wall thickness (mm)')
    inside = ref.kind == 0
    ax.plot(ref.vertices[inside, 0], ref.vertices[inside, 1], 'r.', markersize=3)
    ax.plot(ref.vertices[~inside, 0], ref.vertices[~inside, 1], 'b.', markersize=3)
    ax.plot([x0, x1, x1, x0, x0], [y0, y0, y1, y1, y0], 'k-', linewidth=2)
    ax.set_xlim(x0 - 0.05*(x1 - x0), x1 + 0.05*(x1 - x0))
    ax.set_ylim(y0 - 0.05*(y1 - y0), y1 + 0.05*(y1 - y0))
    ax.set_aspect('equal', adjustable='box')
    ax.set_title('reference cross-section (slice {:d}, z = {:.4g} mm)\n{:d} cell corners, {:d} cell walls'
                 .format(k_ref, origin[2] + (k_ref + 0.5)*voxel[2], ref.nvertex, ref.nridge))
    ax.set_xlabel('x (mm)')
    ax.set_ylabel('y (mm)')

    ax = axes[1]
    rows = np.nonzero(ref.kind < 2)[0]
    tracks = np.transpose(layer_vertices[:, rows, :], (1, 0, 2))  # (nvertex, nlayers, 2)
    if len(rows):
        paths = LineCollection(tracks, array=inclination, cmap='plasma', linewidths=1.0)
        ax.add_collection(paths)
        fig.colorbar(paths, ax=ax, fraction=0.046, pad=0.04, label='grain inclination (degrees)')
        ax.plot(tracks[:, 0, 0], tracks[:, 0, 1], 'k.', markersize=2)
    ax.plot([x0, x1, x1, x0, x0], [y0, y0, y1, y1, y0], 'k-', linewidth=2)
    ax.set_xlim(x0 - 0.05*(x1 - x0), x1 + 0.05*(x1 - x0))
    ax.set_ylim(y0 - 0.05*(y1 - y0), y1 + 0.05*(y1 - y0))
    ax.set_aspect('equal', adjustable='box')
    title = 'grains followed from z = {:.4g} to {:.4g} mm (top view)'.format(z_coord[0], z_coord[-1])
    if field is not None:
        title += '\nfiber direction field: grain angle {:.2f} deg, dispersion {:.2f} deg'.format(field.grain_angle_deg,
                                                                                                 field.dispersion_deg)
    ax.set_title(title)
    ax.set_xlabel('x (mm)')
    ax.set_ylabel('y (mm)')
    fig.tight_layout()
    fig.savefig(file, dpi=150)
    plt.close(fig)


def _write_log(result, k_ref, segment_length, nctrlpt_per_beam, nconnector_total, nlabels, times, fibers=None):
    """The report log geoName-report.log of a reconstruction"""
    c = result.config
    geoName = c.geoName
    prefix = str(Path(c.path, geoName, geoName))
    current_time = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    startTime, loadTime, labelTime, networkTime, trackTime, meshTime, fileTime = times
    lines = []
    w = lines.append
    w('################################################################################')
    w('##        Log Created with RingsPy Mesh Generation Tool V' + __version__ + ' (reconstruct)    ##')
    w('##             Files generated at local time: ' + current_time + '             ##')
    w('################################################################################')
    w('')
    w('GEOMETRY:')
    w('geoName:                                 ' + geoName)
    w('')
    w('IMAGE:')
    w('volume:                                  ' + str(c.volume))
    w('shape (x, y, z voxels):                  ' + str(list(result.shape)))
    w('voxel_size:                              ' + str(list(result.voxel_size)))
    w('origin:                                  ' + str(list(result.origin)))
    w('longitudinal_axis:                       ' + str(result.longitudinal_axis))
    w('solid:                                   ' + str(c.solid))
    w('threshold:                               ' + str(result.threshold))
    w('smooth_sigma:                            ' + str(c.smooth_sigma))
    w('crop:                                    ' + str(c.crop))
    w('min_cell_area:                           ' + str(c.min_cell_area))
    w('min_wall_area:                           ' + str(c.min_wall_area))
    w('image porosity:                          ' + str(result.image_porosity))
    w('')
    w('RECONSTRUCTION PARAMETERS:')
    w('reference_slice:                         ' + str(c.reference_slice) + ' (slice ' + str(k_ref) + ')')
    w('tracking:                                ' + str(c.tracking))
    w('link_stride:                             ' + str(c.link_stride))
    w('direction_scale:                         ' + str(c.direction_scale))
    w('max_jump:                                ' + str(c.max_jump))
    w('merge_operation:                         ' + rpgen._flag_text(c.merge_operation))
    if c.merge_operation:
        w('merge_tol:                               ' + str(c.merge_tol))
    w('boundaryFlag:                            ' + rpgen._flag_text(c.boundaryFlag))
    w('wall_thickness:                          ' + ('measured' if c.wall_thickness is None else str(c.wall_thickness)))
    w('')
    w('FIBERS:')
    w('grain angle (mean fiber direction):      ' + str(result.grain_angle) + ' deg')
    w('grain azimuth:                           ' + str(result.grain_azimuth) + ' deg')
    w('fiber dispersion:                        ' + str(result.grain_dispersion) + ' deg')
    w('grain inclination (mean):                ' + str(result.grain_inclination_mean) + ' deg')
    w('grain inclination (max):                 ' + str(result.grain_inclination_max) + ' deg')
    w('median cell diameter:                    ' + str(result.median_cell_diameter))
    w('median wall thickness:                   ' + str(result.median_wall_thickness))
    w('')
    w('CELL CORNERS IN THE LAYERS:')
    for name, count in result.tracking_counts.items():
        w('{:<41s}'.format(name + ':') + str(count))
    w('')
    w('GRAIN EXTRUSION PARAMETERS:')
    w('nsegments:                               ' + str(c.nsegments))
    w('segment_length:                          ' + str(segment_length))
    w('long_connector_ratio:                    ' + str(c.long_connector_ratio))
    w('NURBS_degree:                            ' + str(NURBS_DEGREE))
    w('nctrlpt_per_beam:                        ' + str(nctrlpt_per_beam))
    w('')
    w('FIBERS:')
    if fibers is None:
        w('fiber_ends:                              none (all grains cut into nsegments segments)')
    else:
        for key, value in fibers.items():
            w('{:<41s}'.format(key + ':') + str(value))
    w('')
    w('GENERATION:')
    w('nParticles:                              ' + str(result.nParticles))
    w('nlabels (cells in all cross-sections):   ' + str(nlabels))
    w('nvertex:                                 ' + str(len(result.vertices)))
    w('nridge:                                  ' + str(len(result.ridges)))
    w('nbeamElem:                               ' + str(result.nbeamElem))
    w('nconnector:                              ' + str(nconnector_total))
    w('MODEL PROPERTIES:')
    w('skeleton density:                        ' + str(c.skeleton_density))
    w('total mass:                              ' + str(result.mass))
    w('bulk volume:                             ' + str(result.bulk_volume))
    w('bulk density:                            ' + str(result.bulk_density))
    w('porosity:                                ' + str(result.porosity))
    w('')
    w('PERFORMANCE:')
    w('Loading and Thresholding Time:           ' + str(loadTime - startTime))
    w('Cell Labelling Time:                     ' + str(labelTime - loadTime))
    w('Network Extraction Time:                 ' + str(networkTime - labelTime))
    w('Corner Tracking Time:                    ' + str(trackTime - networkTime))
    w('Grain Extrusion Time:                    ' + str(meshTime - trackTime))
    w('File Writing Time:                       ' + str(fileTime - meshTime))
    w('Total Time:                              ' + str(fileTime - startTime))
    w('')
    w('FILE CREATION:')
    w('Cross-sectional Image File:              ' + prefix + '.png')
    w('NURBS Beam File:                         ' + prefix + 'IGA.txt')
    w('Connector Data File:                     ' + prefix + '-mesh.txt')
    w('Grain-ridge Data File:                   ' + prefix + '-vertex.mesh')
    w('Ridge Data File:                         ' + prefix + '-ridge.mesh')
    w('Paraview Vertices File:                  ' + prefix + '_vertices.vtu')
    w('Paraview Beams File:                     ' + prefix + '_beams.vtu')
    w('Paraview Connectors File:                ' + prefix + '_conns.vtu')
    w('Paraview Connectors (Volume) File:       ' + prefix + '_conns_vol.vtu')
    if c.save_tracks:
        w('Grain Tracks File:                       ' + prefix + '-grains.npz')
    if c.stlFlag:
        w('3D Model File:                           ' + prefix + '.stl')
    if c.inpFlag:
        w('Abaqus Input File:                       ' + prefix + '.inp')
    w('Log File:                                ' + prefix + '-report.log')
    with open(Path(c.path, geoName, geoName + '-report.log'), 'w') as logfile:
        logfile.write('\n'.join(lines) + '\n')
