# -*- coding: utf-8 -*-
"""
Figures of the 2D cell mesh of a Result: the cell walls coloured by growth ring, by earlywood and
latewood, by wall thickness or by cell type, with the ring boundaries and the clipping outline
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.lines import Line2D

COLORINGS = ['ring', 'dense', 'thickness', 'type', None]
COLORS = {'heart': '#7f7f7f', 'sparse': '#deb887', 'dense': '#6b3e1e', 'fiber': '#9e9e9e',
          'vessel': '#d62728', 'ray': '#ff7f0e', 'boundary': 'k', 'rings': '#555555'}
RING_CMAP = 'viridis'       # one hue per growth ring (color='ring'), the dense part darker
THICKNESS_CMAP = 'magma_r'  # color='thickness'
MAX_LINEWIDTH = 2.0         # points, for the thickest wall (linewidth='thickness')
MIN_LINEWIDTH = 0.4
MARGIN = 0.03               # of the box size, around the clipping outline


def plot_mesh(result, ax=None, color='ring', linewidth='thickness', types=True, rings=True, boundary=True,
              sites=False, legend=True, title=True):
    """Draw the clipped 2D cell mesh of a Result (of generate(), preview() or a reconstruction) in a
    Matplotlib axes, and return the axes.

    color :: what the colour of the walls shows: 'ring' (one hue per growth ring, the heart grey, the
        dense part of a ring darker than its sparse part), 'dense' (earlywood light, latewood dark, as
        in a cross-section of wood), 'thickness' (the wall thickness, with a colour bar), 'type' (fiber,
        vessel and ray walls of a hardwood), or None (all black).
    linewidth :: 'thickness' (line widths proportional to the wall thickness, the thickest MAX_LINEWIDTH
        points wide), or a number of points for all walls.
    types :: also draw the vessel and ray walls of a hardwood in their colours with color 'ring',
        'dense' or 'thickness' (no effect for a softwood).
    rings :: draw the growth ring boundaries (dashed).
    boundary :: draw the clipping outline.
    sites :: also mark the cell sites (small dots, coloured by cell type).
    legend, title :: add a legend of the colours, and a title with the numbers of cells and walls.
    """
    if color not in COLORINGS:
        raise ValueError("color must be one of {}, not {!r}".format(COLORINGS, color))
    if ax is None:
        _, ax = plt.subplots(figsize=(8, 8))
    vertices = np.asarray(result.vertices, dtype=float)
    ridges = np.asarray(result.ridges, dtype=int)
    thickness = np.asarray(result.vertex_thickness, dtype=float)[ridges[:, 1]]  # the thickness of a wall
    segments = vertices[ridges]
    hardwood = bool(np.any(result.ridge_vessel) or np.any(result.ridge_ray))

    if linewidth == 'thickness':
        widths = np.maximum(MAX_LINEWIDTH*thickness/max(thickness.max(), 1e-12), MIN_LINEWIDTH)
    else:
        widths = np.full(len(ridges), float(linewidth))

    handles = []
    mappable = None
    if color == 'ring':
        colors = ring_colors(result.ridge_ring, result.ridge_dense, int(max(result.radii.size//2, 1)))
        nrings = int(result.ridge_ring.max()) if len(ridges) else 0
        handles.append(Line2D([], [], color=COLORS['heart'], lw=2, label='heart'))
        for k in range(1, nrings + 1):
            handles.append(Line2D([], [], color=ring_hue(k, nrings), lw=2, label='ring {:d}'.format(k)))
        handles.append(Line2D([], [], color=_darker(ring_hue(max(nrings, 1), nrings)), lw=2, label='dense part (darker)'))
    elif color == 'dense':
        colors = np.where(result.ridge_dense, COLORS['dense'], COLORS['sparse'])
        handles.append(Line2D([], [], color=COLORS['sparse'], lw=2, label='sparse (earlywood)'))
        handles.append(Line2D([], [], color=COLORS['dense'], lw=2, label='dense (latewood, heart)'))
    elif color == 'thickness':
        cmap = plt.get_cmap(THICKNESS_CMAP)
        norm = plt.Normalize(thickness.min(), thickness.max())
        colors = cmap(norm(thickness))
        mappable = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    elif color == 'type':
        colors = np.where(result.ridge_vessel, COLORS['vessel'], np.where(result.ridge_ray, COLORS['ray'], COLORS['fiber']))
        handles.append(Line2D([], [], color=COLORS['fiber'], lw=2, label='fibers'))
    else:
        colors = np.full(len(ridges), 'k', dtype=object)
    colors = np.asarray(colors, dtype=object) if colors.dtype.kind in 'OU' else np.asarray(colors)
    if hardwood and (types or color == 'type'):
        colors = np.array([_to_rgba(c) for c in colors])
        colors[np.asarray(result.ridge_vessel, dtype=bool)] = _to_rgba(COLORS['vessel'])
        colors[np.asarray(result.ridge_ray, dtype=bool)] = _to_rgba(COLORS['ray'])
        handles.append(Line2D([], [], color=COLORS['vessel'], lw=2, label='vessel walls'))
        handles.append(Line2D([], [], color=COLORS['ray'], lw=2, label='ray walls'))
    elif color == 'type':
        handles = []

    walls = LineCollection(segments, colors=list(colors) if colors.dtype == object else colors, linewidths=widths,
                           capstyle='round', joinstyle='round', zorder=2)
    ax.add_collection(walls)

    lo, hi = vertices.min(axis=0), vertices.max(axis=0)
    outline = getattr(result, 'boundary_points', None)
    if outline is not None and len(outline):
        outline = np.asarray(outline, dtype=float)
        lo, hi = np.minimum(lo, outline.min(axis=0)), np.maximum(hi, outline.max(axis=0))
        if boundary:
            ax.add_patch(plt.Polygon(outline, closed=True, fill=False, edgecolor=COLORS['boundary'], linewidth=1.2,
                                     zorder=3))
    margin = MARGIN*max(hi - lo)
    ax.set_xlim(lo[0] - margin, hi[0] + margin)
    ax.set_ylim(lo[1] - margin, hi[1] + margin)

    shape = getattr(result, 'ring_shape', None)
    if rings and shape is not None and result.radii is not None:
        for r in np.asarray(result.radii, dtype=float)[1:]:
            ring = shape.boundary(r)
            ax.plot(ring[:, 0], ring[:, 1], color=COLORS['rings'], linestyle='--', linewidth=0.6, zorder=1)
    if sites and result.sites is not None:
        pts = np.asarray(result.sites, dtype=float)
        site_type = np.asarray(getattr(result, 'site_type', np.zeros(len(pts))), dtype=int)
        site_colors = np.array([COLORS['fiber'], COLORS['vessel'], COLORS['ray']])[np.clip(site_type, 0, 2)]
        ax.scatter(pts[:, 0], pts[:, 1], s=3, c=site_colors, linewidths=0, zorder=4)

    ax.set_aspect('equal', adjustable='box')
    ax.set_xlabel('x (mm)')
    ax.set_ylabel('y (mm)')
    if mappable is not None:
        plt.colorbar(mappable, ax=ax, fraction=0.046, pad=0.04, label='wall thickness (mm)')
    if legend and handles:
        ax.legend(handles=handles, loc='upper left', bbox_to_anchor=(1.01, 1.0), fontsize='small', frameon=False)
    if title:
        name = result.config.geoName if getattr(result, 'config', None) is not None else ''
        counts = '{:d} cells placed, {:d} vertices and {:d} walls in the box'.format(
            int(result.nParticles), len(vertices), len(ridges))
        if getattr(result, 'cell_stats', None):  # the cells (partly) in the box have an area
            counts = '{:d} cells, {:d} vertices and {:d} walls in the box'.format(
                int(np.isfinite(np.asarray(result.cell_stats['area'], dtype=float)).sum()), len(vertices), len(ridges))
        ax.set_title((name + ': ' if name else '') + counts, fontsize='medium')
    return ax


def ring_hue(ring, nrings):
    """The colour of a growth ring (1 ... nrings) with color='ring'"""
    cmap = plt.get_cmap(RING_CMAP)
    return cmap(0.15 + 0.7*(ring - 1)/max(nrings - 1, 1))


def ring_colors(ring, dense, nrings):
    """RGBA colours of walls from their growth ring (0: the heart) and dense flags"""
    ring = np.asarray(ring, dtype=int)
    dense = np.asarray(dense, dtype=bool)
    nrings = max(int(ring.max()) if len(ring) else nrings, 1)
    colors = np.empty((len(ring), 4))
    colors[:] = _to_rgba(COLORS['heart'])
    for k in range(1, nrings + 1):
        hue = np.array(ring_hue(k, nrings))
        colors[(ring == k) & ~dense] = hue
        colors[(ring == k) & dense] = _darker(hue)
    return colors


def _darker(rgba, factor=0.55):
    rgba = np.array(_to_rgba(rgba), dtype=float)
    rgba[:3] *= factor
    return tuple(rgba)


def _to_rgba(color):
    from matplotlib.colors import to_rgba
    return to_rgba(color)
