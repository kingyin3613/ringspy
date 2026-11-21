# -*- coding: utf-8 -*-
"""
Figures of the 2D cell mesh of a Result: the cell walls coloured by growth ring, by earlywood and
latewood, by wall thickness or by cell type, with the ring boundaries and the clipping outline.

All RingsPy figures share one palette (also the README figures and animation made in contents/): wood
tones for the earlywood (sand) and the latewood (walnut), a teal accent for the vessels and a violet one
for the rays, warm greys for the heart, the fibers and the lines around them, and four muted hues for
the growth rings. The hues were checked with simulated colour vision deficiencies (protanopia and
deuteranopia): neighbouring growth rings, the two accents, and the accents against the wood tones stay
apart for colour-blind readers too.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import patheffects
from matplotlib.collections import LineCollection
from matplotlib.colors import LinearSegmentedColormap, to_rgba
from matplotlib.legend_handler import HandlerTuple
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator

COLORINGS = ['ring', 'dense', 'thickness', 'type', None]
COLORS = {'sparse': '#c9ab87',     # earlywood: the sparse part of a ring (sand)
          'dense': '#6a4225',      # latewood: the dense part of a ring, and the heart (walnut)
          'heart': '#9a948d',      # the heart with color='ring' (warm grey)
          'fiber': '#b9b3ac',      # fiber walls with color='type' (light warm grey)
          'vessel': '#00747a',     # vessel walls (deep teal)
          'ray': '#9a7dbb',        # ray walls (soft violet)
          'boundary': '#2b2521',   # the clipping outline (ink)
          'rings': '#b5aca3',      # the growth ring boundaries
          'walls': '#2b2521',      # all walls with color=None
          'text': '#534b45', 'muted': '#7a736b', 'frame': '#d3cdc6'}
# color='ring': one hue per growth ring, repeated after the fourth (copper, sky, olive, lavender); the
# sparse part of a ring takes the lighter step, its dense part the darker one
RING_COLORS = [('#c4854f', '#7d460b'), ('#539dcf', '#085b87'), ('#8d9f51', '#4e5d0e'), ('#9f86cc', '#5e4785')]
# color='thickness': light sand to dark walnut (one hue, lightness falling monotonically); a second
# quantity in the same figure (the grain inclination of a reconstruction) takes the slate ramp
THICKNESS_CMAP = LinearSegmentedColormap.from_list(
    'ringspy_wall', ['#cdb797', '#bb9b76', '#a57d57', '#875f3f', '#68432c', '#482b1c'])
ANGLE_CMAP = LinearSegmentedColormap.from_list(
    'ringspy_angle', ['#a6c6db', '#7ea7c7', '#5988b4', '#3f689b', '#2f497b', '#222d56'])
MAX_LINEWIDTH = 1.4         # points, for the thickest wall (linewidth='thickness')
MIN_LINEWIDTH = 0.35
PRECRACK_LINEWIDTH = 2.0    # points: a pre-crack is a cut in the ink of the outline, with a white edge
MARGIN = 0.03               # of the box size, around the clipping outline
# rcParams of the figures RingsPy draws itself (the preview PNG, the reconstruction figure, the README
# figures): quiet frame and ticks, so that the cell walls carry the figure; the first installed of FONTS
FONTS = ['Lato', 'Roboto', 'Segoe UI', 'Helvetica Neue', 'Arial', 'DejaVu Sans']
STYLE = {'font.family': 'sans-serif', 'font.sans-serif': FONTS,
         'axes.edgecolor': COLORS['frame'], 'axes.linewidth': 0.7, 'axes.labelcolor': COLORS['text'],
         'axes.titlecolor': COLORS['boundary'], 'axes.titlesize': 'medium', 'axes.titlepad': 8.0,
         'xtick.color': COLORS['frame'], 'ytick.color': COLORS['frame'], 'xtick.labelcolor': COLORS['muted'],
         'ytick.labelcolor': COLORS['muted'], 'xtick.major.width': 0.7, 'ytick.major.width': 0.7,
         'xtick.major.size': 3.0, 'ytick.major.size': 3.0, 'legend.labelcolor': COLORS['text'],
         'legend.frameon': False, 'legend.handlelength': 1.6, 'figure.facecolor': 'white',
         'savefig.facecolor': 'white'}


def plot_mesh(result, ax=None, color='ring', linewidth='thickness', types=True, rings=True, boundary=True,
              sites=False, precracks=None, legend=True, title=True):
    """Draw the clipped 2D cell mesh of a Result (of generate(), preview() or a reconstruction) in a
    Matplotlib axes, and return the axes.

    color :: what the colour of the walls shows: 'ring' (one hue per growth ring, the heart grey, the
        dense part of a ring darker than its sparse part), 'dense' (earlywood light, latewood dark, as
        in a cross-section of wood), 'thickness' (the wall thickness, with a colour bar), 'type' (fiber,
        vessel and ray walls of a hardwood), or None (all dark).
    linewidth :: 'thickness' (line widths proportional to the wall thickness, the thickest MAX_LINEWIDTH
        points wide), or a number of points for all walls.
    types :: also draw the vessel and ray walls of a hardwood in their colours with color 'ring',
        'dense' or 'thickness' (no effect for a softwood).
    rings :: draw the growth ring boundaries (dashed, inside the clipping outline).
    boundary :: draw the clipping outline.
    sites :: also mark the cell sites (small dots, coloured by cell type).
    precracks :: line segments [[x1, y1, x2, y2], ...] of pre-cracks to draw across the walls (the figure
        geoName.png of a run with precrackFlag has those of the run).
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

    handles, labels = [], []
    legend_title = None
    mappable = None
    if color == 'ring':
        colors = ring_colors(result.ridge_ring, result.ridge_dense, int(max(result.radii.size//2, 1)))
        nrings = int(result.ridge_ring.max()) if len(ridges) else 0
        handles.append(_swatch(COLORS['heart']))
        labels.append('heart')
        for k in range(1, nrings + 1):
            handles.append((_swatch(ring_hue(k, nrings)), _swatch(ring_hue(k, nrings, dense=True))))
            labels.append('ring {:d}'.format(k))
        legend_title = 'sparse | dense'
    elif color == 'dense':
        colors = np.where(result.ridge_dense, COLORS['dense'], COLORS['sparse'])
        handles += [_swatch(COLORS['sparse']), _swatch(COLORS['dense'])]
        labels += ['sparse (earlywood)', 'dense (latewood, heart)']
    elif color == 'thickness':
        norm = plt.Normalize(thickness.min(), thickness.max())
        colors = THICKNESS_CMAP(norm(thickness))
        mappable = plt.cm.ScalarMappable(norm=norm, cmap=THICKNESS_CMAP)
    elif color == 'type':
        colors = np.where(result.ridge_vessel, COLORS['vessel'], np.where(result.ridge_ray, COLORS['ray'], COLORS['fiber']))
        handles.append(_swatch(COLORS['fiber']))
        labels.append('fiber walls')
    else:
        colors = np.full(len(ridges), COLORS['walls'], dtype=object)
    colors = np.asarray(colors, dtype=object) if colors.dtype.kind in 'OU' else np.asarray(colors)
    if hardwood and (types or color == 'type'):
        colors = np.array([_to_rgba(c) for c in colors])
        colors[np.asarray(result.ridge_vessel, dtype=bool)] = _to_rgba(COLORS['vessel'])
        colors[np.asarray(result.ridge_ray, dtype=bool)] = _to_rgba(COLORS['ray'])
        handles += [_swatch(COLORS['vessel']), _swatch(COLORS['ray'])]
        labels += ['vessel walls', 'ray walls']
    elif color == 'type':
        handles, labels = [], []

    walls = LineCollection(segments, colors=list(colors) if colors.dtype == object else colors, linewidths=widths,
                           capstyle='round', joinstyle='round', zorder=2)
    ax.add_collection(walls)

    lo, hi = vertices.min(axis=0), vertices.max(axis=0)
    outline = getattr(result, 'boundary_points', None)
    clip = None
    if outline is not None and len(outline):
        outline = np.asarray(outline, dtype=float)
        lo, hi = np.minimum(lo, outline.min(axis=0)), np.maximum(hi, outline.max(axis=0))
        clip = plt.Polygon(outline, closed=True, transform=ax.transData)   # the ring boundaries stay inside
        if boundary:
            ax.add_patch(plt.Polygon(outline, closed=True, fill=False, edgecolor=COLORS['boundary'], linewidth=1.0,
                                     zorder=3))
    margin = MARGIN*max(hi - lo)
    ax.set_xlim(lo[0] - margin, hi[0] + margin)
    ax.set_ylim(lo[1] - margin, hi[1] + margin)

    shape = getattr(result, 'ring_shape', None)
    if rings and shape is not None and result.radii is not None:
        for r in np.asarray(result.radii, dtype=float)[1:]:
            ring = shape.boundary(r)
            line, = ax.plot(ring[:, 0], ring[:, 1], color=COLORS['rings'], linestyle=(0, (4, 3)), linewidth=0.7,
                            zorder=1)
            if clip is not None:
                line.set_clip_path(clip)
    if sites and result.sites is not None:
        pts = np.asarray(result.sites, dtype=float)
        site_type = np.asarray(getattr(result, 'site_type', np.zeros(len(pts))), dtype=int)
        site_colors = np.array([COLORS['muted'], COLORS['vessel'], COLORS['ray']])[np.clip(site_type, 0, 2)]
        ax.scatter(pts[:, 0], pts[:, 1], s=2, c=site_colors, linewidths=0, zorder=4)
    if precracks is not None and len(precracks):
        cuts = np.asarray(precracks, dtype=float).reshape(-1, 2, 2)
        edge = [patheffects.withStroke(linewidth=PRECRACK_LINEWIDTH + 2.4, foreground='white')]
        for cut in cuts:
            ax.plot(cut[:, 0], cut[:, 1], color=COLORS['boundary'], linewidth=PRECRACK_LINEWIDTH,
                    solid_capstyle='round', path_effects=edge, zorder=5)
        handles.append(Line2D([], [], color=COLORS['boundary'], linewidth=PRECRACK_LINEWIDTH, solid_capstyle='round'))
        labels.append('pre-crack' if len(cuts) == 1 else 'pre-cracks')

    ax.set_aspect('equal', adjustable='box')
    ax.set_xlabel('x (mm)')
    ax.set_ylabel('y (mm)')
    if mappable is not None:
        # a slim bar beside the drawing, as tall as the drawing (and as wide for a wide one as for a square one)
        x0, x1 = ax.get_xlim()
        y0, y1 = ax.get_ylim()
        cax = ax.inset_axes([1.03, 0.0, 0.035*min(1.0, (y1 - y0)/max(x1 - x0, 1e-12)), 1.0])
        bar = ax.figure.colorbar(mappable, cax=cax, label='wall thickness (mm)')
        bar.outline.set_visible(False)
        bar.ax.tick_params(length=0)
        bar.ax.yaxis.set_major_locator(MaxNLocator(nbins=4, steps=[1, 2, 5, 10]))   # a few round values
    if legend and handles:
        key = ax.legend(handles, labels, loc='upper left', bbox_to_anchor=(1.01, 1.0), fontsize='small', frameon=False,
                        title=legend_title, title_fontsize='small',
                        handler_map={tuple: HandlerTuple(ndivide=None, pad=0.0)})
        key.get_title().set_color(COLORS['muted'])
    if title:
        name = result.config.geoName if getattr(result, 'config', None) is not None else ''
        counts = '{:d} cells placed, {:d} vertices and {:d} walls in the box'.format(
            int(result.nParticles), len(vertices), len(ridges))
        if getattr(result, 'cell_stats', None):  # the cells (partly) in the box have an area
            counts = '{:d} cells, {:d} vertices and {:d} walls in the box'.format(
                int(np.isfinite(np.asarray(result.cell_stats['area'], dtype=float)).sum()), len(vertices), len(ridges))
        ax.set_title((name + ': ' if name else '') + counts, fontsize='medium')
    return ax


def ring_hue(ring, nrings, dense=False):
    """The colour of the sparse part (or, dense=True, of the dense part) of growth ring `ring`
    (1 ... nrings) with color='ring' (RING_COLORS, repeated after the fourth ring)"""
    sparse_color, dense_color = RING_COLORS[(int(ring) - 1) % len(RING_COLORS)]
    return _to_rgba(dense_color if dense else sparse_color)


def ring_colors(ring, dense, nrings):
    """RGBA colours of walls from their growth ring (0: the heart) and dense flags"""
    ring = np.asarray(ring, dtype=int)
    dense = np.asarray(dense, dtype=bool)
    nrings = max(int(ring.max()) if len(ring) else nrings, 1)
    colors = np.empty((len(ring), 4))
    colors[:] = _to_rgba(COLORS['heart'])
    for k in range(1, nrings + 1):
        colors[(ring == k) & ~dense] = ring_hue(k, nrings)
        colors[(ring == k) & dense] = ring_hue(k, nrings, dense=True)
    return colors


def _swatch(color):
    """A legend handle: a short thick line in `color`"""
    return Line2D([], [], color=color, linewidth=3.0, solid_capstyle='butt')


def _to_rgba(color):
    return to_rgba(color)
