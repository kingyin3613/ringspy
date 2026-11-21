# -*- coding: utf-8 -*-
"""
Rebuild the README figures in contents/ (the animation has its own script, make_hero_gif.py):

    Model2DPreview.png       previews of the 2D cross-sections of seven examples, drawn with result.plot()
    wood_cube_preview.png    the figure of `ringspy run wood_cube.toml --preview`
    ModelVisualization.png   three off-screen VTK renders of the ParaView files: the components of wood_cube
                             coloured by Connector_width, the connector volumes of the hardwood example with
                             the Vessel and Ray cell data, and the tapered, spiral-grained volumes of
                             wood_spiral coloured by Ring
    wood_ct_preview.png      the figure of a reconstruction from a synthetic CT scan (examples/reconstruction)
    lloyd.gif                Lloyd's relaxation of the cells of the heart (radial_growth_rule = "binary_lloyd")
    growth_rule_binary.png   the binary radial growth rule

    python contents/make_screenshots.py --work <folder> [--only previews views ct lloyd growth]

The meshes, the synthetic scan and the reconstruction are made in <folder> (some minutes; the scan and its
reconstruction take about 2 GB of memory). All figures use the RingsPy palette of ringspy/plot.py: the 2D
figures on white, the 3D renders on the charcoal of the animation (the palette stepped for a dark
background, make_hero_gif.py). Needs vtk and Pillow, an OpenGL context for the off-screen renders (on a
machine without a GPU: the vtk-osmesa wheel from https://wheels.vtk.org), and hexalattice for the
honeycomb example. The text uses the first installed font of ringspy.plot.FONTS.
"""

import argparse
import shutil
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
EXAMPLES = HERE.parent / 'examples'
OUT = HERE
PREVIEWS = [('wood_cube', 'ring'), ('wood_cube_lloyd', 'ring'), ('hardwood_ring_porous', 'type'),
            ('wood_eccentric', 'dense'), ('wood_spiral', 'ring'), ('honeycomb_hex', 'ring'), ('dogbone_svg', 'thickness')]
MESHES = ['wood_cube', 'hardwood_ring_porous', 'wood_spiral']
STEPS = ['previews', 'views', 'ct', 'lloyd', 'growth']

# Connector_width on the charcoal background of the 3D renders (as in the animation): sand (thin) to copper
# (thick), lightness falling monotonically; the grains are drawn in a dark warm grey so that the walls
# carry the colour. The growth rings take the colours of color='ring' (ringspy.plot.RING_COLORS).
WIDTH_DARK = ['#dac09d', '#d1ad84', '#ca996b', '#c08252', '#b36c3e', '#9d5633']
GRAIN_DARK = (0x7A, 0x74, 0x6C)


def _hex(rgb):
    return '#{:02x}{:02x}{:02x}'.format(*rgb)


def _pyplot():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    return plt


# ============================================================================ 2D previews
def previews_2d(work):
    """Model2DPreview.png: six square panels and the wide dog-bone panel, each the preview of an example"""
    plt = _pyplot()
    from matplotlib import gridspec
    import ringspy
    from ringspy.plot import COLORS, STYLE
    with plt.rc_context(dict(STYLE, **{'font.size': 12.5})):    # larger text: the README shows it at half size
        fig = plt.figure(figsize=(19, 16.4))
        gs = gridspec.GridSpec(3, 3, figure=fig, height_ratios=[1, 1, 0.72], left=0.05, right=0.925, top=0.9,
                               bottom=0.04, hspace=0.42, wspace=0.72)
        axes = [fig.add_subplot(gs[i // 3, i % 3]) for i in range(6)] + [fig.add_subplot(gs[2, :])]
        for ax, (name, color) in zip(axes, PREVIEWS):
            r = ringspy.preview(EXAMPLES / (name + '.toml'), path=str(work / 'previews'))
            r.plot(ax=ax, color=color, legend=(name != 'honeycomb_hex'))
            ax.set_title('')
            ax.text(0.0, 1.075, name + '.toml', transform=ax.transAxes, ha='left', va='bottom', fontsize=16,
                    fontweight='semibold', color=COLORS['boundary'])
            ax.text(0.0, 1.022, 'color="{}"   {:d} cells, {:d} walls'.format(color, r.summary()['ncells'], r.nridge),
                    transform=ax.transAxes, ha='left', va='bottom', fontsize=12.5, color=COLORS['muted'])
        fig.text(0.05, 0.966, 'Previews of the examples', ha='left', va='bottom', fontsize=22, fontweight='semibold',
                 color=COLORS['boundary'])
        fig.text(0.05, 0.945, 'ringspy run <example>.toml --preview, drawn with result.plot(color=...)', ha='left',
                 va='bottom', fontsize=15, color=COLORS['muted'])
        fig.savefig(OUT / 'Model2DPreview.png', dpi=100)
        plt.close(fig)


def preview_figure(work):
    """wood_cube_preview.png: the figure `ringspy run wood_cube.toml --preview` writes, scaled to 1000 px"""
    _pyplot()
    from PIL import Image
    import ringspy
    r = ringspy.preview(EXAMPLES / 'wood_cube.toml', path=str(work / 'cli'))
    image = Image.open(r.output_dir / 'wood_cube.png').convert('RGB')
    width = 1000
    image.resize((width, round(image.height*width/image.width)), Image.LANCZOS).save(OUT / 'wood_cube_preview.png',
                                                                                    optimize=True)


# ============================================================================ 3D renders
def views_3d(work):
    """ModelVisualization.png: the three VTK renders side by side on the charcoal of the animation, with
    their keys drawn by Matplotlib"""
    import vtk
    from vtkmodules.util.numpy_support import vtk_to_numpy, numpy_to_vtk
    from PIL import Image
    plt = _pyplot()
    from matplotlib.colors import LinearSegmentedColormap, Normalize
    from matplotlib.legend_handler import HandlerTuple
    from matplotlib.lines import Line2D
    import ringspy
    from ringspy.plot import COLORS, RING_COLORS, STYLE
    sys.path.insert(0, str(HERE))
    from make_hero_gif import BG, CAPTION, LIGHT, RAY, VESSEL

    for name in MESHES:
        if not (work / 'meshes' / name / (name + '_conns_vol.vtu')).is_file():
            ringspy.generate(EXAMPLES / (name + '.toml'), path=str(work / 'meshes'), inpFlag=False, stlFlag=False)
    m = str(work / 'meshes') + '/'
    W, H, SS = 900, 760, 2          # panel size; rendered SS times larger and scaled down (antialiasing)
    width_cmap = LinearSegmentedColormap.from_list('ringspy_wall_dark', WIDTH_DARK)

    def read(path):
        r = vtk.vtkXMLUnstructuredGridReader()
        r.SetFileName(path)
        r.Update()
        return r.GetOutput()

    def cell_array(grid, name):
        return vtk_to_numpy(grid.GetCellData().GetArray(name))

    def actor(grid, rgb=None, color=None, line_width=1.0):
        mapper = vtk.vtkDataSetMapper()
        mapper.SetInputData(grid)
        if rgb is not None:
            colors = numpy_to_vtk(np.ascontiguousarray(rgb, dtype=np.uint8), deep=True)
            colors.SetName('rgb')
            grid.GetCellData().AddArray(colors)
            mapper.SetScalarModeToUseCellFieldData()
            mapper.SelectColorArray('rgb')
            mapper.SetColorModeToDirectScalars()
        else:
            mapper.ScalarVisibilityOff()
        a = vtk.vtkActor()
        a.SetMapper(mapper)
        p = a.GetProperty()
        if color is not None:
            p.SetColor(*[v/255 for v in color])
        p.SetLineWidth(line_width*SS)
        p.SetInterpolationToFlat()
        p.SetAmbient(0.28)
        p.SetDiffuse(0.78)
        p.SetSpecular(0.0)
        return a

    def render(actors):
        ren = vtk.vtkRenderer()
        ren.SetBackground(*[v/255 for v in BG])
        ren.AutomaticLightCreationOff()
        for a in actors:
            ren.AddActor(a)
        rw = vtk.vtkRenderWindow()
        rw.SetOffScreenRendering(1)
        rw.SetSize(W*SS, H*SS)
        rw.AddRenderer(ren)
        cam = ren.GetActiveCamera()
        b = np.array(actors[0].GetBounds()).reshape(3, 2)
        c = b.mean(axis=1)
        d = 3.2*np.ptp(b, axis=1).max()
        # seen from 40 degrees above, so that the cross-section on top shows; the specimen a little smaller
        # than the panel and right of its centre, the text of the panel on the left
        cam.SetFocalPoint(*c)
        cam.SetPosition(*(c + d*np.array([0.95, -1.15, 1.25])/1.7))
        cam.SetViewUp(0, 0, 1)
        cam.SetViewAngle(24)
        ren.ResetCamera()
        cam.Zoom(0.86)
        cam.SetWindowCenter(-0.17, 0.03)
        ren.ResetCameraClippingRange()
        # a light from above (the top face brightest), and a softer one from the camera
        top = vtk.vtkLight()
        top.SetLightTypeToSceneLight()
        top.SetPosition(c[0] - 0.4*d, c[1] - 0.6*d, c[2] + 2.0*d)
        top.SetFocalPoint(*c)
        top.SetIntensity(0.75)
        head = vtk.vtkLight()
        head.SetLightTypeToHeadlight()
        head.SetIntensity(0.45)
        ren.AddLight(top)
        ren.AddLight(head)
        rw.Render()
        grab = vtk.vtkWindowToImageFilter()
        grab.SetInput(rw)
        grab.ReadFrontBufferOff()
        grab.Update()
        im = grab.GetOutput()
        w, h, _ = im.GetDimensions()
        px = vtk_to_numpy(im.GetPointData().GetScalars()).reshape(h, w, -1)[::-1, :, :3]
        return Image.fromarray(np.ascontiguousarray(px)).resize((W, H), Image.LANCZOS)

    def width_rgb(width, lo, hi):
        return np.rint(255*width_cmap(Normalize(lo, hi)(width))[:, :3])

    # (a) wood_cube components: beams and connectors, the connectors coloured by Connector_width
    conns = read(m + 'wood_cube/wood_cube_conns.vtu')
    beams = read(m + 'wood_cube/wood_cube_beams.vtu')
    w = cell_array(conns, 'Connector_width')
    rng_a = (float(w.min()), float(w.max()))
    panel_a = render([actor(conns, rgb=width_rgb(w, *rng_a), line_width=1.0),
                      actor(beams, color=GRAIN_DARK, line_width=0.5)])
    # (b) hardwood connector volumes: Connector_width, the vessels and rays in their colours
    vol = read(m + 'hardwood_ring_porous/hardwood_ring_porous_conns_vol.vtu')
    w = cell_array(vol, 'Connector_width')
    rng_b = (float(w.min()), float(w.max()))
    rgb = width_rgb(w, *rng_b)
    rgb[cell_array(vol, 'Vessel') > 0] = VESSEL
    rgb[cell_array(vol, 'Ray') > 0] = RAY
    panel_b = render([actor(vol, rgb=rgb)])
    # (c) wood_spiral volumes, coloured by Ring (the dense part of a ring in the darker step), with the taper
    # and the spiral grain
    vol3 = read(m + 'wood_spiral/wood_spiral_conns_vol.vtu')
    ring = cell_array(vol3, 'Ring').astype(int)
    dense = cell_array(vol3, 'Dense') > 0
    rgb = np.zeros((len(ring), 3))
    rgb[:] = _rgb(COLORS['heart'])
    nrings = int(ring.max())
    ring_pair = lambda k: RING_COLORS[(k - 1) % len(RING_COLORS)]
    for k in range(1, nrings + 1):
        rgb[(ring == k) & ~dense] = _rgb(ring_pair(k)[0])
        rgb[(ring == k) & dense] = _rgb(ring_pair(k)[1])
    panel_c = render([actor(vol3, rgb=rgb)])

    bg, light, caption = _hex(BG), _hex(LIGHT), _hex(CAPTION)
    with plt.rc_context(STYLE):
        width = 18.1                                      # inches at 100 dpi: about as wide as before
        fig = plt.figure(figsize=(width, width*H/(3*W)), dpi=100, facecolor=bg)
        panels = [(panel_a, '(a)  wood_cube', '_beams.vtu + _conns.vtu, coloured by Connector_width'),
                  (panel_b, '(b)  hardwood_ring_porous', '_conns_vol.vtu, Connector_width, Vessel and Ray'),
                  (panel_c, '(c)  wood_spiral', '_conns_vol.vtu, coloured by Ring: taper and spiral grain')]
        for i, (image, title, subtitle) in enumerate(panels):
            ax = fig.add_axes([i/3, 0, 1/3, 1])
            ax.imshow(np.asarray(image), interpolation='antialiased')
            ax.set_axis_off()
            ax.text(0.045, 0.955, title, transform=ax.transAxes, ha='left', va='top', fontsize=17.5, fontweight='semibold',
                    color=light)
            ax.text(0.045, 0.905, subtitle, transform=ax.transAxes, ha='left', va='top', fontsize=13, color=caption)
            if i < 2:                  # Connector_width: a slim horizontal bar at the bottom left
                lo, hi = (rng_a, rng_b)[i]
                cax = ax.inset_axes([0.05, 0.075, 0.33, 0.022])
                bar = fig.colorbar(plt.cm.ScalarMappable(norm=Normalize(lo, hi), cmap=width_cmap), cax=cax,
                                   orientation='horizontal')
                bar.outline.set_visible(False)
                bar.set_ticks([])
                for x, value, align in ((0.05, lo, 'left'), (0.38, hi, 'right')):   # the ends, inside the panel
                    ax.text(x, 0.062, '{:.3g} mm'.format(value), transform=ax.transAxes, ha=align, va='top',
                            fontsize=12, color=caption)
                ax.text(0.05, 0.112, 'Connector_width', transform=ax.transAxes, ha='left', va='bottom',
                        fontsize=13, color=caption)
            if i == 1:
                swatch = lambda color: Line2D([], [], color=_hex(color), linewidth=5, solid_capstyle='butt')
                ax.legend([swatch(VESSEL), swatch(RAY)], ['Vessel', 'Ray'], loc='lower left',
                          bbox_to_anchor=(0.035, 0.16), frameon=False, fontsize=13, labelcolor=caption,
                          handlelength=1.4)
            if i == 2:
                swatch = lambda color: Line2D([], [], color=color, linewidth=6, solid_capstyle='butt')
                handles = [swatch(COLORS['heart'])] + [(swatch(ring_pair(k)[0]), swatch(ring_pair(k)[1]))
                                                       for k in range(1, nrings + 1)]
                labels = ['heart'] + ['ring {:d}'.format(k) for k in range(1, nrings + 1)]
                key = ax.legend(handles, labels, loc='upper left', bbox_to_anchor=(0.035, 0.86), frameon=False,
                                fontsize=13, labelcolor=caption, handlelength=1.8, title='sparse | dense',
                                title_fontsize=12, alignment='left', handler_map={tuple: HandlerTuple(ndivide=None, pad=0)})
                key.get_title().set_color(caption)
        fig.savefig(OUT / 'ModelVisualization.png', dpi=100, facecolor=bg)
        plt.close(fig)


def _rgb(color):
    from matplotlib.colors import to_rgb
    return np.rint(255*np.array(to_rgb(color)))


# ============================================================================ reconstruction
def ct_preview(work):
    """wood_ct_preview.png: the figure of the reconstruction of examples/reconstruction (its synthetic scan
    is made in <work>, as wood_ct.py makes it)"""
    _pyplot()
    import ringspy
    scan_file = work / 'wood_ct.npy'
    if not scan_file.is_file():
        cells = ringspy.load_config(EXAMPLES / 'wood_cube.toml').replace(
            r_max=2.0, nrings=3, box_size=0.6, box_center=[1.1, 0.0], z_max=0.6, iter_max=4000, seed=2023)
        scan = ringspy.voxelize(cells, voxel_size=0.002, tilt=[0.05, 0.02], twist=0.02, wobble_amplitude=0.006,
                                wobble_wavelength=0.4, noise=0.15, blur=0.8, seed=2023, fiber_length=1.5,
                                fiber_length_cov=0.4)
        np.save(scan_file, np.clip(255*scan.volume, 0, 255).astype(np.uint8))
    config = ringspy.ReconstructionConfig.from_toml(EXAMPLES / 'reconstruction' / 'wood_ct.toml').replace(
        volume=str(scan_file), path=str(work / 'ct'), stlFlag=False, inpFlag=False)
    result = ringspy.reconstruct(config)
    shutil.copyfile(Path(result.output_dir) / (config.geoName + '.png'), OUT / 'wood_ct_preview.png')


# ============================================================================ Lloyd's relaxation
def lloyd_gif(work, size=720, delay_ms=100):
    """lloyd.gif: Lloyd's relaxation of the cells of the heart as CellPlacement_Binary_Lloyd does it for
    wood_cube_lloyd.toml: the cells on the perimeter of the heart stay, the cells inside start in a tiny
    square at the pith and move to the centroids of their Voronoi cells, iter_max times"""
    plt = _pyplot()
    from matplotlib.collections import LineCollection
    from PIL import Image
    from scipy.spatial import Voronoi
    import ringspy
    from ringspy import MeshGenTools as rp
    from ringspy.plot import COLORS, RING_COLORS, STYLE

    c = ringspy.load_config(EXAMPLES / 'wood_cube_lloyd.toml')
    rng = rp._rng(c.seed)
    widths = np.concatenate(([c.width_heart], np.tile([c.width_sparse, c.width_dense], c.nrings)))
    R = float(widths[0]*rng.normal(1, 0.25, len(widths))[0])     # the radius of the heart
    cell = float(np.atleast_1d(c.cellsize_dense)[0])
    npoint = int(np.ceil(2*np.pi*R/(2*cell)))
    t = np.linspace(0, 2*np.pi, npoint, endpoint=False)
    perimeter = np.c_[R*np.cos(t), R*np.sin(t)]
    inner = 1e-3*(rng.random((int(2*np.floor((R/cell)**2)), 2)) - [0.5, 0.5])
    states = [inner]
    for _ in range(c.iter_max):
        inner = rp.relax_points(Voronoi(np.vstack((perimeter, inner))), 1.0)
        states.append(inner)

    reach = 1.32*R                  # the drawing: a disc a little larger than the heart
    ink, guide = COLORS['boundary'], COLORS['rings']
    outer = '#e2ddd6'               # the walls between the cells on the perimeter, faint, out to 1.18 R
    site_color = RING_COLORS[0][0]
    frames = []
    with plt.rc_context(STYLE):
        for k, inner in enumerate(states):
            points = np.vstack((perimeter, inner))
            vor = Voronoi(points)
            segments, outside = [], []
            center = points.mean(axis=0)
            for (p, q), ridge in zip(vor.ridge_points, vor.ridge_vertices):
                ridge = np.asarray(ridge)
                if np.all(ridge >= 0):
                    seg = vor.vertices[ridge]
                else:                        # an infinite ridge (between two perimeter cells): out to `reach`
                    i = ridge[ridge >= 0][0]
                    tangent = points[q] - points[p]
                    normal = np.array([-tangent[1], tangent[0]])/np.linalg.norm(tangent)
                    midpoint = points[[p, q]].mean(axis=0)
                    direction = np.sign(np.dot(midpoint - center, normal))*normal
                    seg = np.array([vor.vertices[i], vor.vertices[i] + direction*2*reach])
                both_perimeter = p < npoint and q < npoint
                (outside if both_perimeter else segments).append(seg)
            fig = plt.figure(figsize=(size/100, size/100), dpi=100)
            ax = fig.add_axes([0, 0, 1, 1])
            ax.set_xlim(-1.08*reach, 1.08*reach)
            ax.set_ylim(-1.08*reach, 1.08*reach)
            ax.set_aspect('equal')
            ax.set_axis_off()
            for segs, color, width, clip in ((outside, outer, 0.9, 1.18*R), (segments, ink, 1.1, reach)):
                lines = LineCollection(segs, colors=color, linewidths=width, capstyle='round')
                lines.set_clip_path(plt.Circle((0, 0), clip, transform=ax.transData))
                ax.add_collection(lines)
            ax.add_patch(plt.Circle((0, 0), R, fill=False, edgecolor=guide, linewidth=0.9, linestyle=(0, (4, 3))))
            ax.scatter(perimeter[:, 0], perimeter[:, 1], s=9, color=COLORS['muted'], linewidths=0, zorder=3)
            ax.scatter(inner[:, 0], inner[:, 1], s=9, color=site_color, linewidths=0, zorder=3)
            ax.text(0.04, 0.955, "Lloyd's relaxation of the heart", transform=ax.transAxes, ha='left', va='top',
                    fontsize=15, fontweight='semibold', color=ink)
            ax.text(0.04, 0.912, 'iteration {:d} of {:d}'.format(k, len(states) - 1), transform=ax.transAxes,
                    ha='left', va='top', fontsize=12.5, color=COLORS['muted'])
            bar = 0.05                       # scale bar, mm
            x0, y0 = 1.08*reach - 0.06*reach - bar, -0.98*reach
            ax.plot([x0, x0 + bar], [y0, y0], color=ink, linewidth=1.6, solid_capstyle='butt')
            ax.text(x0 + bar/2, y0 + 0.025*reach, '{:g} mm'.format(bar), ha='center', va='bottom', fontsize=11.5,
                    color=COLORS['muted'])
            fig.canvas.draw()
            frames.append(Image.fromarray(np.asarray(fig.canvas.buffer_rgba())[:, :, :3].copy()))
            plt.close(fig)
    durations = [delay_ms]*len(frames)
    durations[0], durations[-1] = 700, 2000           # a pause on the start and on the relaxed cells
    _write_gif(frames, durations, OUT / 'lloyd.gif', ncolors=48)


def _write_gif(frames, durations, path, ncolors):
    """GIF with one palette for all frames (a median cut of a sample of their pixels), no dithering"""
    from PIL import Image
    sample = np.concatenate([np.asarray(f).reshape(-1, 3)[::7] for f in frames[::3] + frames[-1:]])
    strip = Image.fromarray(sample.reshape(1, -1, 3).astype(np.uint8))
    palette = strip.quantize(colors=ncolors, method=Image.Quantize.MEDIANCUT, dither=Image.Dither.NONE)
    indexed = [f.quantize(palette=palette, dither=Image.Dither.NONE) for f in frames]
    indexed[0].save(path, save_all=True, append_images=indexed[1:], duration=durations, loop=0, optimize=True,
                    disposal=1)
    print('{}: {} frames, {:.1f} s, {:,} bytes'.format(path.name, len(frames), sum(durations)/1000,
                                                       path.stat().st_size))


# ============================================================================ growth rule
def growth_rule():
    """growth_rule_binary.png: the radii of the binary growth rule (heart, then sparse and dense parts of
    every ring) and the rings they make, with the clipping box"""
    plt = _pyplot()
    from matplotlib.patches import Circle, FancyArrowPatch, Rectangle, Wedge
    from ringspy.plot import COLORS, STYLE
    ink, muted, text = COLORS['boundary'], COLORS['muted'], COLORS['text']
    sparse_fill, dense_fill = '#f2e9dc', '#a7835f'       # the earlywood and latewood of the palette, as fills
    edge = '#8c7f72'
    rc = dict(STYLE, **{'mathtext.fontset': 'custom'})   # math in the same sans as the text
    with plt.rc_context(rc):
        fig = plt.figure(figsize=(13.4, 4.1), dpi=140)
        ax = fig.add_axes([0, 0, 1, 1])
        ax.set_xlim(0, 13.4)
        ax.set_ylim(0, 4.1)
        ax.set_axis_off()

        def arrow(a, b, both=False, color=ink, lw=1.1):
            ax.add_patch(FancyArrowPatch(a, b, arrowstyle='<|-|>' if both else '-|>', mutation_scale=11,
                                         color=color, linewidth=lw, shrinkA=0, shrinkB=0))

        # the radii along a line: heart, then sparse and dense parts of two rings
        x = [1.55, 2.45, 4.15, 4.95, 6.65, 7.45]          # 0, r1, ..., r5
        y0, y1 = 1.72, 2.02                                 # the band
        fills = [dense_fill, sparse_fill, dense_fill, sparse_fill, dense_fill]
        names = ['HW', 'SW', 'DW', 'SW', 'DW']
        for a, b, fill, name in zip(x[:-1], x[1:], fills, names):
            ax.add_patch(Rectangle((a, y0), b - a, y1 - y0, facecolor=fill, edgecolor=edge, linewidth=0.8))
            ax.text((a + b)/2, y1 + 0.1, name, ha='center', va='bottom', fontsize=12.5, fontweight='semibold',
                    color=text)
        ax.add_patch(Rectangle((x[-1], y0), 0.55, y1 - y0, facecolor=sparse_fill, edgecolor='none'))
        ax.plot([x[-1], x[-1] + 0.55], [y0, y0], color=edge, linewidth=0.8)
        ax.plot([x[-1], x[-1] + 0.55], [y1, y1], color=edge, linewidth=0.8)
        ax.text(x[-1] + 0.75, (y0 + y1)/2, '...', ha='left', va='center', fontsize=13, color=muted)
        for (a, b), label in (((x[1], x[3]), 'ring 1'), ((x[3], x[5]), 'ring 2')):
            yb = y1 + 0.62
            ax.plot([a + 0.03, a + 0.03, b - 0.03, b - 0.03], [yb - 0.1, yb, yb, yb - 0.1], color=ink, linewidth=1.0)
            ax.text((a + b)/2, yb + 0.08, label, ha='center', va='bottom', fontsize=13, color=ink)
        labels = ['0', '$r_1$', '$r_2$', '$r_3$', '$r_4$', '$r_5$']
        for xi, label in zip(x, labels):
            ax.plot([xi, xi], [y0 - 0.1, y0], color=ink, linewidth=1.0)
            ax.text(xi, y0 - 0.17, label, ha='center', va='top', fontsize=14, color=ink)
        ax.text(x[0] - 0.22, y0 - 0.17, 'radii = [', ha='right', va='top', fontsize=14, color=ink)
        ax.text(x[-1] + 0.75, y0 - 0.17, '...]', ha='left', va='top', fontsize=14, color=ink)
        arrow((x[0], 0.72), (x[-1] + 1.15, 0.72))
        ax.plot([x[0], x[0]], [0.64, 0.80], color=ink, linewidth=1.0)
        ax.text(x[0], 0.55, '$r_{min} = 0$', ha='center', va='top', fontsize=13.5, color=ink)
        ax.text(x[-1] + 1.15, 0.55, '$r_{max}$', ha='center', va='top', fontsize=13.5, color=ink)
        # the key and the abbreviations
        for i, (fill, label) in enumerate(((sparse_fill, 'sparse part (earlywood)'), (dense_fill, 'dense part (latewood)'))):
            ax.add_patch(Rectangle((0.4, 3.68 - 0.30*i), 0.42, 0.2, facecolor=fill, edgecolor=edge, linewidth=0.8))
            ax.text(0.95, 3.78 - 0.30*i, label, ha='left', va='center', fontsize=12.5, color=text)
        for i, (name, label) in enumerate((('HW', 'heart width'), ('SW', 'sparse width'), ('DW', 'dense width'))):
            ax.text(5.0, 3.78 - 0.30*i, name, ha='left', va='center', fontsize=12.5, fontweight='semibold', color=text)
            ax.text(5.53, 3.78 - 0.30*i, '= ' + label, ha='left', va='center', fontsize=12.5, color=text)

        # the rings they make, and a clipping box
        O = np.array([10.75, 2.1])
        radii = [0.33, 0.8, 1.0, 1.5, 1.7]
        ax.add_patch(Circle(O, radii[-1], facecolor=dense_fill, edgecolor=edge, linewidth=0.8))
        for r, fill in zip(radii[::-1][1:], [sparse_fill, dense_fill, sparse_fill, dense_fill]):
            ax.add_patch(Circle(O, r, facecolor=fill, edgecolor=edge, linewidth=0.8))
        box_center, box_size = O + np.array([0.78, 0.1]), 1.1
        ax.add_patch(Rectangle(box_center - box_size/2, box_size, box_size, fill=False, edgecolor=ink, linewidth=1.6))
        ax.plot(*box_center, 'o', color=ink, markersize=4)
        ax.text(box_center[0] + 0.08, box_center[1] + 0.08, '$O_{box}$', ha='left', va='bottom', fontsize=13.5, color=ink)
        xd = box_center[0] + box_size/2 + 0.2
        arrow((xd, box_center[1] - box_size/2), (xd, box_center[1] + box_size/2), both=True)
        ax.text(xd + 0.1, box_center[1], '$d_{box}$', ha='left', va='center', fontsize=13.5, color=ink)
        ax.plot(*O, 'o', color=ink, markersize=4)
        ax.text(O[0] + 0.06, O[1] - 0.1, '$O$', ha='left', va='top', fontsize=13.5, color=ink)
        # each radius labelled where the background is light: r1 and r2 inside the sparse ring before them,
        # r3 in the sparse ring after it, r_max outside
        for r, angle, label, offset, align in ((radii[0], 62, '$r_1$', (0.1, 0.0), ('left', 'center')),
                                               (radii[1], 104, '$r_2$', (-0.1, -0.16), ('right', 'center')),
                                               (radii[2], 132, '$r_3$', (-0.04, 0.07), ('right', 'bottom')),
                                               (radii[4], 215, '$r_{max}$', (-0.05, -0.08), ('center', 'top'))):
            end = O + r*np.array([np.cos(np.radians(angle)), np.sin(np.radians(angle))])
            arrow(O, end)
            ax.text(*(end + np.array(offset)), label, ha=align[0], va=align[1], fontsize=13.5, color=ink)
        ax.text(O[0] + 1.65, O[1] - 1.72, '$r_{min} = 0$', ha='left', va='center', fontsize=13.5, color=ink)
        ax.text(6.7, 0.12, 'Radial growth rule: binary', ha='center', va='bottom', fontsize=13, color=muted)
        fig.savefig(OUT / 'growth_rule_binary.png', dpi=140, facecolor='white')
        plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    parser.add_argument('--work', type=Path, required=True, help='folder for the meshes, the scan and the previews')
    parser.add_argument('--only', nargs='+', choices=STEPS, default=STEPS, help='figures to make (default: all)')
    args = parser.parse_args()
    args.work.mkdir(parents=True, exist_ok=True)
    if 'previews' in args.only:
        previews_2d(args.work)
        preview_figure(args.work)
    if 'views' in args.only:
        views_3d(args.work)
    if 'ct' in args.only:
        ct_preview(args.work)
    if 'lloyd' in args.only:
        lloyd_gif(args.work)
    if 'growth' in args.only:
        growth_rule()
    print('written to', OUT)


if __name__ == '__main__':
    main()
