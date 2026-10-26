# -*- coding: utf-8 -*-
"""
Build the README animation contents/ringspy_hero.gif from real RingsPy output.

    python contents/make_hero_gif.py prepare --data <folder>   # generate the mesh, keep what the frames need
    python contents/make_hero_gif.py render  --data <folder>   # render the frames and write the GIF

prepare runs ringspy.generate() with contents/hero_log.toml (37k cells; about two minutes and 3 GB of
memory, the temporary mesh files take ~1 GB) and writes ~250 MB of compact data to <folder>. render needs
vtk, scipy and Pillow and an OpenGL context for off-screen rendering (about a minute with a GPU); it writes
the frames to <folder>/frames/ and the GIF next to this script, and checks the GIF (decoded frames, size).

Story of the loop: the finished twisted log orbits while it is wiped from the .vtu colours to the .stl
solid and the .inp mesh, where the orbit comes to rest; cell seeds are placed ring by ring; Voronoi cell
walls replace them, the camera zooms down to single cells and back; the log is clipped with a dog-bone, a
star and an irregular polygon, then the hexagonal specimen is clipped; the camera tilts over it while it
is extruded layer by layer with the twist, and eases back into the orbit of the first frame.
"""

import argparse
import io
import json
import math
import tempfile
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
TOML = HERE / 'hero_log.toml'
GIF = HERE / 'ringspy_hero.gif'

# Polygons the log is clipped with in turn before the hexagon of hero_log.toml (mm, pith at the origin):
# a dog-bone tensile specimen, a star and an irregular polygon; all non-convex (box_shape = 'polygon')
_star_angles = np.pi/2 + np.arange(10)*np.pi/5
CLIP_POLYGONS = {
    'dogbone': [[-5.0, -2.6], [-3.0, -2.6], [-2.0, -1.3], [2.0, -1.3], [3.0, -2.6], [5.0, -2.6],
                [5.0, 2.6], [3.0, 2.6], [2.0, 1.3], [-2.0, 1.3], [-3.0, 2.6], [-5.0, 2.6]],
    'star': np.column_stack((np.where(np.arange(10) % 2 == 0, 5.2, 2.3)*np.cos(_star_angles),
                             np.where(np.arange(10) % 2 == 0, 5.2, 2.3)*np.sin(_star_angles))).round(4).tolist(),
    'irregular': [[-4.6, -2.2], [-2.0, -4.3], [1.2, -3.4], [3.0, -4.6], [4.8, -1.6], [2.6, 0.4], [4.4, 3.0],
                  [1.4, 4.6], [-0.8, 2.2], [-3.2, 4.2], [-4.9, 1.0]],
}

W, H = 800, 450           # GIF size (the layout below scales with W)
K = W/800
SS = 2                    # supersampling of the 2D drawings
SS3D = 4                  # supersampling of the 3D renders (sub-pixel cell walls shimmer less)
DELAY_MS = 80             # 12.5 frames per second
CENTER_PX = (round(535*K), round(240*K))    # where the log centre goes on the canvas
CAPTION_BOX = tuple(round(v*K) for v in (12, 6, 322, 172))     # wordmark and captions (x0, y0, x1, y1)
LEGEND_BOX = tuple(round(v) for v in (12*K, H - 66*K, 322*K, H - 2*K))   # colour key; the object
                                                                   # must stay clear of both

BG = (0x2B, 0x2F, 0x3D)
EARLY = (0xA9, 0xC1, 0xEE)        # thin cell walls (Connector_width = smaller value)
LATE = (0xE0, 0x55, 0x45)         # thick cell walls (Connector_width = larger value)
GUIDE = (0x4A, 0x50, 0x68)
LIGHT = (0xE6, 0xE9, 0xF0)
CAPTION = (0xB8, 0xBF, 0xD0)
IVORY = (0xEC, 0xE3, 0xCF)
BEAM = (0x8E, 0x96, 0xAB)
KEY_COLORS = [BG, EARLY, LATE, GUIDE, LIGHT, CAPTION, IVORY, BEAM]


# ============================================================================ prepare
def prepare(data, work):
    import ringspy
    from ringspy import MeshGenTools as rp

    data.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    config = ringspy.Config.from_toml(TOML).replace(path=str(work))
    res = ringspy.generate(config)
    c, out, name = res.config, res.output_dir, res.config.geoName
    print('generate: {} cells in {:.1f} s'.format(res.nParticles, time.time() - t0))

    polygon = rp.Clipping_Box(c.box_shape, c.box_center, c.box_size, c.boundaryFlag)[5]
    np.savez_compressed(data / 'cells.npz', sites=res.sites, radii=res.radii, polygon=polygon,
                        thickness=np.array([c.cellwallthickness_sparse, c.cellwallthickness_dense]),
                        center=np.array(c.generation_center, dtype=float))
    check_latewood_walls(out / (name + '-vertex.mesh'), res.radii, c)
    clipped_meshes(res, c, data / 'clips.npz')
    nfaces, levels = connector_volume_surface(out / (name + '_conns_vol.vtu'), data / 'conns_vol.vtp')
    nbeams, nconns = inp_lines(out / (name + '.inp'), data / 'inp_lines.vtp')
    info = dict(cells=int(res.nParticles), beam_elements=int(res.nbeamElem), porosity=float(res.porosity),
                surface_faces=nfaces, zc_levels=levels, inp_beams=nbeams, inp_top_connectors=nconns,
                z_min=c.z_min, z_max=c.z_max, theta_max=c.theta_max)
    (data / 'info.json').write_text(json.dumps(info, indent=1))
    print('prepare done in {:.1f} s: {}'.format(time.time() - t0, info))


def clipped_meshes(res, c, path):
    """RingsPy's 2D mesh of the same log clipped with each of CLIP_POLYGONS and with the box of
    hero_log.toml (the same functions as ringspy.generate: Clipping_Box and RebuildVoronoi(_merge))"""
    import matplotlib.pyplot as plt
    from scipy.spatial import Voronoi
    from ringspy import MeshGenTools as rp
    vor = Voronoi(res.sites[:, 0:2])
    shapes = [(key, 'polygon', polygon) for key, polygon in CLIP_POLYGONS.items()] + [('box', c.box_shape, None)]
    arrays = {}
    for key, shape, polygon in shapes:
        x_min, x_max, y_min, y_max, boundaries, corners, _ = \
            rp.Clipping_Box(shape, c.box_center, c.box_size, c.boundaryFlag, polygon=polygon)
        if c.merge_operation:
            out = rp.RebuildVoronoi_merge(vor, res.sites, boundaries, c.generation_center, x_min, x_max, y_min, y_max,
                                          c.box_center, shape, c.merge_tol, c.boundaryFlag)
        else:
            out = rp.RebuildVoronoi(vor, res.sites, boundaries, c.generation_center, x_min, x_max, y_min, y_max,
                                    c.box_center, shape, c.boundaryFlag)
        plt.close('all')
        arrays['polygon_' + key], arrays['vertices_' + key], arrays['ridges_' + key] = corners, out[0], out[11]
        print('clipped with {}: {} vertices, {} walls'.format(key, len(out[0]), len(out[11])))
    np.savez_compressed(path, order=np.array(list(CLIP_POLYGONS) + ['box']), **arrays)


def ring_bins(points, radii, center):
    """0 = heart, odd = rings placed with cellsize_sparse, even >= 2 = rings placed with cellsize_dense"""
    r = np.linalg.norm(np.asarray(points)[:, :2] - center, axis=1)
    return np.searchsorted(radii, r, side='right') - 1


def check_latewood_walls(vertex_mesh, radii, c):
    """The colours of the animation assume that the dense (latewood) rings got the thicker walls."""
    lines = vertex_mesh.read_text().splitlines()
    max_wings = int(lines[4])
    table = np.loadtxt(io.StringIO('\n'.join(lines[6:])))
    width = table[:, 3 + 3*max_wings]
    bins = ring_bins(table[:, :2], radii, np.asarray(c.generation_center, dtype=float))
    dense = width[(bins >= 2) & (bins % 2 == 0)].mean()
    sparse = width[bins % 2 == 1].mean()
    print('mean wall thickness: dense rings {:.4f}, sparse rings {:.4f}'.format(dense, sparse))
    if not dense > sparse:
        raise RuntimeError('The dense (latewood) rings do not have the thicker walls: check '
                           'cellwallthickness_sparse/_dense in contents/hero_log.toml.')


def connector_volume_surface(vtu, vtp):
    """Outer surface of the connector volumes (the solid that StlModelFile exports), keeping
    Connector_width and adding zc, the mean z of each volume, for the layer-by-layer build-up."""
    from vtkmodules.util.numpy_support import numpy_to_vtk, vtk_to_numpy
    from vtkmodules.vtkFiltersGeometry import vtkDataSetSurfaceFilter
    from vtkmodules.vtkIOXML import vtkXMLUnstructuredGridReader

    reader = vtkXMLUnstructuredGridReader()
    reader.SetFileName(str(vtu))
    reader.Update()
    grid = reader.GetOutput()
    z = vtk_to_numpy(grid.GetPoints().GetData())[:, 2]
    conn = vtk_to_numpy(grid.GetCells().GetConnectivityArray()).reshape(-1, 8)   # all hexahedra
    zc = np.round(z[conn].mean(axis=1), 6).astype(np.float32)
    array = numpy_to_vtk(zc, deep=True)
    array.SetName('zc')
    grid.GetCellData().AddArray(array)
    surface = vtkDataSetSurfaceFilter()
    surface.SetInputData(grid)
    surface.Update()
    write_vtp(surface.GetOutput(), vtp)
    return surface.GetOutput().GetNumberOfCells(), int(len(np.unique(zc)))


def inp_lines(inp, vtp):
    """Line drawing of the Abaqus input file: all B32 beams (grains) and the transverse
    connectors (T3D2) of the top layer."""
    from vtkmodules.util.numpy_support import numpy_to_vtk, numpy_to_vtkIdTypeArray
    from vtkmodules.vtkCommonCore import vtkPoints
    from vtkmodules.vtkCommonDataModel import vtkCellArray, vtkPolyData

    blocks, current = {}, None
    with open(inp) as f:
        for line in f:
            if line.startswith('*'):
                key = line.lower().replace(' ', '').strip()
                current = None
                if key.startswith('*node,nset=allnodes'):
                    current = 'nodes'
                elif key.startswith('*element,type=b32,elset=allbeams'):
                    current = 'beams'
                elif key.startswith('*element,type=t3d2,elset=alltopconns'):
                    current = 'top'
                if current:
                    blocks[current] = []
            elif current:
                blocks[current].append(line)
    parse = lambda name: np.loadtxt(io.StringIO(''.join(blocks[name])), delimiter=',')
    nodes = parse('nodes')
    coords = np.zeros((int(nodes[:, 0].max()) + 1, 3))
    coords[nodes[:, 0].astype(int)] = nodes[:, 1:4]
    beams = parse('beams')[:, 1:4].astype(np.int64)
    top = parse('top')[:, 1:3].astype(np.int64)
    zmax = coords[:, 2].max()
    top = top[np.all(coords[top, 2] > zmax - 1e-9, axis=1)]

    points = vtkPoints()
    points.SetData(numpy_to_vtk(coords, deep=True))
    offsets = np.concatenate(([0], np.cumsum([3]*len(beams) + [2]*len(top))))
    connectivity = np.concatenate((beams.ravel(), top.ravel()))
    cells = vtkCellArray()
    cells.SetData(numpy_to_vtkIdTypeArray(offsets, deep=True), numpy_to_vtkIdTypeArray(connectivity, deep=True))
    poly = vtkPolyData()
    poly.SetPoints(points)
    poly.SetLines(cells)
    kind = numpy_to_vtk(np.concatenate((np.zeros(len(beams)), np.ones(len(top)))).astype(np.uint8), deep=True)
    kind.SetName('kind')    # 0: beam, 1: top-layer connector
    poly.GetCellData().AddArray(kind)
    write_vtp(poly, vtp)
    return len(beams), len(top)


def write_vtp(polydata, path):
    from vtkmodules.vtkIOXML import vtkXMLPolyDataWriter
    writer = vtkXMLPolyDataWriter()
    writer.SetFileName(str(path))
    writer.SetInputData(polydata)
    writer.SetDataModeToBinary()
    writer.SetCompressorTypeToZLib()
    writer.Write()


# ============================================================================ render: 3D
class Renderer3D:
    """Off-screen VTK renders of the connector volumes (.vtu colours), the same solid as the STL,
    and the line drawing of the .inp, all with a parallel camera at the 2D scale."""

    def __init__(self, data, scale_px_per_mm, polygon):
        from vtkmodules.util.numpy_support import numpy_to_vtk, vtk_to_numpy
        from vtkmodules.vtkIOXML import vtkXMLPolyDataReader
        from vtkmodules.vtkRenderingCore import vtkLight, vtkRenderer, vtkRenderWindow
        import vtkmodules.vtkRenderingOpenGL2  # noqa: F401  (OpenGL backend)

        self.vtk_to_numpy = vtk_to_numpy
        read = lambda name: self._read(vtkXMLPolyDataReader, data / name)
        self.solid = read('conns_vol.vtp')
        self.lines = read('inp_lines.vtp')
        width = vtk_to_numpy(self.solid.GetCellData().GetArray('Connector_width'))
        thick = width > (width.min() + width.max()) / 2
        colors = numpy_to_vtk(np.where(thick[:, None], LATE, EARLY).astype(np.uint8), deep=True)
        colors.SetName('rgb')
        self.solid.GetCellData().AddArray(colors)
        self.lines = self._lattice(self.lines, polygon)
        self.zc = vtk_to_numpy(self.solid.GetCellData().GetArray('zc'))
        self.levels = np.unique(self.zc)
        b = self.solid.GetBounds()
        self.zmid = (b[4] + b[5]) / 2
        self.size = max(b[1] - b[0], b[3] - b[2], b[5] - b[4])

        self.renderer = vtkRenderer()
        self.renderer.SetBackground(*[v/255 for v in BG])
        self.renderer.AutomaticLightCreationOff()
        key = vtkLight()
        key.SetLightTypeToSceneLight()
        key.SetPosition(0.6*self.size, -0.8*self.size, 3.0*self.size)
        key.SetFocalPoint(0, 0, 0)
        key.SetIntensity(0.85)
        head = vtkLight()
        head.SetLightTypeToHeadlight()
        head.SetIntensity(0.35)
        self.renderer.AddLight(key)
        self.renderer.AddLight(head)
        self.window = vtkRenderWindow()
        self.window.SetOffScreenRendering(1)
        self.window.SetMultiSamples(0)
        self.window.AddRenderer(self.renderer)
        self.window.SetSize(W*SS3D, H*SS3D)
        self.parallel_scale = (H/2) / scale_px_per_mm
        self.solid_actor = self._actor(self.solid, flat_color=None)
        self.stl_actor = self._actor(self.solid, flat_color=IVORY)
        # the .inp lines; an invisible lid (the clipping polygon turned by the twist of the top layer,
        # in the background colour) hides the beams behind the top face, as a solid top would (polygon
        # offset lets the top-layer connectors win against the lid)
        from vtkmodules.vtkRenderingCore import vtkMapper
        vtkMapper.SetResolveCoincidentTopologyToPolygonOffset()
        occluder = self._actor(self._lid(polygon), flat_color=BG)
        occluder.GetProperty().LightingOff()
        self.line_actor = [occluder, self._actor(self.lines, flat_color=None, lines=True),
                           self._outline_actor(self.bottom_outline)]
        self.layers = {}

    def _lid(self, polygon):
        """The clipping polygon at the top of the model, turned by the twist of the top layer (measured
        on the .inp nodes: the median turn of each grain from the bottom to the top layer)"""
        from vtkmodules.util.numpy_support import numpy_to_vtk, numpy_to_vtkIdTypeArray
        from vtkmodules.vtkCommonCore import vtkPoints
        from vtkmodules.vtkCommonDataModel import vtkCellArray, vtkPolyData
        nodes = self.vtk_to_numpy(self.lines.GetPoints().GetData())[1:]     # point 0 is unused
        npl = int(np.isclose(nodes[:, 2], nodes[:, 2].min()).sum())
        bottom, top = nodes[:npl], nodes[-npl:]
        turn = np.angle(np.exp(1j*(np.arctan2(top[:, 1], top[:, 0]) - np.arctan2(bottom[:, 1], bottom[:, 0]))))
        turn = np.median(turn[np.linalg.norm(bottom[:, :2], axis=1) > 0.5])
        c, s_ = np.cos(turn), np.sin(turn)
        corners = np.column_stack((polygon @ np.array([[c, s_], [-s_, c]]), np.full(len(polygon), top[:, 2].max())))
        points = vtkPoints()
        points.SetData(numpy_to_vtk(corners, deep=True))
        cells = vtkCellArray()
        cells.SetData(numpy_to_vtkIdTypeArray(np.array([0, len(polygon)], dtype=np.int64), deep=True),
                      numpy_to_vtkIdTypeArray(np.arange(len(polygon), dtype=np.int64), deep=True))
        lid = vtkPolyData()
        lid.SetPoints(points)
        lid.SetPolys(cells)
        self.top_turn = turn
        return lid

    def _lattice(self, lines, polygon):
        """The .inp line drawing reduced to what can be seen: all top-layer connectors (the cell mesh
        of the top face, drawn a little dimmer) and, on the side faces, the beams of the outermost
        grain (vertex column) in every 0.2 mm of the perimeter, which trace the twisted grain; every 5th
        of them and the grains at the 6 corners are drawn light. Node numbers follow the .inp:
        layer*npl + vertex + 1."""
        from vtkmodules.util.numpy_support import numpy_to_vtk, numpy_to_vtkIdTypeArray
        from vtkmodules.vtkCommonDataModel import vtkCellArray, vtkPolyData
        nodes = self.vtk_to_numpy(lines.GetPoints().GetData())
        npl = int(np.isclose(nodes[1:, 2], nodes[1:, 2].min()).sum())    # nodes per layer (0 is unused)
        kind = self.vtk_to_numpy(lines.GetCellData().GetArray('kind'))
        offsets = self.vtk_to_numpy(lines.GetLines().GetOffsetsArray())
        conn = self.vtk_to_numpy(lines.GetLines().GetConnectivityArray())
        grain = (conn[offsets[:-1]] - 1) % npl
        beam = kind == 0

        # position of every grain along the clipping polygon's perimeter, and its distance to it
        grains = np.unique(grain[beam])
        xy = nodes[grains + 1, :2]                                       # grain position at z_min
        a, b = polygon, np.roll(polygon, -1, axis=0)
        lengths = np.linalg.norm(b - a, axis=1)
        u = np.clip(np.einsum('nmk,mk->nm', xy[:, None] - a, b - a)/lengths**2, 0, 1)
        dist = np.linalg.norm(xy[:, None] - (a + u[..., None]*(b - a)), axis=2)
        edge = dist.argmin(axis=1)
        to_box = dist.min(axis=1)
        arc = np.concatenate(([0], np.cumsum(lengths)))[edge] + u[np.arange(len(grains)), edge]*lengths[edge]
        near = to_box < 0.09
        slot = np.floor(arc/0.2).astype(int)
        order = np.lexsort((to_box, slot))                               # per slot, the outermost first
        first = order[near[order]]
        _, pick = np.unique(slot[first], return_index=True)
        chosen = grains[first[pick]]
        corners = grains[np.linalg.norm(xy[:, None] - polygon, axis=2).argmin(axis=0)]
        light = np.concatenate((chosen[::5], corners))
        side = beam & np.isin(grain, chosen)
        bright = beam & np.isin(grain, light)
        keep = ~beam | side | bright

        sizes = np.diff(offsets)[keep]
        pieces = [conn[i:j] for i, j, k in zip(offsets[:-1], offsets[1:], keep) if k]
        cells = vtkCellArray()
        cells.SetData(numpy_to_vtkIdTypeArray(np.concatenate(([0], np.cumsum(sizes))).astype(np.int64), deep=True),
                      numpy_to_vtkIdTypeArray(np.concatenate(pieces).astype(np.int64), deep=True))
        subset = vtkPolyData()
        subset.SetPoints(lines.GetPoints())
        subset.SetLines(cells)
        rgb = np.where(bright[:, None], LIGHT, np.where(beam[:, None], BEAM, CAPTION))[keep].astype(np.uint8)
        colors = numpy_to_vtk(np.ascontiguousarray(rgb), deep=True)
        colors.SetName('rgb')
        subset.GetCellData().AddArray(colors)
        self.bottom_outline = np.column_stack((polygon, np.full(len(polygon), nodes[1:, 2].min())))
        return subset

    def _outline_actor(self, corners):
        """A closed light line through `corners` (the bottom edges of the prism in the .inp drawing)"""
        from vtkmodules.util.numpy_support import numpy_to_vtk, numpy_to_vtkIdTypeArray
        from vtkmodules.vtkCommonCore import vtkPoints
        from vtkmodules.vtkCommonDataModel import vtkCellArray, vtkPolyData
        points = vtkPoints()
        points.SetData(numpy_to_vtk(np.ascontiguousarray(corners, dtype=float), deep=True))
        loop = np.append(np.arange(len(corners)), 0)
        cells = vtkCellArray()
        cells.SetData(numpy_to_vtkIdTypeArray(np.array([0, len(loop)], dtype=np.int64), deep=True),
                      numpy_to_vtkIdTypeArray(loop.astype(np.int64), deep=True))
        poly = vtkPolyData()
        poly.SetPoints(points)
        poly.SetLines(cells)
        return self._actor(poly, flat_color=LIGHT, lines=True)

    @staticmethod
    def _read(reader_class, path):
        reader = reader_class()
        reader.SetFileName(str(path))
        reader.Update()
        return reader.GetOutput()

    def _actor(self, dataset, flat_color, lines=False):
        from vtkmodules.vtkRenderingCore import vtkActor, vtkDataSetMapper
        mapper = vtkDataSetMapper()
        mapper.SetInputData(dataset)
        actor = vtkActor()
        actor.SetMapper(mapper)
        prop = actor.GetProperty()
        if flat_color is None:
            mapper.SetScalarModeToUseCellFieldData()
            mapper.SelectColorArray('rgb')
            mapper.SetColorModeToDirectScalars()
            mapper.ScalarVisibilityOn()
        else:
            mapper.ScalarVisibilityOff()
            prop.SetColor(*[v/255 for v in flat_color])
        if lines:
            prop.LightingOff()
            prop.SetLineWidth(SS3D/2)
        else:
            prop.SetInterpolationToFlat()
            prop.SetAmbient(0.25)
            prop.SetDiffuse(0.8)
            prop.SetSpecular(0.0)
        return actor

    def layer_actor(self, level):
        """Connector volumes up to (and including) z level `level` (only the last one is kept: the
        full model has millions of faces)"""
        if level not in self.layers:
            from vtkmodules.util.numpy_support import numpy_to_vtk, numpy_to_vtkIdTypeArray
            from vtkmodules.vtkCommonDataModel import vtkCellArray, vtkPolyData
            polys = self.solid.GetPolys()
            offsets = self.vtk_to_numpy(polys.GetOffsetsArray())
            conn = self.vtk_to_numpy(polys.GetConnectivityArray())
            assert np.all(np.diff(offsets) == 4), 'expected quads only'
            keep = self.zc <= level + 1e-5
            cells = vtkCellArray()
            cells.SetData(numpy_to_vtkIdTypeArray(np.arange(keep.sum() + 1, dtype=np.int64)*4, deep=True),
                          numpy_to_vtkIdTypeArray(conn.reshape(-1, 4)[keep].ravel().astype(np.int64), deep=True))
            subset = vtkPolyData()
            subset.SetPoints(self.solid.GetPoints())
            subset.SetPolys(cells)
            rgb = self.vtk_to_numpy(self.solid.GetCellData().GetArray('rgb'))[keep]
            colors = numpy_to_vtk(np.ascontiguousarray(rgb), deep=True)
            colors.SetName('rgb')
            subset.GetCellData().AddArray(colors)
            self.layers = {level: self._actor(subset, flat_color=None)}
        return self.layers[level]

    def render(self, actor, elevation, azimuth, focal_z):
        """Render one actor (or a list of actors); returns an (H*SS, W*SS, 3) uint8 image."""
        from vtkmodules.vtkRenderingCore import vtkWindowToImageFilter
        self.renderer.RemoveAllViewProps()
        for a in (actor if isinstance(actor, list) else [actor]):
            self.renderer.AddActor(a)
        el, az = math.radians(elevation), math.radians(azimuth)
        cam = self.renderer.GetActiveCamera()
        cam.ParallelProjectionOn()
        cam.SetParallelScale(self.parallel_scale)
        d = 4*self.size
        cam.SetFocalPoint(0, 0, focal_z)
        cam.SetPosition(d*math.cos(el)*math.cos(az), d*math.cos(el)*math.sin(az), focal_z + d*math.sin(el))
        cam.SetViewUp(-math.sin(el)*math.cos(az), -math.sin(el)*math.sin(az), math.cos(el))
        # put the focal point at CENTER_PX instead of the window centre
        cam.SetWindowCenter((W/2 - CENTER_PX[0]) / (W/2), (CENTER_PX[1] - H/2) / (H/2))
        self.renderer.ResetCameraClippingRange()
        self.window.Render()
        grab = vtkWindowToImageFilter()
        grab.SetInput(self.window)
        grab.ReadFrontBufferOff()
        grab.Update()
        image = grab.GetOutput()
        w, h, _ = image.GetDimensions()
        pixels = self.vtk_to_numpy(image.GetPointData().GetScalars()).reshape(h, w, -1)[::-1, :, :3]
        from PIL import Image
        return np.asarray(Image.fromarray(np.ascontiguousarray(pixels)).reduce(SS3D // SS))


# ============================================================================ render: 2D
class Drawer2D:
    """Top view drawings (cell nuclei, Voronoi walls, clipping outline) with the same mapping as the
    3D renders: x right, y up, pith at CENTER_PX, `scale` pixels per mm (before supersampling)."""

    def __init__(self, cells, scale):
        from scipy.spatial import Voronoi
        self.sites, self.radii, self.polygon = cells['sites'], cells['radii'], cells['polygon']
        self.center = cells['center']
        self.scale = scale
        thick_sparse, thick_dense = cells['thickness']
        self.bins = ring_bins(self.sites, self.radii, self.center)
        assert np.all(np.diff(self.bins) >= 0), 'sites are not in placement order'

        vor = Voronoi(self.sites[:, :2])
        ridges = np.array([r for r in vor.ridge_vertices if -1 not in r])
        segments = vor.vertices[ridges]                       # (n, 2, 2)
        keep = np.all(np.linalg.norm(segments - self.center, axis=2) < self.radii[-1], axis=1)
        self.segments = segments[keep]
        mid = self.segments.mean(axis=1)
        self.mid_r = np.linalg.norm(mid - self.center, axis=1)
        self.thick_sparse, self.thick_dense = thick_sparse, thick_dense
        self.thickness, self.thick = self.wall_thickness(self.segments)
        from matplotlib.path import Path as MplPath
        self.mid_inside = MplPath(self.polygon).contains_points(mid)

    def wall_thickness(self, segments):
        """RingsPy's rule: each vertex gets the thickness of its ring, bisect(radii, r) even (a ring placed
        with cellsize_sparse) -> cellwallthickness_sparse, else cellwallthickness_dense (VertexandRidgeinfo),
        and a wall takes the thickness of its second vertex (ConnectorMeshFile); scipy's ridge order is
        kept by RingsPy for walls inside the box. Returns the thickness of each wall and whether it is
        the thicker one."""
        r1 = np.linalg.norm(segments[:, 1] - self.center, axis=1)
        sparse = (np.searchsorted(self.radii, r1, side='right') % 2) == 0
        thickness = np.where(sparse, self.thick_sparse, self.thick_dense)
        return thickness, thickness == max(self.thick_sparse, self.thick_dense)

    def to_px(self, xy, view):
        """world (mm) -> supersampled pixels for a view (center_world, scale, center_px)"""
        (cx, cy), s, (px, py) = view
        xy = np.asarray(xy, dtype=float)
        return np.stack(((px + s*(xy[..., 0] - cx))*SS, (py - s*(xy[..., 1] - cy))*SS), axis=-1)

    def full_view(self):
        return (tuple(self.center), self.scale, CENTER_PX)

    def canvas(self):
        from PIL import Image
        return Image.new('RGB', (W*SS, H*SS), BG)

    def nuclei(self, image, indices, view, size=1.0):
        from PIL import ImageDraw
        draw = ImageDraw.Draw(image)
        p = self.to_px(self.sites[indices, :2], view)
        r = np.maximum(self.sites[indices, 2]*view[1]*SS*size, 0.8)
        dense = self.bins[indices] % 2 == 0
        for (x, y), rr, d in zip(p, r, dense):
            draw.ellipse((x - rr, y - rr, x + rr, y + rr), fill=LATE if d else EARLY)

    def guide_circle(self, image, radius, view):
        from PIL import ImageDraw
        draw = ImageDraw.Draw(image)
        (x, y), r = self.to_px(self.center, view), radius*view[1]*SS
        draw.ellipse((x - r, y - r, x + r, y + r), outline=GUIDE, width=SS)

    def walls(self, image, view, max_r=np.inf, only_inside=False, tint=1.0):
        """The walls of the whole log's Voronoi tessellation (with midpoint within max_r)"""
        sel = self.mid_r <= max_r
        if only_inside:
            sel &= self.mid_inside
        self.draw_walls(image, self.segments[sel], self.thickness[sel], self.thick[sel], view, tint)

    def mesh(self, image, vertices, ridges, view, tint=1.0):
        """The walls of a mesh made by RingsPy (vertices and ridges of RebuildVoronoi(_merge))"""
        segments = vertices[ridges]
        thickness, thick = self.wall_thickness(segments)
        self.draw_walls(image, segments, thickness, thick, view, tint)

    def draw_walls(self, image, segments, thickness, thick, view, tint):
        from PIL import ImageDraw
        draw = ImageDraw.Draw(image)
        p = self.to_px(segments, view)
        # cull to the canvas
        lo, hi = p.min(axis=1), p.max(axis=1)
        on = (hi[:, 0] >= 0) & (lo[:, 0] <= W*SS) & (hi[:, 1] >= 0) & (lo[:, 1] <= H*SS)
        # walls thinner than a pixel are drawn 1 px wide; with tint = 1 their colour is mixed with the
        # background by the fraction of the pixel they cover, so the drawing has the true wall density
        # (as the 3D renders); with tint = 0 they are drawn in full colour, which reads better far away
        px = thickness*view[1]*SS
        widths = np.maximum(np.rint(px), 1)
        cover = tint*np.minimum(px/widths, 1)[:, None] + (1 - tint)
        colors = np.where(thick[:, None], LATE, EARLY)
        fills = np.rint(np.array(BG) + cover*(colors - np.array(BG))).astype(int)
        for (a, b), w, c in zip(p[on], widths[on].astype(int), fills[on]):
            draw.line((a[0], a[1], b[0], b[1]), fill=tuple(c), width=int(w))

    def outline(self, image, view, nedges, polygon=None):
        """The first nedges edges of the clipping polygon (default: the box of hero_log.toml)"""
        from PIL import ImageDraw
        draw = ImageDraw.Draw(image)
        p = self.to_px(self.polygon if polygon is None else polygon, view)
        for i in range(nedges):
            a, b = p[i], p[(i + 1) % len(p)]
            draw.line((a[0], a[1], b[0], b[1]), fill=LIGHT, width=2*SS)


# ============================================================================ render: story
def smoothstep(u):
    u = min(max(u, 0.0), 1.0)
    return u*u*(3 - 2*u)


def ease_in_out_cubic(u):
    u = min(max(u, 0.0), 1.0)
    return 4*u**3 if u < 0.5 else 1 - (-2*u + 2)**3/2


def render(data, gif_path, frames_only=False, encode_only=False):
    from PIL import Image

    frames_dir = data / 'frames'
    if not encode_only:
        render_frames(data, frames_dir)
    frames = [Image.open(p).convert('RGB') for p in sorted(frames_dir.glob('frame_*.png'))]
    if not frames_only:
        encode(frames, gif_path, preview=data / 'first_frame.png')


def render_frames(data, frames_dir):
    from PIL import Image

    frames_dir.mkdir(parents=True, exist_ok=True)
    for old in frames_dir.glob('frame_*.png'):
        old.unlink()
    cells = np.load(data / 'cells.npz')
    polygon = cells['polygon']
    scale = 0.80*H / (polygon[:, 1].max() - polygon[:, 1].min())    # px per mm
    r3d = Renderer3D(data, scale, polygon)
    d2d = Drawer2D(cells, scale)
    fonts = load_fonts()
    frames = []
    overlaps = []      # frames whose object reaches into the caption or legend boxes (zoomed frames excepted)

    def emit(image, caption, fills_canvas=False, legend=True):
        """image: supersampled image of the object only; caption: list of text lines.
        fills_canvas: the object covers the whole canvas (zoom), the caption gets a backdrop.
        legend: show the colour key (frames coloured by wall thickness)."""
        if isinstance(image, np.ndarray):
            image = Image.fromarray(image)
        small = image.reduce(SS)
        if not fills_canvas:
            pixels = np.asarray(small)
            for x0, y0, x1, y1 in (CAPTION_BOX, LEGEND_BOX):
                n = int(np.any(pixels[y0:y1, x0:x1] != BG, axis=2).sum())
                if n:
                    overlaps.append((len(frames), n))
        frames.append((small, caption, fills_canvas, legend))

    def blend(a, b, t):
        a = a if isinstance(a, Image.Image) else Image.fromarray(a)
        b = b if isinstance(b, Image.Image) else Image.fromarray(b)
        return Image.blend(a, b, t)

    def wipe(a, b, u):
        """b replaces a from the left edge of the object to its right edge; a light line at the front"""
        a, b = np.asarray(a), np.asarray(b)
        on = np.any((a != BG) | (b != BG), axis=2)
        cols, rows = np.nonzero(on.any(axis=0))[0], np.nonzero(on.any(axis=1))[0]
        x = int(round(cols[0] + (cols[-1] - cols[0])*smoothstep(u)))
        out = a.copy()
        out[:, :x] = b[:, :x]
        r0, r1 = max(rows[0] - 8*SS, 0), rows[-1] + 8*SS
        out[r0:r1, max(x - SS, 0):x + SS] = LIGHT
        return out

    az0, el_hero, step = -60.0, 40.0, 1.2
    t = time.time()

    # S1: the finished log orbits; wipes from the .vtu colours to the .stl solid and the .inp drawing.
    # The caption switches when the wipe line passes the middle. The orbit slows to a stop on the .inp
    # drawing, which then stays still (still frames cost next to nothing in the GIF).
    views = {'vtu': r3d.solid_actor, 'stl': r3d.stl_actor, 'inp': r3d.line_actor}
    labels = {'vtu': '.vtu  for ParaView', 'stl': '.stl  for 3D printing', 'inp': '.inp  for Abaqus'}
    tagline = ['Generates 3D meshes', 'of cellular solids']
    plan = [('vtu', 10, None), ('stl', 6, 'vtu'), ('stl', 8, None), ('inp', 6, 'stl'), ('inp', 24, None)]
    n_slow, n_still = 10, 14                    # the last 24 frames: slowing down, then still
    speed = np.full(sum(n for _, n, _ in plan), step)
    speed[-n_slow - n_still:-n_still] = [step*(1 - smoothstep((j + 1)/(n_slow + 1))) for j in range(n_slow)]
    speed[-n_still:] = 0
    azimuths = az0 + np.concatenate(([0.0], np.cumsum(speed[:-1])))
    k = 0
    rendered = {}
    for kind, n, previous in plan:
        for i in range(n):
            az = azimuths[k]
            key = (kind, az)
            image = rendered[key] if key in rendered else r3d.render(views[kind], el_hero, az, r3d.zmid)
            rendered = {key: image}
            shown = kind
            if previous:
                old = r3d.render(views[previous], el_hero, az, r3d.zmid)
                image = wipe(old, image, (i + 1)/(n + 1))
                shown = kind if i >= n//2 else previous
            emit(image, tagline + [labels[shown]], legend=(shown == 'vtu'))
            k += 1
    print('S1 done ({} frames, {:.0f} s)'.format(len(frames), time.time() - t))

    # S2: the still .inp drawing fades into the ring layout of the next scene
    r_log = np.linalg.norm(d2d.sites[:, :2] - d2d.center, axis=1).max()
    wide = (tuple(d2d.center), min(scale, 0.43*H/r_log), CENTER_PX)
    rings = d2d.canvas()
    for r in d2d.radii[1:]:
        d2d.guide_circle(rings, r, wide)
    for i in range(4):
        emit(blend(image, rings, (i + 1)/4), [], legend=False)

    # S3: RingsPy draws the ring radii first, then places the cell seeds ring by ring (real order)
    caption = ['1  Cell seeds placed', 'ring by ring']
    nsites, nframes = len(d2d.sites), 32
    shown = 0
    base = rings.copy()
    for f in range(nframes):
        upto = int(round(nsites*((f + 4)/(nframes + 3))**1.3))
        new = np.arange(shown, upto)
        image = base.copy()
        d2d.nuclei(image, new, wide, size=0.6)
        d2d.nuclei(base, new, wide)
        shown = upto
        emit(image, caption)

    # S4: Voronoi walls replace the seeds from the pith outward, zoom down to single cells and back
    # out to the specimen scale, clip the hexagon
    caption = ['2  Voronoi', 'tessellation']
    close_caption = ['2  Thick cell walls', 'in latewood']
    site_r = np.linalg.norm(d2d.sites[:, :2] - d2d.center, axis=1)
    rmax = max(d2d.mid_r.max(), site_r.max())
    for f in range(12):
        front = rmax*smoothstep((f + 1)/12)
        image = rings.copy()
        d2d.nuclei(image, np.nonzero(site_r > front)[0], wide)
        d2d.walls(image, wide, max_r=front, tint=0)
        emit(image, caption)
    full = d2d.full_view()
    walls_full = d2d.canvas()
    d2d.walls(walls_full, full)

    # zoom target: the outer edge of a dense (latewood) ring, lower right of the pith
    ring = d2d.radii[2*(len(d2d.radii)//4) + 1]
    target = d2d.center + ring*np.array([math.cos(math.radians(-25)), math.sin(math.radians(-25))])
    s_close = 14*scale
    def zoom_view(u, s_far):
        """u = 0: the pith at CENTER_PX at s_far px/mm; u = 1: the target at the canvas centre, s_close"""
        s = s_far*(s_close/s_far)**u
        w = (1 - s_far/s)/(1 - s_far/s_close)
        c = d2d.center + (target - d2d.center)*w
        px = (CENTER_PX[0] + (W/2 - CENTER_PX[0])*w, CENTER_PX[1] + (H/2 - CENTER_PX[1])*w)
        return (tuple(c), s, px)
    for f in range(10):
        image = d2d.canvas()
        d2d.walls(image, zoom_view(ease_in_out_cubic((f + 1)/10), wide[1]), tint=0)
        emit(image, caption if f < 7 else close_caption, fills_canvas=True)
    close = d2d.canvas()
    d2d.walls(close, zoom_view(1.0, scale))
    for f in range(18):
        emit(close, close_caption, fills_canvas=True)
    for f in range(7):
        image = d2d.canvas()
        d2d.walls(image, zoom_view(1 - ease_in_out_cubic((f + 1)/7), scale), tint=smoothstep((f + 1)/7))
        emit(image, close_caption, fills_canvas=True)
    # Clipping: the log is clipped with several polygons in turn, each time RingsPy's clipped mesh is shown
    # over the dimmed log; the last one is the hexagon of hero_log.toml, whose outside then fades out.
    # (the walls of the whole log reach the caption corner, so these frames get the caption backdrop)
    caption = ['3  Clipped to', 'any polygon']
    clips = np.load(data / 'clips.npz')
    dimmed = blend(walls_full, d2d.canvas(), 0.8)
    shown = walls_full                         # what the outline of the next polygon is traced over
    n_trace, n_reveal, n_hold = 5, 3, 7
    for key in clips['order']:
        polygon = clips['polygon_' + key]
        clipped = dimmed.copy()
        d2d.mesh(clipped, clips['vertices_' + key], clips['ridges_' + key], full)
        d2d.outline(clipped, full, len(polygon), polygon)
        for f in range(n_trace):                # the previous clip fades while the new outline is traced
            image = blend(shown, dimmed, (f + 1)/n_trace)
            d2d.outline(image, full, int(math.ceil(len(polygon)*(f + 1)/n_trace)), polygon)
            emit(image, caption, fills_canvas=True)
        traced = blend(shown, dimmed, 1.0)
        d2d.outline(traced, full, len(polygon), polygon)
        for f in range(n_reveal):
            emit(blend(traced, clipped, (f + 1)/n_reveal), caption, fills_canvas=True)
        for f in range(n_hold):
            emit(clipped, caption, fills_canvas=True)
        shown = clipped
    inside = d2d.canvas()
    d2d.mesh(inside, clips['vertices_box'], clips['ridges_box'], full)
    d2d.outline(inside, full, len(clips['polygon_box']), clips['polygon_box'])
    for f in range(4):
        emit(blend(shown, inside, (f + 1)/4), caption, fills_canvas=f < 3)
    slab = r3d.render(r3d.layer_actor(r3d.levels[0]), 90.0, -90.0, 0.0)
    for f in range(4):
        emit(blend(inside, slab, (f + 1)/4), caption)
    for f in range(2):
        emit(slab, caption)
    print('S3-S4 done ({} frames, {:.0f} s)'.format(len(frames), time.time() - t))

    # S5: the camera tilts to the S1 elevation over the bottom slab, the layers are added at a fixed
    # azimuth (so all the turning seen is RingsPy's twist), then the camera eases into the S1 orbit:
    # the last frame is one orbit step before frame 0
    caption = ['4  Extruded with', 'a twist']
    n_tilt, n_build, n_ease = 9, 17, 6
    azimuth = np.zeros(n_ease)
    azimuth[-1] = az0 - step
    for j in range(n_ease - 1, 0, -1):          # speed ramps up to `step` per frame
        azimuth[j - 1] = azimuth[j] - step*(j + 1)/n_ease
    az_build = azimuth[0] - step/n_ease
    # the build steps through the z levels that end in a flat top (the lower halves of the thin
    # longitudinal connectors between segments end in a sheared plane, which would flash dark)
    build_levels = r3d.levels[np.append(np.diff(r3d.levels) > 0.02, True)]
    for f in range(n_tilt):
        s = smoothstep((f + 1)/n_tilt)
        image = r3d.render(r3d.layer_actor(r3d.levels[0]), 90 - (90 - el_hero)*s, -90 + (az_build + 90)*s,
                           r3d.zmid*s)
        emit(image, caption)
    for f in range(n_build):
        level = build_levels[int(round((len(build_levels) - 1)*ease_in_out_cubic((f + 1)/n_build)))]
        emit(r3d.render(r3d.layer_actor(level), el_hero, az_build, r3d.zmid), caption)
    for f in range(n_ease):
        emit(r3d.render(r3d.solid_actor, el_hero, azimuth[f], r3d.zmid), caption)
    print('S5 done ({} frames, {:.0f} s)'.format(len(frames), time.time() - t))

    for i, (small, caption, fills_canvas, legend) in enumerate(frames):
        annotate(small, caption, fonts, backdrop=fills_canvas, legend=legend).save(
            frames_dir / 'frame_{:04d}.png'.format(i))
    if overlaps:
        raise RuntimeError('the object overlaps the captions or the legend in frames (index, pixels) {}; '
                           'move CENTER_PX or reduce the scale'.format(overlaps))


def load_fonts():
    from matplotlib import font_manager
    from PIL import ImageFont
    bold = font_manager.findfont(font_manager.FontProperties(family='DejaVu Sans', weight='bold'))
    regular = font_manager.findfont(font_manager.FontProperties(family='DejaVu Sans'))
    return dict(wordmark=ImageFont.truetype(bold, round(30*K)), caption=ImageFont.truetype(regular, round(26*K)),
                label=ImageFont.truetype(bold, round(25*K)), legend=ImageFont.truetype(regular, round(20*K)))


def annotate(image, caption, fonts, backdrop=False, legend=False):
    """Wordmark and caption lines at the top left (the last line of a 3-line caption is the file type,
    in bold), the colour key at the bottom left; each group on a darkened backdrop when the object
    fills the canvas. All text must fit in CAPTION_BOX / LEGEND_BOX."""
    from PIL import Image, ImageDraw
    image = image.copy()
    draw = ImageDraw.Draw(image)
    x = round(22*K)
    texts = [((x, round(14*K)), 'RingsPy', fonts['wordmark'], LIGHT)]
    for i, line in enumerate(caption):
        bold = len(caption) == 3 and i == 2
        texts.append(((x, round((56 + 32*i + (6 if bold else 0))*K)), line,
                      fonts['label' if bold else 'caption'], LIGHT if bold else CAPTION))
    swatches = []
    if legend:
        for i, (color, text) in enumerate([(LATE, 'latewood: thick walls'), (EARLY, 'earlywood: thin walls')]):
            y = H - round((60 - 28*i)*K)
            texts.append(((x + round(24*K), y), text, fonts['legend'], CAPTION))
            swatches.append(((x, y + round(4*K), x + round(15*K), y + round(19*K)), color))
    boxes = np.array([draw.textbbox(xy, text, font=font) for xy, text, font, _ in texts] +
                     [box for box, _ in swatches]).reshape(-1, 4)
    lower = boxes[:, 1] > H/2
    for box, is_lower in zip(boxes, lower):
        limits = LEGEND_BOX if is_lower else CAPTION_BOX
        assert box[0] >= limits[0] and box[2] <= limits[2] and box[1] >= limits[1] and box[3] <= limits[3], \
            'text outside its box: {}'.format(box)
    if backdrop:
        for group in (boxes[~lower], boxes[lower]):
            if len(group):
                x0, y0 = group[:, :2].min(axis=0) - 10
                x1, y1 = group[:, 2:].max(axis=0) + 10
                region = image.crop((x0, y0, x1, y1))
                image.paste(Image.blend(region, Image.new('RGB', region.size, BG), 0.85), (int(x0), int(y0)))
    for box, color in swatches:
        draw.rectangle(box, fill=color)
    for xy, text, font, color in texts:
        draw.text(xy, text, font=font, fill=color)
    return image


# ============================================================================ GIF encoding
def shade_palette(frames, ncolors):
    """The flat colours, their shades under the flat-shaded lighting, and their mixes with the
    background (sub-pixel walls, fades), completed by a median cut of the frames' other colours."""
    def ramp(color, factors):
        return [tuple(np.clip(np.rint(np.array(color)*f), 0, 255).astype(int)) for f in factors]
    def mixes(color, ts):
        return [tuple(np.rint(np.array(BG) + (np.array(color) - np.array(BG))*t).astype(int)) for t in ts]
    entries = [BG, GUIDE, LIGHT, CAPTION]
    entries += ramp(EARLY, (.45, .55, .65, .75, .88, 1)) + ramp(LATE, (.45, .6, .75, .88, 1))
    entries += ramp(IVORY, (.33, .37, .41, .45, .55, .62, .68, .75, .88, 1)) + ramp(BEAM, (.8, 1))
    for color in (EARLY, LATE, IVORY, LIGHT):
        entries += mixes(color, (.25, .5, .75))
    entries += mixes(LATE, (.12,))
    entries = list(dict.fromkeys(entries))
    # the rest from the frames' pixels that are not close to these yet
    from PIL import Image
    sample = np.concatenate([np.asarray(f).reshape(-1, 3)[::97] for f in frames[::3]])
    known = np.array(entries)
    far = np.min(np.abs(sample[:, None, :].astype(int) - known[None]).sum(axis=2), axis=1) > 24
    extra = sample[far]
    if len(extra) and len(entries) < ncolors:
        strip = Image.fromarray(extra.reshape(1, -1, 3).astype(np.uint8))
        q = strip.quantize(colors=ncolors - len(entries), method=Image.Quantize.MEDIANCUT, dither=Image.Dither.NONE)
        entries += [tuple(c) for c in np.array(q.getpalette()[:3*(ncolors - len(entries))]).reshape(-1, 3)]
    return entries[:ncolors]


def encode(frames, path, preview=None, ncolors=48, budget_bytes=5_000_000):
    """GIF with one global palette and no dithering. Identical consecutive frames are merged; each
    stored frame keeps only its changed pixels (the others are transparent) when that is smaller."""
    from PIL import Image, ImageSequence

    entries = shade_palette(frames, ncolors)
    transparent = len(entries)
    flat = np.zeros((256, 3), dtype=np.uint8)
    flat[:len(entries)] = entries
    palette_image = Image.new('P', (1, 1))
    palette_image.putpalette(flat.ravel().tolist())
    indices = [np.asarray(f.quantize(palette=palette_image, dither=Image.Dither.NONE)) for f in frames]

    def to_image(a):
        image = Image.fromarray(a.astype(np.uint8), mode='P')
        image.putpalette(flat.ravel().tolist())
        return image

    def gif_size(a):
        buf = io.BytesIO()
        to_image(a).save(buf, format='GIF', transparency=transparent, optimize=False)
        return buf.tell()

    stored, durations, full = [], [], []
    for a in indices:
        if full and np.array_equal(a, full[-1]):
            durations[-1] += DELAY_MS
            continue
        if full:
            changed = a != full[-1]
            rows, cols = np.nonzero(changed.any(axis=1))[0], np.nonzero(changed.any(axis=0))[0]
            box = (slice(rows[0], rows[-1] + 1), slice(cols[0], cols[-1] + 1))
            holes = np.where(changed, a, transparent)
            frame = holes if gif_size(holes[box]) < gif_size(a[box]) else a
        else:
            frame = a
        stored.append(to_image(frame))
        full.append(a)
        durations.append(DELAY_MS)
    stored[0].save(path, save_all=True, append_images=stored[1:], duration=durations, loop=0, disposal=1,
                   transparency=transparent, optimize=False)

    # check the file: every decoded frame equals its quantized source, the durations add up
    size = path.stat().st_size
    with Image.open(path) as gif:
        decoded = [(np.asarray(f.convert('RGB')), f.info['duration']) for f in ImageSequence.Iterator(gif)]
    assert len(decoded) == len(full), (len(decoded), len(full))
    for (pixels, duration), a, d in zip(decoded, full, durations):
        assert np.array_equal(pixels, flat[a]) and duration == d
    assert sum(durations) == DELAY_MS*len(frames)
    print('{}: {} frames ({} stored), {:.1f} s, {:,} bytes'.format(path.name, len(frames), len(full),
                                                                  sum(durations)/1000, size))
    if preview:
        frames[0].save(preview)
    if size > budget_bytes:
        raise RuntimeError('{:,} bytes is over the budget of {:,}'.format(size, budget_bytes))
    return size


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    parser.add_argument('step', choices=['prepare', 'render', 'encode'])
    parser.add_argument('--data', type=Path, required=True, help='folder for the compact data and frames')
    parser.add_argument('--work', type=Path, help='folder for the full mesh files (prepare; default: a temporary folder)')
    parser.add_argument('--gif', type=Path, default=GIF)
    args = parser.parse_args()
    if args.step == 'prepare':
        if args.work:
            prepare(args.data, args.work)
        else:
            with tempfile.TemporaryDirectory() as work:
                prepare(args.data, Path(work))
    else:
        render(args.data, args.gif, encode_only=(args.step == 'encode'))


if __name__ == '__main__':
    main()
