# -*- coding: utf-8 -*-
"""
Rebuild the README screenshots contents/Model2DPreview.png (previews of the 2D cross-sections of seven
examples, drawn with result.plot()) and contents/ModelVisualization.png (three off-screen VTK renders in
the style of ParaView: the components of wood_cube coloured by Connector_width, the connector volumes of
the hardwood example with the Vessel and Ray cell data, and the tapered, spiral-grained volumes of
wood_spiral coloured by Ring).

    python contents/make_screenshots.py --work <folder>

The meshes of wood_cube, hardwood_ring_porous and wood_spiral are generated in <folder> (a few minutes),
the previews are drawn from the examples, and the two images are written next to this script. Needs
vtk and Pillow, and an OpenGL context for the off-screen renders (OSMesa or Xvfb will do).
"""

import argparse
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
EXAMPLES = HERE.parent / 'examples'
OUT = HERE
PREVIEWS = [('wood_cube', 'ring'), ('wood_cube_lloyd', 'ring'), ('hardwood_ring_porous', 'type'),
            ('wood_eccentric', 'dense'), ('wood_spiral', 'ring'), ('honeycomb_hex', 'ring'), ('dogbone_svg', 'thickness')]
MESHES = ['wood_cube', 'hardwood_ring_porous', 'wood_spiral']


def previews_2d(work):
    """Model2DPreview.png: six square panels and the wide dog-bone panel, each the preview of an example"""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib import gridspec
    import ringspy
    fig = plt.figure(figsize=(19, 17))
    gs = gridspec.GridSpec(3, 3, figure=fig, height_ratios=[1, 1, 0.75])
    axes = [fig.add_subplot(gs[i // 3, i % 3]) for i in range(6)] + [fig.add_subplot(gs[2, :])]
    for ax, (name, color) in zip(axes, PREVIEWS):
        r = ringspy.preview(EXAMPLES / (name + '.toml'), path=str(work / 'previews'))
        r.plot(ax=ax, color=color, legend=(name != 'honeycomb_hex'))
        ax.set_title(ax.get_title().replace(name + ': ', ''), fontsize=10)
        ax.text(0.5, 1.06, '%s.toml   (color="%s")' % (name, color), transform=ax.transAxes, ha='center',
                fontsize=13, fontweight='bold')
    fig.suptitle('Previews of 2D cross-sections of the examples: ringspy run <example>.toml --preview, '
                 'drawn with result.plot()', fontsize=15, y=0.995)
    fig.tight_layout(rect=(0, 0, 1, 0.975))
    fig.savefig(OUT / 'Model2DPreview.png', dpi=72)
    plt.close(fig)


def views_3d(work):
    """ModelVisualization.png: the three VTK renders side by side"""
    import vtk
    from vtkmodules.util.numpy_support import vtk_to_numpy, numpy_to_vtk
    from PIL import Image
    import ringspy
    for name in MESHES:
        if not (work / 'meshes' / name / (name + '_conns_vol.vtu')).is_file():
            ringspy.generate(EXAMPLES / (name + '.toml'), path=str(work / 'meshes'), inpFlag=False, stlFlag=False)
    m = str(work / 'meshes') + '/'
    import vtk
    from vtkmodules.util.numpy_support import vtk_to_numpy, numpy_to_vtk
    from PIL import Image
    W, H = 900, 760
    BG = (0.20, 0.22, 0.30)

    def read(path):
        r = vtk.vtkXMLUnstructuredGridReader(); r.SetFileName(path); r.Update(); return r.GetOutput()

    def coolwarm():
        lut = vtk.vtkColorTransferFunction()
        lut.AddRGBPoint(0.0, 0.23, 0.30, 0.75); lut.AddRGBPoint(0.5, 0.86, 0.86, 0.86); lut.AddRGBPoint(1.0, 0.70, 0.015, 0.15)
        return lut

    def actor(grid, array, lut=None, rng=None, line_width=1.0, rgb=None):
        m = vtk.vtkDataSetMapper(); m.SetInputData(grid)
        if rgb is not None:
            colors = numpy_to_vtk(rgb.astype(np.uint8), deep=True); colors.SetName('rgb'); grid.GetCellData().AddArray(colors)
            m.SetScalarModeToUseCellFieldData(); m.SelectColorArray('rgb'); m.SetColorModeToDirectScalars()
        elif array:
            m.SetScalarModeToUseCellFieldData(); m.SelectColorArray(array); m.SetLookupTable(lut); m.SetScalarRange(*rng); m.SetColorModeToMapScalars()
            m.InterpolateScalarsBeforeMappingOff()
        else:
            m.ScalarVisibilityOff()
        a = vtk.vtkActor(); a.SetMapper(m); p = a.GetProperty(); p.SetLineWidth(line_width); p.SetInterpolationToFlat(); p.SetAmbient(0.3); p.SetDiffuse(0.75); p.SetSpecular(0.0)
        return a

    def render(actors, label, bar=None, az=35, el=28, name='out'):
        ren = vtk.vtkRenderer(); ren.SetBackground(*BG)
        for a in actors: ren.AddActor(a)
        if bar is not None:
            lut, title, rng = bar
            sb = vtk.vtkScalarBarActor(); sb.SetLookupTable(lut); sb.SetTitle(title); sb.SetNumberOfLabels(5)
            sb.SetOrientationToVertical(); sb.SetWidth(0.1); sb.SetHeight(0.42); sb.SetPosition(0.085, 0.54)
            for tp in (sb.GetTitleTextProperty(), sb.GetLabelTextProperty()): tp.SetColor(1, 1, 1); tp.SetFontSize(26); tp.BoldOff(); tp.ItalicOff(); tp.ShadowOff()
            sb.SetUnconstrainedFontSize(True); ren.AddViewProp(sb)
        txt = vtk.vtkTextActor(); txt.SetInput(label); tp = txt.GetTextProperty(); tp.SetFontSize(27); tp.SetColor(1, 1, 1); txt.SetPosition(20, 20); ren.AddViewProp(txt)
        axes = vtk.vtkAxesActor(); axes.SetTotalLength(1, 1, 1)
        om = vtk.vtkOrientationMarkerWidget()
        rw = vtk.vtkRenderWindow(); rw.SetOffScreenRendering(1); rw.SetSize(W, H); rw.SetMultiSamples(8); rw.AddRenderer(ren)
        iren = vtk.vtkRenderWindowInteractor(); iren.SetRenderWindow(rw)
        om.SetOrientationMarker(axes); om.SetInteractor(iren); om.SetViewport(0.0, 0.0, 0.18, 0.18); om.SetEnabled(1); om.InteractiveOff()
        cam = ren.GetActiveCamera(); b = np.array(actors[0].GetBounds()).reshape(3, 2); c = b.mean(axis=1); d = 3.2*np.ptp(b, axis=1).max()
        cam.SetFocalPoint(*c); cam.SetPosition(*(c + d*np.array([0.95, -1.15, 0.8])/1.7)); cam.SetViewUp(0, 0, 1); cam.SetViewAngle(24)
        ren.ResetCameraClippingRange()
        light = vtk.vtkLight(); light.SetLightTypeToCameraLight(); light.SetPosition(1, 1, 1); light.SetIntensity(0.9); ren.AddLight(light)
        rw.Render()
        g = vtk.vtkWindowToImageFilter(); g.SetInput(rw); g.ReadFrontBufferOff(); g.Update()
        im = g.GetOutput(); w, h, _ = im.GetDimensions()
        px = vtk_to_numpy(im.GetPointData().GetScalars()).reshape(h, w, -1)[::-1, :, :3]
        Image.fromarray(np.ascontiguousarray(px)).save(str(work / ('panel_%s.png' % name)))

    # (a) wood_cube components: vertices, beams and connectors, connectors coloured by Connector_width
    conns = read(m + 'wood_cube/wood_cube_conns.vtu'); beams = read(m + 'wood_cube/wood_cube_beams.vtu'); verts = read(m + 'wood_cube/wood_cube_vertices.vtu')
    w = vtk_to_numpy(conns.GetCellData().GetArray('Connector_width')); rng = (float(w.min()), float(w.max()))
    lut = coolwarm(); lut.SetRange(*rng) if hasattr(lut, 'SetRange') else None
    lutc = vtk.vtkColorTransferFunction()
    for x, c in [(rng[0], (0.23, 0.30, 0.75)), (0.5*(rng[0]+rng[1]), (0.86, 0.86, 0.86)), (rng[1], (0.70, 0.015, 0.15))]: lutc.AddRGBPoint(x, *c)
    b = actor(beams, None); b.GetProperty().SetColor(0.55, 0.62, 0.85); b.GetProperty().SetLineWidth(1.2)
    render([actor(conns, 'Connector_width', lutc, rng, line_width=2.0), b], '(a) wood_cube: _beams.vtu + _conns.vtu, Connector_width', bar=(lutc, 'Connector_width', rng), name='a')
    # (b) hardwood connector volumes: Connector_width, vessels amber, rays green
    vol = read(m + 'hardwood_ring_porous/hardwood_ring_porous_conns_vol.vtu')
    w = vtk_to_numpy(vol.GetCellData().GetArray('Connector_width')); rng2 = (float(w.min()), float(w.max()))
    lut2 = vtk.vtkColorTransferFunction()
    for x, c in [(rng2[0], (0.23, 0.30, 0.75)), (0.5*(rng2[0]+rng2[1]), (0.86, 0.86, 0.86)), (rng2[1], (0.70, 0.015, 0.15))]: lut2.AddRGBPoint(x, *c)
    rgb = np.array([lut2.GetColor(x) for x in w])*255
    vessel = vtk_to_numpy(vol.GetCellData().GetArray('Vessel')) > 0; ray = vtk_to_numpy(vol.GetCellData().GetArray('Ray')) > 0
    rgb[vessel] = (242, 183, 60); rgb[ray] = (111, 208, 142)
    render([actor(vol, None, rgb=rgb)], '(b) hardwood_ring_porous: _conns_vol.vtu, Vessel + Ray', bar=(lut2, 'Connector_width', rng2), name='b')
    # (c) wood_spiral volumes: taper and spiral grain, coloured by Ring
    vol3 = read(m + 'wood_spiral/wood_spiral_conns_vol.vtu')
    ring = vtk_to_numpy(vol3.GetCellData().GetArray('Ring')); rng3 = (0.0, float(ring.max()))
    lut3 = vtk.vtkLookupTable(); lut3.SetNumberOfTableValues(int(ring.max()) + 1); lut3.SetRange(*rng3)
    import matplotlib; cmap = matplotlib.colormaps['viridis']
    for k in range(int(ring.max()) + 1): lut3.SetTableValue(k, *cmap(0.1 + 0.8*k/max(ring.max(), 1))[:3], 1.0)
    lut3.Build()
    render([actor(vol3, 'Ring', lut3, rng3)], '(c) wood_spiral: _conns_vol.vtu, Ring (taper, spiral grain)', bar=(lut3, 'Ring', rng3), az=30, el=22, name='c')

    ims = [Image.open(str(work / ('panel_%s.png' % k))) for k in 'abc']
    sheet = Image.new('RGB', (W*3 + 20, H), (51, 56, 77))
    for i, im in enumerate(ims): sheet.paste(im, (i*(W + 10), 0))
    sheet = sheet.resize((1818, round(1818*H/(W*3 + 20))), Image.LANCZOS); sheet.save(str(OUT / 'ModelVisualization.png')); print(sheet.size)

def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    parser.add_argument('--work', type=Path, required=True, help='folder for the meshes and previews')
    args = parser.parse_args()
    args.work.mkdir(parents=True, exist_ok=True)
    previews_2d(args.work)
    views_3d(args.work)
    print('written', OUT / 'Model2DPreview.png', OUT / 'ModelVisualization.png')


if __name__ == '__main__':
    main()
