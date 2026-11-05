# -*- coding: utf-8 -*-
"""
Clipping outlines from SVG and DXF files (ringspy.outline), and box_polygon naming a file
"""

import numpy as np
import pytest

import ringspy
from ringspy import outline
from ringspy.cli import main
from helpers import EXAMPLES_DIR, assert_same_files

DOG = [[-1.2, -0.8], [-0.5, -0.8], [-0.1, -0.4], [1.3, -0.4], [1.7, -0.8], [2.4, -0.8],
       [2.4, 0.8], [1.7, 0.8], [1.3, 0.4], [-0.1, 0.4], [-0.5, 0.8], [-1.2, 0.8]]


def dxf_text(entries):
    """A minimal DXF file: the header's $INSUNITS and the given entity (code, value) pairs"""
    lines = ['0', 'SECTION', '2', 'HEADER', '9', '$INSUNITS', '70', '4', '0', 'ENDSEC', '0', 'SECTION', '2', 'ENTITIES']
    for code, value in entries:
        lines += [str(code), str(value)]
    lines += ['0', 'ENDSEC', '0', 'EOF']
    return '\n'.join(lines) + '\n'


def same_polygon(a, b, tol=1e-9):
    """Two closed polygons with the same vertices, up to a cyclic shift and the direction"""
    a, b = np.asarray(a, dtype=float), np.asarray(b, dtype=float)
    if len(a) != len(b):
        return False
    for direction in (b, b[::-1]):
        for shift in range(len(b)):
            if np.allclose(a, np.roll(direction, shift, axis=0), atol=tol):
                return True
    return False


# ------------------------------------------------------------------------------------------ SVG

def test_svg_path_commands_and_transforms(tmp_path):
    svg = tmp_path / "shape.svg"
    svg.write_text('<svg xmlns="http://www.w3.org/2000/svg" width="40mm" height="20mm" viewBox="0 0 40 20">'
                   '<g transform="translate(5,5)"><path d="M 0 0 l 30 0 v 10 H 0 z"/></g>'
                   '<rect x="1" y="1" width="2" height="2"/></svg>')
    p = outline.read_outline(svg)
    # the largest shape, translated, in mm, with y flipped (y down in SVG)
    assert same_polygon(p, [[5, -5], [35, -5], [35, -15], [5, -15]])
    # curves are flattened, relative commands and shorthands work
    svg.write_text('<svg xmlns="http://www.w3.org/2000/svg"><path d="M0,0 C 10,-5 20,-5 30,0 S 50,-5 60,0 Q 30,20 0,0 Z"/></svg>')
    p = outline.read_outline(svg)
    assert len(p) > 40 and np.ptp(p[:, 0]) == pytest.approx(60*25.4/96, rel=1e-6)  # px at 96 dpi without units
    # arcs: a half disc of radius 10
    svg.write_text('<svg xmlns="http://www.w3.org/2000/svg" width="100mm" viewBox="0 0 100 100">'
                   '<path d="M 0 0 A 10 10 0 0 1 20 0 Z"/></svg>')
    p = outline.read_outline(svg)
    area = 0.5*abs(np.sum(p[:, 0]*np.roll(p[:, 1], -1) - np.roll(p[:, 0], -1)*p[:, 1]))
    assert area == pytest.approx(np.pi*10**2/2, rel=0.01)
    # rotate and scale transforms, a polygon element
    svg.write_text('<svg xmlns="http://www.w3.org/2000/svg" width="10mm" viewBox="0 0 10 10">'
                   '<polygon transform="rotate(90) scale(2)" points="0,0 1,0 1,1 0,1"/></svg>')
    p = outline.read_outline(svg)
    assert np.ptp(p, axis=0) == pytest.approx([2, 2])
    with pytest.raises(ValueError, match="no closed"):
        svg.write_text('<svg xmlns="http://www.w3.org/2000/svg"><line x1="0" y1="0" x2="1" y2="1"/></svg>')
        outline.read_outline(svg)


def test_svg_subpaths_and_shapes():
    loops = outline.svg_path_subpaths("M 0 0 L 1 0 L 1 1 Z M 5 5 h 3 v 3 h -3 z")
    assert len(loops) == 2 and len(loops[0]) == 3 and len(loops[1]) == 4
    with pytest.raises(ValueError, match="malformed"):
        outline.svg_path_subpaths("M 1 2 L 3")


# ------------------------------------------------------------------------------------------ DXF

def test_dxf_polyline_lines_arcs_and_units(tmp_path):
    dxf = tmp_path / "shape.dxf"
    # a closed LWPOLYLINE with a bulge (a half disc on its right side), and a loop of lines (smaller)
    dxf.write_text(dxf_text([(0, 'LWPOLYLINE'), (70, 1), (10, 0), (20, 0), (42, 0), (10, 4), (20, 0), (42, 1),
                             (10, 4), (20, 2), (42, 0), (10, 0), (20, 2), (42, 0),
                             (0, 'LINE'), (10, 10), (20, 10), (11, 12), (21, 10),
                             (0, 'LINE'), (10, 12), (20, 10), (11, 12), (21, 11),
                             (0, 'LINE'), (10, 12), (20, 11), (11, 10), (21, 10)]))
    p = outline.read_outline(dxf)
    assert np.ptp(p, axis=0) == pytest.approx([5, 2], abs=1e-6)
    area = 0.5*abs(np.sum(p[:, 0]*np.roll(p[:, 1], -1) - np.roll(p[:, 0], -1)*p[:, 1]))
    assert area == pytest.approx(8 + np.pi/2, rel=0.01)
    # a loop of lines and an arc, in inches
    text = dxf_text([(0, 'LINE'), (10, 0), (20, 0), (11, 2), (21, 0), (0, 'ARC'), (10, 2), (20, 1), (40, 1), (50, 270), (51, 90),
                     (0, 'LINE'), (10, 2), (20, 2), (11, 0), (21, 2), (0, 'LINE'), (10, 0), (20, 2), (11, 0), (21, 0)])
    dxf.write_text(text.replace('70\n4\n', '70\n1\n'))
    p = outline.read_outline(dxf)
    assert np.ptp(p, axis=0) == pytest.approx([3*25.4, 2*25.4], rel=1e-6)
    # a circle, and the polyline chosen over it when larger
    dxf.write_text(dxf_text([(0, 'CIRCLE'), (10, 0), (20, 0), (40, 1)]))
    p = outline.read_outline(dxf)
    assert len(p) >= 16 and np.allclose(np.linalg.norm(p, axis=1), 1.0)
    with pytest.raises(ValueError, match="no closed"):
        dxf.write_text(dxf_text([(0, 'LINE'), (10, 0), (20, 0), (11, 1), (21, 1)]))
        outline.read_outline(dxf)
    with pytest.raises(ValueError, match=".svg or .dxf"):
        outline.read_outline(tmp_path / "shape.stl")


# --------------------------------------------------------------------------------- simplification

def test_simplify_polygon():
    t = np.linspace(0, 2*np.pi, 721, endpoint=False)
    circle = np.column_stack((np.cos(t), np.sin(t)))
    for tol in [0.01, 0.001]:
        p = outline.simplify_polygon(circle, tol)
        assert 3 < len(p) < len(circle)
        # within the tolerance of the circle: every vertex on it, and the chords' sagitta below tol
        np.testing.assert_allclose(np.linalg.norm(p, axis=1), 1.0, atol=1e-12)
        chord = np.linalg.norm(np.roll(p, -1, axis=0) - p, axis=1)
        assert np.all(1 - np.sqrt(1 - (chord/2)**2) <= tol*1.01)
    assert len(outline.simplify_polygon(circle, 0.001)) > len(outline.simplify_polygon(circle, 0.01))
    # a polygon's corners are kept
    square = np.array([[0, 0], [0.5, 0], [1, 0], [1, 1], [0, 1], [0, 0.5]], dtype=float)
    assert same_polygon(outline.simplify_polygon(square, 0.01), [[0, 0], [1, 0], [1, 1], [0, 1]])
    assert len(outline.simplify_polygon(square, 0)) == 6


# ------------------------------------------------------------------------------------------ Config

def test_config_reads_an_outline_file(tmp_path):
    svg = tmp_path / "dog.svg"
    svg.write_text('<svg xmlns="http://www.w3.org/2000/svg" width="20mm" height="10mm" viewBox="0 0 20 10"><path d="M '
                   + ' L '.join('{} {}'.format(x + 10, -y + 5) for x, y in DOG) + ' Z"/></svg>')
    c = ringspy.Config(box_shape="polygon", box_polygon=str(svg), box_center=[0.6, 0.0], box_polygon_tolerance=0.01)
    assert isinstance(c.box_polygon, str)
    vertices = c.polygon_vertices()
    assert same_polygon(vertices, DOG, tol=1e-9)
    r = c.resolved()
    assert isinstance(r.box_polygon, list) and same_polygon(r.box_polygon, DOG) and r.z_max == pytest.approx(3.6)
    # from a TOML file, relative to the file's folder
    toml = tmp_path / "dog.toml"
    toml.write_text('box_shape = "polygon"\nbox_polygon = "dog.svg"\nbox_center = [0.6, 0.0]\n')
    c2 = ringspy.load_config(toml)
    assert c2.box_polygon == str(svg) and same_polygon(c2.polygon_vertices(), DOG)
    with pytest.raises(ValueError, match="does not exist"):
        ringspy.Config(box_shape="polygon", box_polygon=str(tmp_path / "missing.svg"))
    with pytest.raises(ValueError, match=".svg or .dxf"):
        ringspy.Config(box_shape="polygon", box_polygon=str(toml))
    with pytest.raises(ValueError, match="box_polygon_tolerance"):
        ringspy.Config(box_shape="polygon", box_polygon=str(svg), box_polygon_tolerance=-1)
    with pytest.raises(ValueError, match="box_polygon_scale"):
        ringspy.Config(box_shape="polygon", box_polygon=str(svg), box_polygon_scale=0)


def test_svg_dogbone_gives_the_same_mesh_as_the_vertex_list(tmp_path):
    # examples/dogbone.svg holds the vertices of dogbone.toml: the same files, the log naming the file
    a = ringspy.generate(EXAMPLES_DIR / "dogbone.toml", geoName="dog", path=str(tmp_path / "a"), stlFlag=False)
    b = ringspy.generate(EXAMPLES_DIR / "dogbone_svg.toml", geoName="dog", path=str(tmp_path / "b"), stlFlag=False)
    assert same_polygon(b.config.box_polygon, a.config.box_polygon)
    log_b = (b.output_dir / "dog-report.log").read_text()
    assert "box_polygon:" in log_b and "dogbone.svg (12 vertices)" in log_b
    for name in ["dog-vertex.mesh", "dog-ridge.mesh", "dog-mesh.txt", "dogIGA.txt", "dog.inp"]:
        assert (a.output_dir / name).read_bytes() == (b.output_dir / name).read_bytes(), name


def test_cli_outline(tmp_path, capsys):
    dxf = tmp_path / "shape.dxf"
    dxf.write_text(dxf_text([(0, 'LWPOLYLINE'), (70, 1), (10, 0), (20, 0), (10, 4), (20, 0), (10, 4), (20, 2), (10, 0), (20, 2)]))
    assert main(["outline", str(dxf), "--center", "0", "0"]) == 0
    out = capsys.readouterr().out
    assert "4 vertices" in out and out.strip().splitlines()[-1].startswith("box_polygon = [[-2.0, -1.0]")
    with pytest.raises(SystemExit):
        main(["outline", str(tmp_path / "none.svg")])
